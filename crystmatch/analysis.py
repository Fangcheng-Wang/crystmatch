"""
Analyze and visualize CSMs.
"""

from .utils import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "pgf.texsystem": "xelatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.dpi': 150,
})

np.set_printoptions(suppress=True)
Cryst = Tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]

def multiplicity(crystA: Cryst, crystB: Cryst, slist: Union[SLM, List[SLM], NDArray[np.int32]]) -> Union[int, NDArray[np.int32]]:
    """Return multiplicities of elements in `slist`.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slist : {(3, 3, 3), (..., 3, 3, 3)} array_like of ints
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.

    Returns
    -------
    mu : int or (...,) array of ints
        Multiplicities of each SLM in `slist`.
    """
    zA = crystA[2].shape[0]
    zB = crystB[2].shape[0]
    dA = np.lcm(zA,zB) // zA
    if len(slist.shape) == 3:
        return la.det(slist[0]).round().astype(int) // dA
    else:
        return la.det(slist[:,0,:,:]).round().astype(int) // dA

def sing_val(crystA: Cryst, crystB: Cryst, slist: Union[List[SLM], NDArray[np.int32]]) -> NDArray[np.float64]:
    """Return singular values of elements in `slist`.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slist : list of slm
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.

    Returns
    -------
    sv : (..., 3) array
        Contains singular values of each SLM in `slist`.
    """
    cA = crystA[0].T
    cB = crystB[0].T
    hA, hB, q = zip(*slist)
    s_mat = cB @ np.array(hB) @ np.array(q) @ la.inv(cA @ np.array(hA))
    sv = la.svd(s_mat, compute_uv=False)
    return sv

def distance_slm(slist: ArrayLike, s0: ArrayLike, crystA: Cryst, crystB: Cryst) -> NDArray[np.float64]:
    """The Frobenius distance between SLMs.

    Parameters
    ----------
    slist : list of slm
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.
    s0 : slm
        `(hA, hB, q)`, representing a SLM.
    crystA, crystB : cryst
        `(lattice, species, positions)`, representing the crystal structure, usually obtained by `load_poscar`.
    
    Returns
    -------
    dlist : (...,) array
        Contains Frobenius distances from `slist` to `s0`, where equivalent SLMs coincide.
    """
    cA = crystA[0].T
    cB = crystB[0].T
    gA = get_pure_rotation(crystA)
    gB = get_pure_rotation(crystB)
    hA0, hB0, q0 = s0
    x0 = np.transpose(np.dot((gB @ hB0) @ q0, la.inv(gA @ hA0)), axes=[2,0,1,3]).reshape(-1,9)
    _, i = np.unique(x0.round(decimals=4), axis=0, return_index=True)
    cl0 = cB @ x0[i,:].reshape(-1,3,3) @ la.inv(cA)
    ss = cB @ slist[:,1,:,:] @ slist[:,2,:,:] @ la.inv(cA @ slist[:,0,:,:])
    dlist = np.amin(la.norm(ss.reshape(-1,1,9), cl0.reshape(1,-1,9), axis=2),axis=1)
    return dlist

def distance_shuffle(crystA: Cryst, j1: Tuple[Cryst, Cryst], j2: Tuple[Cryst, Cryst]) -> float:
    """The RMSD between shuffles.
    
    Parameters
    ----------
    crystA : cryst
        `(lattice, species, positions)`, representing the crystal structure, usually obtained by `load_poscar`.
    j1, j2 : 2-tuple
        `(crystA_sup, crystB_sup)`, representing a CSM.
    
    Returns
    -------
    d : float
        The RMSD between `j10` and `j20`, the shuffles obtained by omitting the lattice deformation in `j1` and `j2`.
    """
    # cell geometry
    cA = crystA[0].T
    cA_sup_1 = j1[0][0].T
    cA_sup_2 = j2[0][0].T
    mA1 = (la.inv(cA) @ cA_sup_1).round().astype(int)
    mA2 = (la.inv(cA) @ cA_sup_2).round().astype(int)
    mA_hyp = matrix_lcm(mA1, mA2)
    nA1 = (la.inv(mA1) @ mA_hyp).round().astype(int)
    nA2 = (la.inv(mA2) @ mA_hyp).round().astype(int)
    # fractional coordinates in hypercells
    pA_sup_1 = j1[0][2].T
    pB_sup_1 = j1[1][2].T
    pA_hyp_1 = la.inv(nA1) @ (pA_sup_1.reshape(3,1,-1) + int_vectors_inside(nA1).reshape(3,-1,1)).reshape(3,-1)
    pB_hyp_1 = la.inv(nA1) @ (pB_sup_1.reshape(3,1,-1) + int_vectors_inside(nA1).reshape(3,-1,1)).reshape(3,-1)
    pA_sup_2 = j2[0][2].T
    pB_sup_2 = j2[1][2].T
    pA_hyp_2 = la.inv(nA2) @ (pA_sup_2.reshape(3,1,-1) + int_vectors_inside(nA2).reshape(3,-1,1)).reshape(3,-1)
    pB_hyp_2 = la.inv(nA2) @ (pB_sup_2.reshape(3,1,-1) + int_vectors_inside(nA2).reshape(3,-1,1)).reshape(3,-1)
    # shift each coordinate of `pA_hyp_2` into [0,1)
    pA_hyp_2 = pA_hyp_2.round(decimals=7)
    pB_hyp_2 -= np.floor(pA_hyp_2)
    pA_hyp_2 -= np.floor(pA_hyp_2)
    pA_hyp_2, ind = np.unique(pA_hyp_2.round(decimals=7), axis=1, return_index=True)
    pB_hyp_2 = pB_hyp_2[:,ind]
    # try every translation element in L_A
    tA = la.inv(mA_hyp) @ int_vectors_inside(matrix_gcd(mA1, mA2))
    l1d_list = []
    for i in range(tA.shape[1]):
        t = tA[:,i]
        pA_hyp_1_t = (pA_hyp_1 + t.reshape(3,1)).round(decimals=7)
        pB_hyp_1_t = pB_hyp_1 + t.reshape(3,1)
        pB_hyp_1_t -= np.floor(pA_hyp_1_t)
        pA_hyp_1_t -= np.floor(pA_hyp_1_t)
        pA_hyp_1_t, ind = np.unique(pA_hyp_1_t.round(decimals=7), axis=1, return_index=True)
        pB_hyp_1_t = pB_hyp_1_t[:,ind]
        if not ((pA_hyp_1_t - pA_hyp_2).round(decimals=7) == 0).all():
            print(pA_hyp_1_t.mean(axis=1), pA_hyp_2.mean(axis=1))
            assert False
        # print(la.norm(cA @ mA_hyp @ (pB_hyp_1_t - pB_hyp_2), axis=0))
        l1d_list.append(la.norm(cA @ mA_hyp @ (pB_hyp_1_t - pB_hyp_2), axis=0).mean())
        '''
        for j in range(pA_hyp_1.shape[1]):
            diff = (pA_hyp_1[:,j] + t).reshape(-1,1) - pA_hyp_2
            ind = np.nonzero((diff.round() == diff.round(decimals=7-1)).all(axis=0))[0]
            assert len(ind) == 1
            ind = ind[0]
            d += np.sum((cA @ mA_hyp @ (pB_hyp_1[:,j] + t - diff[:,ind].round() - pB_hyp_2[:,ind])) ** 2)
        rmsdlist.append((d / pA_hyp_1.shape[1]) ** 0.5)
        '''
    i0 = np.argmin(l1d_list)
    return l1d_list[i0]

def orientation_relation(vA1: ArrayLike, vB1: ArrayLike, vA2: ArrayLike, vB2: ArrayLike) -> NDArray[np.float64]:
    """Rotation matrix `r` such that `r @ vA1` parallel to `vB1` and `r @ vA2` parallel to `vB2`.

    Parameters
    ----------
    vA1, vB1 : (3,) array_like
        Vectors (in cartesian coordinates) such that `r @ vA1` parallel to `vB1`.
    vA2, vB2 : (3,) array_like
        The same as `vA1` and `vB1`.
        
    Returns
    -------
    r : (3, 3) array
        A rotation matrix representing the given orientation relationship, which rotates `vA1` to `vB1` and `vA2` to `vB2`.
    """
    b = np.array([vA1, np.cross(vA1, vA2), np.cross(vA1, np.cross(vA1, vA2))]).T
    c = np.array([vB1, np.cross(vB1, vB2), np.cross(vB1, np.cross(vB1, vB2))]).T
    b = b * la.norm(b, axis=0).reshape(1,3) ** -1
    c = c * la.norm(c, axis=0).reshape(1,3) ** -1
    r = c @ b.T
    return r

def compare_orientation(
    crystA: Cryst, crystB: Cryst, slist: Union[List[SLM], NDArray[np.int32]], r: NDArray[np.float64], fix_usp: bool = False
) -> NDArray[np.float64]:
    """Calculate how much each SLM in `slist` differ from a given orientation relationship.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slist : (..., 3, 3, 3) array_like of ints
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.
    r : (3, 3) array
        A rotation matrix representing the given orientation relationship.
    fix_usp : bool, optional
        Whether to fix the uniformed scaled plane. Default is False.

    Returns
    -------
    anglelist : (...,) array
        Contains rotation angles that measure the difference of each SLM and the given orientation.
    """
    assert la.det(r).round(decimals=4) == 1
    cA = crystA[0].T
    cB = crystB[0].T
    rA = cA @ get_pure_rotation(crystA) @ la.inv(cA)
    rB = cB @ get_pure_rotation(crystB) @ la.inv(cB)
    r_equiv = np.transpose(np.dot(np.dot(rB, r), rA), axes=(2,0,1,3)).reshape(-1,3,3)
    hA, hB, q = zip(*slist)
    s = cB @ np.array(hB) @ np.array(q) @ la.inv(cA @ np.array(hA))
    u, sigma, vT = la.svd(s)
    rS = u @ vT
    if fix_usp:
        eps = np.array([[[0,0,0],[0,0,-1],[0,1,0]],[[0,0,1],[0,0,0],[-1,0,0]],[[0,-1,0],[1,0,0],[0,0,0]]])
        v_cross = np.tensordot(vT[:,1,:], eps, axes=(1,0))
        s1 = sigma[:,0]
        s2 = sigma[:,1]
        s3 = sigma[:,2]
        theta = np.arctan(np.sqrt((s1**2 - s2**2) * (s2**2 - s3**2)) / (s1 * s3 + s2**2))
        rH = np.eye(3).reshape(1,3,3) + np.sin(theta).reshape(-1,1,1) * v_cross + (1 - np.cos(theta)).reshape(-1,1,1) * (v_cross @ v_cross)
        rS = np.array([rS @ rH, rS @ la.inv(rH)])       # rS.shape = (2, ..., 3, 3)
        anglelist = np.arccos(np.clip(0.5 * (-1 + np.amax(np.trace(np.dot(la.inv(r_equiv), rS), axis1=1, axis2=4), axis=(0,1))), -1, 1))
    else:
        anglelist = np.arccos(np.clip(0.5 * (-1 + np.amax(np.trace(np.dot(la.inv(r_equiv), rS), axis1=1, axis2=3), axis=0)), -1, 1))
    return anglelist.round(decimals=7)

def scatter_colored(
    kappalist: ArrayLike, rmsdlist: ArrayLike, colorlist: ArrayLike, cbarlabel: str = 'Benchmark', save: Union[str, None] = None
) -> None:
    """Scatter plot of the CSMs with colorbar.

    Parameters
    ----------
    kappalist : (N,) array_like
        The root-mean-square strain of each CSM.
    rmsdlist : (N,) array_like
        The root-mean-square distance of each CSM.
    colorlist : (N,) array_like
        Some quantity of each CSM, which is to be colored.
    cbarlabel : str
        The label of the colorbar. Default is `Benchmark`.
    save : str, optional
        The filename to save. Directly show the figure if not given.
    """
    kappalist = np.array(kappalist)
    rmsdlist = np.array(rmsdlist)
    colorlist = np.array(colorlist)
    plt.figure()
    ax = plt.subplot()
    ind = np.argsort(colorlist)[::-1]
    sc = plt.scatter(kappalist[ind], rmsdlist[ind], marker='d', c=colorlist[ind], cmap=plt.cm.get_cmap('rainbow'), s=20)
    ind0 = colorlist==0
    n0 = np.sum(ind0)
    if n0 >= 1:
        print(f"\nThere are {n0:.0f} CSMs (indices: {', '.join(np.nonzero(ind0)[0].astype(str).tolist())}) with {cbarlabel}=0, emphasized by pink stars in the plot.")
        plt.scatter(kappalist[ind0], rmsdlist[ind0], marker='*', color=(1.0,0.75,0.95), s=12)
    plt.xlabel('Root-mean-square strain $\kappa$', fontsize=15)
    plt.ylabel('RMSD $d$ (\AA)', fontsize=15)
    plt.xlim(0, np.amax(kappalist) * 1.05)
    plt.ylim(min(0, np.amin(rmsdlist) - 0.1), np.amax(rmsdlist) * 1.05)
    cbar = plt.colorbar(sc, aspect=40)
    cbar.set_label(cbarlabel, fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    if save == None:
        plt.show()
    else:
        plt.savefig(f"{save}")
        print(f"\nScatter plot saved as '{save}'")