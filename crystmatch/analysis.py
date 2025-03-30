"""
Analyze and visualize CSMs.
"""

from .utilities import *
from matplotlib import rcParams, colors
import matplotlib.pyplot as plt
rcParams.update({
    'font.family': 'serif',
    'pgf.rcfonts': False,
    'figure.dpi': 150
})

np.set_printoptions(suppress=True)

def multiplicity(crystA: Cryst, crystB: Cryst, slmlist: Union[SLM, List[SLM], NDArray[np.int32]]) -> Union[int, NDArray[np.int32]]:
    """Return multiplicities of elements in `slmlist`.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slmlist : list of slm
        A list of SLMs, each represented by a triplet of integer matrices like `(hA, hB, q)`.

    Returns
    -------
    mu : int or (...,) array of ints
        Multiplicities of each SLM in `slmlist`.
    """
    slmlist = np.array(slmlist)
    zA = crystA[2].shape[0]
    zB = crystB[2].shape[0]
    dA = np.lcm(zA,zB) // zA
    if len(slmlist.shape) == 3:
        return la.det(slmlist[0]).round().astype(int) // dA
    else:
        return la.det(slmlist[:,0,:,:]).round().astype(int) // dA

def sing_val(crystA: Cryst, crystB: Cryst, slmlist: Union[List[SLM], NDArray[np.int32]]) -> NDArray[np.float64]:
    """Return singular values of elements in `slmlist`.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slmlist : list of slm
        A list of SLMs, each represented by a triplet of integer matrices like `(hA, hB, q)`.

    Returns
    -------
    sv : (..., 3) array
        Contains singular values of each SLM in `slmlist`.
    """
    cA = crystA[0].T
    cB = crystB[0].T
    hA, hB, q = zip(*slmlist)
    deform = cB @ np.array(hB) @ np.array(q) @ la.inv(cA @ np.array(hA))
    sv = la.svd(deform, compute_uv=False)
    return sv

def deform_distance(slmlist: ArrayLike, s0: ArrayLike, crystA: Cryst, crystB: Cryst) -> NDArray[np.float64]:
    """The Frobenius distance between deformation gradients.

    Parameters
    ----------
    slmlist : list of slm
        A list of SLMs, each represented by a triplet of integer matrices like `(hA, hB, q)`.
    s0 : slm
        `(hA, hB, q)`, representing a SLM.
    crystA, crystB : cryst
        `(lattice, species, positions)`, representing the crystal structure, usually obtained by `load_poscar`.
    
    Returns
    -------
    dlist : (...,) array
        Contains Frobenius distances from `slmlist` to `s0`, where equivalent SLMs coincide.
    """
    cA = crystA[0].T
    cB = crystB[0].T
    gA = get_pure_rotation(crystA)
    gB = get_pure_rotation(crystB)
    hA0, hB0, q0 = s0
    x0 = np.transpose(np.dot((gB @ hB0) @ q0, la.inv(gA @ hA0)), axes=[2,0,1,3]).reshape(-1,9)
    _, i = np.unique(x0.round(decimals=4), axis=0, return_index=True)
    cl0 = cB @ x0[i,:].reshape(-1,3,3) @ la.inv(cA)
    ss = cB @ slmlist[:,1,:,:] @ slmlist[:,2,:,:] @ la.inv(cA @ slmlist[:,0,:,:])
    dlist = np.amin(la.norm(ss.reshape(-1,1,9), cl0.reshape(1,-1,9), axis=2),axis=1)
    return dlist

def orient_matrix(vi: ArrayLike, vf: ArrayLike, wi: ArrayLike, wf: ArrayLike) -> NDArray[np.float64]:
    """Rotation matrix `r` such that `r @ vi` || `vf` and `r @ wi` || `wf`.

    Parameters
    ----------
    vi, vf : (3,) array_like
        Vectors (cartesian coordinates) that satisfy `r @ vi` || `vf`.
    wi, wf : (3,) array_like
        Vectors (cartesian coordinates) that satisfy `r @ wi` || `wf`.
        
    Returns
    -------
    r : (3, 3) array
        A rotation matrix representing the given orientation relationship.
    """
    b = np.array([vi, np.cross(vi, wi), np.cross(vi, np.cross(vi, wi))]).T
    c = np.array([vf, np.cross(vf, wf), np.cross(vf, np.cross(vf, wf))]).T
    b = b * la.norm(b, axis=0).reshape(1,3) ** -1
    c = c * la.norm(c, axis=0).reshape(1,3) ** -1
    r = c @ b.T
    return r

def deviation_angle(
    crystA: Cryst,
    crystB: Cryst,
    slmlist: Union[List[SLM], NDArray[np.int32]],
    r: NDArray[np.float64],
    uspfix: bool = False
) -> NDArray[np.float64]:
    """Calculate how much each SLM in `slmlist` differ from a given orientation relationship.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slmlist : list of slm
        A list of SLMs, each represented by a triplet of integer matrices like `(hA, hB, q)`.
    r : (3, 3) array
        A rotation matrix representing the given orientation relationship.
    uspfix : bool, optional
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
    hA, hB, q = zip(*slmlist)
    s = cB @ np.array(hB) @ np.array(q) @ la.inv(cA @ np.array(hA))
    u, sigma, vT = la.svd(s)
    rS = u @ vT
    if uspfix:
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

def save_scatter(
    filename : str,
    rmsslist: ArrayLike,
    rmsdlist: ArrayLike,
    colorlist: ArrayLike,
    cmap : colors.Colormap = plt.cm.get_cmap('viridis'),
    cbarlabel: str = None
) -> None:
    """Scatter plot of the CSMs with colorbar.

    Parameters
    ----------
    filename : str
        The filename of the saved plot.
    rmsslist : (N,) array_like
        The root-mean-square strain of each CSM.
    rmsdlist : (N,) array_like
        The root-mean-square distance of each CSM.
    colorlist : (N,) array_like
        Some quantity of each CSM, which is to be colored.
    cmap : `matplotlib.colors.Colormap`, optional
        The colormap to use. Default is `plt.cm.get_cmap('viridis')`.
    cbarlabel : str, optional
        The label of the colorbar. Default is None, in which case the filename is used.
    """
    rmsslist = np.array(rmsslist)
    rmsdlist = np.array(rmsdlist)
    colorlist = np.array(colorlist)
    plt.figure()
    ax = plt.subplot()
    ind = np.argsort(colorlist)[::-1]
    sc = plt.scatter(rmsslist[ind], rmsdlist[ind], marker='d', c=colorlist[ind], cmap=cmap, s=20)
    ind0 = colorlist==0
    n0 = np.sum(ind0)
    if n0 >= 1:
        print(f"\nThere are {n0:d} CSMs (indices: {', '.join(np.nonzero(ind0)[0].astype(str).tolist())}) with {cbarlabel}=0, emphasized by pink stars in the plot.")
        plt.scatter(rmsslist[ind0], rmsdlist[ind0], marker='*', color=(1.0,0.75,0.95), s=12)
    plt.xlabel("Root-mean-square strain", fontsize=15)
    plt.ylabel("RMSD (Å)", fontsize=15)
    plt.xlim(0, np.amax(rmsslist) * 1.05)
    plt.ylim(min(0, np.amin(rmsdlist) - 0.1), np.amax(rmsdlist) + 0.1)
    if colorlist.dtype == int:
        cbar = plt.colorbar(sc, aspect=40, ticks=np.unique(colorlist))
    else:
        cbar = plt.colorbar(sc, aspect=40)
    cbar.set_label(cbarlabel if cbarlabel != None else filename.split('.')[0], fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    plt.savefig(f"{filename}", bbox_inches='tight')
    return

def save_interpolation(
    filename: str,
    crystA_sup: Cryst,
    crystB_sup: Cryst,
    images: int = 10,
    crystname: Union[str, None] = None
) -> None:
    """
    Save the linear interpolation between `crystA` and `crystB` to a single XDATCAR file.
    
    Parameters
    ----------
    filename : str
        The name of the file to save, must not already exist in current directory.
    crystA_sup, crystB_sup : cryst
        The initial and final crystal structures with specified atomic correspondence, usually obtained by `minimize_rmsd`.
    images : int, optional
        Number of images to generate. Default is 10.
    crystname : str, optional
        A system description to write to the comment line of the POSCAR file. If `crystname = None`, `filename` will be used.
    """
    if not (crystA_sup[1] == crystB_sup[1]).all():
        raise ValueError("Atomic species of crystA and crystB must be the same.")
    if type(images) != int or images < 1:
        raise ValueError("Number of images must be a positive integer.")
    
    cA = crystA_sup[0].T
    cB = crystB_sup[0].T
    pA = crystA_sup[2].T
    pB = crystB_sup[2].T
    s = cB @ la.inv(cA)
    if not ((s.T - s).round(decimals=4) == 0).all():
        print(f"Warning: Extra rotation detected when interpolating crystals, which is removed in {filename}.")
    _, sigma, vT = la.svd(s)
    crystlist = []
    tlist = np.linspace(0, 1, images+2)
    for t in tlist:
        c = vT.T @ np.diag(sigma ** t) @ vT @ cA
        p = pA * (1-t) + pB * t
        crystlist.append((c.T, crystA_sup[1], p.T))
    
    content = crystname
    for i in range(images+2):
        if i > 0: content += '\n'
        if crystname: content += crystname
        else: content += filename.split(sep='.')[0]
        content += '\n1.0\n'
        content += '\n'.join(f'{v[0]:.12f}\t{v[1]:.12f}\t{v[2]:.12f}' for v in crystlist[i][0].tolist())
        species_name, species_counts = species_poscar_format(crystA_sup[1])
        content += '\n' + ' '.join(species_name.tolist())
        content += '\n' + ' '.join(str(n) for n in species_counts.tolist())
        content += f'\nDirect configuration= {i+1:.0f}\n'
        content += '\n'.join(f'{p[0]:.12f}\t{p[1]:.12f}\t{p[2]:.12f}' for p in crystlist[i][2].tolist())
        
    f = open(filename, mode='x')
    f.write(content)
    f.close()
    return