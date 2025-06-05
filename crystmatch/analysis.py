"""
Save, load, analyze and visualize CSMs.
"""

from .utilities import *

__all__ = ["decompose_cryst", "conventional_cryst", "csm_to_cryst", "csm_distance", "primitive_shuffle",
            "orient_matrix", "rot_usp", "deviation_angle", "visualize_slmlist", "visualize_pctlist", "visualize_csm",
            "save_npz", "load_npz", "unzip_csm", "save_poscar", "save_xdatcar"]

np.set_printoptions(suppress=True)
Table = namedtuple("Table", ["data", "header"])

def decompose_cryst(cryst_sup, cryst_prim=None, tol=1e-3):
    """Compute the rotation `r` and the integer transformation matrix `m` such that `c_sup = r @ c_prim @ m`
    """
    if cryst_prim is None:
        cryst_prim = primitive_cryst(cryst_sup, tol=tol)

    # c_conv @ m0 = r0 @ c_prim
    sym0 = get_symmetry_dataset(cryst_to_spglib(cryst_prim), symprec=tol)
    m0 = sym0.transformation_matrix
    r0 = sym0.std_rotation_matrix

    # c_conv @ m1 = r1 @ c_sup
    sym1 = get_symmetry_dataset(cryst_to_spglib(cryst_sup), symprec=tol)
    m1 = sym1.transformation_matrix
    r1 = sym1.std_rotation_matrix

    # c_sup = (r1.inv @ r0) @ c_prim @ (m0.inv @ m1)
    r = r1.T @ r0
    m = la.inv(m0) @ m1
    if (np.abs(m - m.round()) < tol).all(): m = m.round().astype(int)
    return r, m

def conventional_cryst(cryst_sup, return_primitive=False, tol=1e-3):
    """Rotate the crystal structure and change its lattice basis such that `c` is lower triangular.
    """
    c_sup = cryst_sup[0].T
    p_sup = cryst_sup[2].T
    r, m = decompose_cryst(cryst_sup, tol=tol)
    h, q = hnf(m.round().astype(int), return_q=True)
    c_sup_tri = la.inv(r) @ c_sup @ la.inv(q)
    p_sup_tri = (q @ p_sup) % 1.0
    cryst_sup_tri = (c_sup_tri.T, cryst_sup[1], p_sup_tri.T)
    if return_primitive:
        cell_sup, species_dict = cryst_to_spglib(cryst_sup, return_dict=True)
        cryst_tri = spglib_to_cryst(refine_cell(cell_sup, symprec=tol), species_dict)
        if not np.allclose(cryst_tri[0].T @ h, c_sup_tri):
            raise ValueError("Error in triangularization. Please report this bug to wfc@pku.edu.cn.")
        return cryst_sup_tri, cryst_tri
    else:
        return cryst_sup_tri

def csm_to_cryst(crystA: Cryst, crystB: Cryst, slm: SLM, p: NDArray[np.int32], ks: NDArray[np.int32],
    orientation: Literal['norot', 'uspfixed'] = 'norot',
    min_t0: bool = True, weight_func: Union[Dict[str, float], None] = None, l: float = 2
) -> Union[Tuple[Cryst, Cryst], Tuple[Cryst, Cryst, Cryst]]:
    """Convert the integer arrays representation `(slm, p, translations)` of a CSM to a pair of crysts.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slm : slm
        The SLM of the CSM.
    p : (Z, ) array of ints
        The permutaion of the shuffle.
    ks : (3, Z) array of ints
        The lattice-vector translations of the shuffle.
    orientation : str, optional
        The orientation of the final structure, either 'norot' or 'uspfixed', which means that the deformation is rotation-free \
            or fixing the uniformly scaled plane (USP). When 'uspfixed', two final structures are returned since there are two USPs. \
            Default is 'norot'. 
    min_t0 : bool, optional
        Whether to use optimal translation of `crystB_sup_norot` to minimize the shuffle distance. Default is True.
    weight_func : dict, optional
        The weight function used when minimizing the shuffle distance, with keys as atomic species. \
            Default is None, which means all atoms have the same weight.
    l : float, optional
        The l-norm to be used for distance calculation, must not be less than 1. Default is 2.
        
    Returns
    -------
    crystA_sup : cryst
        The initial structure.
    crystB_sup_norot : cryst
        The final structure, whose lattice vectors and atoms are matched to `crystA_sup` according to the CSM, with rotation-free orientation.
    crystB_sup_uspfixed_1, crystB_sup_uspfixed_2 : cryst
        The final structure, whose lattice vectors and atoms are matched to `crystA_sup` according to the CSM, with uniformly scaled plane fixed.
    """
    crystA_sup, crystB_sup, c_sup_half, _, _ = create_common_supercell(crystA, crystB, slm)
    species = crystA_sup[1]
    if not (crystB_sup[1][p] == species).all():
        raise ValueError("The input permutation does not preserve the atomic species.")
    cB = crystB_sup[0].T
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    if min_t0:
        weights = [weight_func[s] for s in species] if weight_func else None
        pB_sup = pB_sup - pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weights, l=l, return_t0=True)[1]
    if orientation == 'norot':
        return crystA_sup, (cB.T, species, pB_sup.T)
    elif orientation == 'uspfixed':
        rH = rot_usp(deformation_gradient(crystA_sup, crystB_sup, slm))
        return crystA_sup, ((rH @ cB).T, species, pB_sup.T), ((rH.T @ cB).T, species, pB_sup.T)
    else: raise ValueError("'deformation' must be 'norot' or 'uspfixed'.")

def csm_distance(crystA, crystB, slm, p, ks, weight_func=None, l=2, min_t0=True, return_t0=False):
    """Return the shuffle distance of a CSM.
    """
    crystA_sup, crystB_sup, c_sup_half, _, _ = create_common_supercell(crystA, crystB, slm)
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    species = crystA_sup[1]
    weights = [weight_func[s] for s in species] if weight_func else None
    return pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weights, l=l, min_t0=min_t0, return_t0=return_t0)

def _is_translation_element(t, pA, pB, p, ks):
    # since pA -> pB[:,p] + ks, we expect pA[:,pAt] + ksAt -> pB[:,p[pAt]] + ks[:,pAt] + ksAt == pB[:,pBt[p]] + ks + ksBt
    pAt = np.nonzero((((pA.reshape(3,1,-1) + t.reshape(3,1,1) - pA.reshape(3,-1,1)).round(7) % 1.0) == 0).all(axis=0))[1]
    ksAt = (pA + t.reshape(3,1) - pA[:,pAt]).round().astype(int)
    pBt = np.nonzero((((pB.reshape(3,1,-1) + t.reshape(3,1,1) - pB.reshape(3,-1,1)).round(7) % 1.0) == 0).all(axis=0))[1]
    ksBt = (pB[:,p] + t.reshape(3,1) - pB[:,pBt[p]]).round().astype(int)
    return (p[pAt] == pBt[p]).all() and (ks[:,pAt] + ksAt == ks + ksBt).all()

def primitive_shuffle(crystA, crystB, slm0, p0, ks0):
    """Identify the primitive (with minimal multiplicity) IMT and PCT representations of a CSM.
    """
    crystA_sup0, crystB_sup0, _, mA0, mB0 = create_common_supercell(crystA, crystB, slm0)
    pA_sup0 = crystA_sup0[2].T
    pB_sup0 = crystB_sup0[2].T
    
    # identify translation elements
    tA = int_vec_inside(mA0)
    tA = la.inv(mA0) @ tA[:,(tA != 0).any(axis=0)]
    tB = int_vec_inside(mB0)
    tB = la.inv(mB0) @ tB[:,(tB != 0).any(axis=0)]
    common_list = np.array([t for t in set([tuple(x) for x in tA.T.round(7)]) & set([tuple(y) for y in tB.T.round(7)])])
    f = np.array([t for t in common_list if _is_translation_element(t, pA_sup0, pB_sup0, p0, ks0)]).T

    # identify primitive IMT
    hA0, hB0, q0 = slm0
    hA = hnf((mA0 @ np.hstack([np.eye(3), f.reshape(3,-1)])).round().astype(int), return_q=False)[:,:3]
    hB = hnf((mB0 @ np.hstack([np.eye(3), f.reshape(3,-1)])).round().astype(int), return_q=False)[:,:3]
    q = (la.inv(hB) @ hB0 @ q0 @ la.inv(hA0) @ hA).round().astype(int)
    slm = np.array([hA, hB, q])

    # calculate new PCT
    crystA_sup, crystB_sup, _, mA, mB = create_common_supercell(crystA, crystB, slm)
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    z = pA_sup.shape[1]
    p = np.ones(z, dtype=int) * z
    ks = np.zeros((3,z), dtype=int)

    # find primitive-cell atoms in supercells
    iA, iA0 = np.nonzero((((la.inv(mA) @ mA0 @ pA_sup0).reshape(3,1,-1) - pA_sup.reshape(3,-1,1)).round(7) % 1.0 == 0).all(axis=0))
    iB, iB0 = np.nonzero((((la.inv(mB) @ mB0 @ pB_sup0).reshape(3,1,-1) - pB_sup.reshape(3,-1,1)).round(7) % 1.0 == 0).all(axis=0))
    iB0_to_iB = np.zeros_like(p0)
    iB0_to_iB[iB0] = iB
    for i in range(z):
        pi = iB0_to_iB[p0[iA0[iA == i]]]
        if np.diff(pi).any(): raise AssertionError("'p0' is not a class function. Please report this bug to wfc@pku.edu.cn if you see this message.")
        p[i] = pi[0]
        i0 = iA0[iA == i][0]
        ksi = (la.inv(mB) @ mB0 @ (pB_sup0[:,p0[i0]] + ks0[:,i0]) - la.inv(mA) @ mA0 @ pA_sup0[:,i0]) - (pB_sup[:,p[i]] - pA_sup[:,i])
        if (ksi.round(7) % 1.0 != 0).any(): raise AssertionError("'ks0' is not integer. Please report this bug to wfc@pku.edu.cn if you see this message.")
        ks[:,i] = ksi.round().astype(int)
    return slm, p, ks
    
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

def rot_usp(s):
    """The rotation that the uniformly scaled plane undergoes.
    """
    _, sigma, vT = la.svd(s, compute_uv=True)
    eps = np.array([[[0,0,0],[0,0,-1],[0,1,0]],[[0,0,1],[0,0,0],[-1,0,0]],[[0,-1,0],[1,0,0],[0,0,0]]])
    v_cross = np.tensordot(vT[:,1,:], eps, axes=(1,0))
    s1 = sigma[:,0]
    s2 = sigma[:,1]
    s3 = sigma[:,2]
    theta = np.arctan(np.sqrt((s1**2 - s2**2) * (s2**2 - s3**2)) / (s1 * s3 + s2**2))
    r = np.eye(3).reshape(1,3,3) + np.sin(theta).reshape(-1,1,1) * v_cross + (1 - np.cos(theta)).reshape(-1,1,1) * (v_cross @ v_cross)
    return r

def deviation_angle(
    crystA: Cryst,
    crystB: Cryst,
    slmlist: Union[List[SLM], NDArray[np.int32]],
    r: NDArray[np.float64],
    orientation: Literal['norot', 'uspfixed']
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
    orientation : str, optional
        The orientation of the final structure, either 'norot' or 'uspfixed', which means that the deformation is rotation-free \
            or fixing the uniformly scaled plane (USP). Default is 'norot'. When 'uspfixed', two final structures are returned since there are two USPs.

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
    u, _, vT = la.svd(s, compute_uv=True)
    rS = u @ vT
    if orientation == 'norot':
        anglelist = np.arccos(np.clip(0.5 * (-1 + np.amax(np.trace(np.dot(la.inv(r_equiv), rS), axis1=1, axis2=3), axis=0)), -1, 1))
    elif orientation == 'uspfixed':
        rH = rot_usp(s)
        rS = np.array([rS @ rH, rS @ rH.T])       # rS.shape = (2, ..., 3, 3)
        anglelist = np.arccos(np.clip(0.5 * (-1 + np.amax(np.trace(np.dot(la.inv(r_equiv), rS), axis1=1, axis2=4), axis=(0,1))), -1, 1))
    else: raise ValueError("'deformation' must be 'norot' or 'uspfixed'.")
    return anglelist.round(decimals=7)

def visualize_slmlist(
    filename : Union[str, None],
    wlist: ArrayLike,
    d0list: ArrayLike,
    colorlist: ArrayLike = None,
    cmap : colors.Colormap = plt.get_cmap('viridis'),
    cbarlabel: str = None
) -> None:
    """Scatter plot of the CSMs with colorbar.

    Parameters
    ----------
    filename : str
        The filename of the saved plot. If None, the plot is shown on screen.
    rmsslist : (N,) array_like
        The root-mean-square strain of each CSM.
    rmsdlist : (N,) array_like
        The root-mean-square distance of each CSM.
    colorlist : (N,) array_like, optional
        Some quantity of each CSM, which is to be colored. If None, do not color the points.
    cmap : `matplotlib.colors.Colormap`, optional
        The colormap to use. Default is `plt.cm.get_cmap('viridis')`.
    cbarlabel : str, optional
        The label of the colorbar. Default is None, in which case the filename is used.
    """
    wlist = np.array(wlist)
    d0list = np.array(d0list)
    if colorlist is None:
        colorlist = np.zeros_like(wlist)
        n0 = 0
    else:
        colorlist = np.array(colorlist)
        ind0 = colorlist==0
        n0 = np.sum(ind0)
    plt.figure()
    ax = plt.subplot()
    ind = np.argsort(colorlist)[::-1]
    sc = plt.scatter(wlist[ind], d0list[ind], marker='d', c=colorlist[ind], cmap=cmap, s=20)
    if n0 >= 1:
        print(f"\nThere are {n0:d} CSMs (indices: {', '.join(np.nonzero(ind0)[0].astype(str).tolist())}) with {cbarlabel}=0, emphasized by pink stars in the plot.")
        plt.scatter(wlist[ind0], d0list[ind0], marker='*', color=(1.0,0.75,0.95), s=12)
    plt.xlabel("w (same units as input)", fontsize=13)
    plt.ylabel("Shuffle distance (Å)", fontsize=13)
    plt.xlim(0, np.amax(wlist) * 1.05)
    plt.ylim(min(0, np.amin(d0list) - 0.1), np.amax(d0list) + 0.1)
    if (colorlist == 0).all(): pass
    elif colorlist.dtype == int: cbar = plt.colorbar(sc, aspect=40, ticks=np.unique(colorlist))
    else: cbar = plt.colorbar(sc, aspect=40)
    if (colorlist != 0).any() and cbarlabel is not None: cbar.set_label(cbarlabel, fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    if filename: plt.savefig(f"{filename}", bbox_inches='tight')
    else: plt.show()

def visualize_pctlist(filename, pctlist, dlist):
    _, ind = np.unique([pct[:,0] for pct in pctlist], axis=0, return_inverse=True)
    d_min = []
    for i in range(ind.max()+1):
        d_min.append(dlist[ind==i].min())
    a = np.argsort(d_min)
    aa = np.zeros_like(a)
    aa[a] = np.arange(len(a))
    _, ax = plt.subplots()
    ax.scatter(aa[ind], dlist, s=20, linewidths=1, c='r', marker='x')
    ax.set_xticks(np.arange(len(a)))
    ax.set_xticklabels([str(i+1) for i in range(len(a))])
    ax.set_xlabel("Permutation index", fontsize=13)
    ax.set_ylabel("Shuffle distance (Å)", fontsize=13)
    ax.grid(True, linestyle=':')
    if filename: plt.savefig(f"{filename}", bbox_inches='tight')
    else: plt.show()

def cell_lines(c):
    line1 = c @ np.array([[0,0,0,0,0,1,1,0], [0,1,1,0,0,0,1,1], [0,0,1,1,0,0,0,0]])
    line2 = c @ np.array([[0,1,1], [0,0,0], [1,1,0]])
    line3 = c @ np.array([[0,1,1], [1,1,1], [1,1,0]])
    line4 = c @ np.array([[1,1], [0,1], [1,1]])
    return line1, line2, line3, line4

def colors_dark(i):
    """Returns the i-th dark color.

    Parameters
    ----------
    i : int
        Index of the color to return; must be between 0 and 8.

    Returns
    -------
    rgba : (4, ) tuple of float
        RGBA values of the i-th dark color.
    """
    if i > 8: raise ValueError('i must be between 0 and 8')
    return plt.get_cmap('Set1')(0.1 * i + 0.05)

def colors_light(i):
    """Returns the i-th light color.

    Parameters
    ----------
    i : int
        Index of the color to return; must be between 0 and 8.

    Returns
    -------
    rgba : (4, ) tuple of float
        RGBA values of the i-th light color.
    """
    if i > 8: raise ValueError('i must be between 0 and 8')
    return plt.get_cmap('Pastel1')(0.1 * i + 0.05)

def visualize_csm(crystA, crystB, slm, p, ks, weights=None, l=2, cluster_size=1.5, show_conventional=True, tol=1e-3):
    """Use with `%matplotlib widget` in Jupyter notebook (need to install `ipympl`) to interactively visualize the shuffling process.
    """
    z = len(p)
    crystA_sup, crystB_sup, c_half, _, _ = create_common_supercell(crystA, crystB, slm)
    species, numbers, counts = np.unique(crystA_sup[1], return_inverse=True, return_counts=True)
    if numbers.max() > 8: raise ValueError("Too many (>8) atoms species to visualize.")
    cA = crystA_sup[0].T
    cB = crystB_sup[0].T
    pA = crystA_sup[2].T
    pB = crystB_sup[2].T
    t_cl = pA // 1.0
    t_cen = np.average(pA % 1.0, weights=weights, axis=1, keepdims=True) - np.ones((3,1)) * 0.5
    t0 = pct_distance(c_half, pA, pB, p, ks, weights=weights, l=l, return_t0=True)[1]
    
    fig = plt.figure(figsize=(8,4))
    for left in [True, False]:
        c = cA if left else cB
        cooA = c @ (pA - t_cl - t_cen)
        cooB = c @ (pB[:,p] + ks - t_cl + t0 - t_cen)
        shuf = cooB - cooA
        center = c @ np.ones((3,1)) * 0.5
        radius = la.norm(c @ np.mgrid[0:2,0:2,0:2].reshape(3,-1) - center, axis=0).max() * cluster_size
        dcoo = c @ np.mgrid[-2:3,-2:3,-2:3].reshape(3,-1)
        dcoo = dcoo[:,(dcoo!=0).any(axis=0)]

        ax = fig.add_subplot(1, 2, 1 if left else 2, projection='3d')
        ax.set_facecolor('white')
        ax.set_proj_type('ortho')
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

        for line in cell_lines(c):
            ax.plot(line[0], line[1], line[2], color='#cccccc', linestyle='--', linewidth=0.8)
        for i in range(numbers.max()+1):
            ind = numbers == i
            # cell content
            coo = cooA[:,ind] if left else cooB[:,ind]
            ax.scatter(coo[0], coo[1], coo[2], color=colors_dark(i), s=32)
            ax.quiver(cooA[0,ind], cooA[1,ind], cooA[2,ind], shuf[0,ind], shuf[1,ind], shuf[2,ind], color='#bbbbbb', arrow_length_ratio=0.5)
            # cluster
            coo = (coo.reshape(3,-1,1) + dcoo.reshape(3,1,-1)).reshape(3,-1)
            coo = coo[:, la.norm(coo - center, axis=0) < radius]
            ax.scatter(coo[0], coo[1], coo[2], color=colors_light(i), s=8)
        if show_conventional:
            cryst = crystA_sup if left else crystB_sup
            sym = get_symmetry_dataset(cryst_to_spglib(cryst), symprec=tol)
            c_conv = sym.std_rotation_matrix.T @ sym.std_lattice.T
            ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), c_conv[0], c_conv[1], c_conv[2], color=['#666666', '#888888', '#aaaaaa'], arrow_length_ratio=0.2)

        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_zlim3d()
        ax.set_box_aspect((xlim[1]-xlim[0], ylim[1]-ylim[0], zlim[1]-zlim[0]))
        ax.text2D(0.05, 0.95, f"{'Initial' if left else 'Final'} structure ({sym.international})", fontsize=10, transform=ax.transAxes)
        if show_conventional: ax.text2D(0.05, 0.91, f"with conventional basis represented by gray arrows", color='#777777', fontsize=8, transform=ax.transAxes)
        ax.text2D(0.05, 0.10, f"Atoms in the dashed supercell (shuffle lattice):", color='#999999', fontsize=9, transform=ax.transAxes)
        for i, s in enumerate(species):
            ax.text2D(0.05 + 0.1 * i, 0.05, f"{counts[i]:d}{s}", color=colors_dark(i), fontsize=9, transform=ax.transAxes)

    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    plt.show()

def _species_to_poscar(species: NDArray[np.str_]) -> Tuple[NDArray[np.str_], NDArray[np.int32]]:
    _, sp_idx, sp_inv, sp_counts = np.unique(species, return_index=True, return_inverse=True, return_counts=True)
    if np.sum(np.diff(sp_inv) != 0) != sp_idx.shape[0] - 1:
        raise ValueError("Species array is ill-sorted. Please report this bug to wfc@pku.edu.cn if you see this message.")
    return species[np.sort(sp_idx)], sp_counts

def save_npz(filename, crystA, crystB, slmlist, slm_ind, pct_arrs, table):
    metadata = np.array([table.header, crystA, crystB], dtype=object)
    if not (type(pct_arrs[0]) == str or pct_arrs[0] is None): raise ValueError("The first element of 'pct_arrs' should be a placeholder.")
    np.savez(unique_filename("Saving CSMs to", filename), *pct_arrs, metadata=metadata, slmlist=slmlist, slm_ind=slm_ind, table=table.data)
    return

def load_npz(filename, verbose=True):
    npz = np.load(filename, allow_pickle=True)
    header, crystA, crystB = npz['metadata']
    slmlist = npz['slmlist']
    slm_ind = npz['slm_ind']
    max_mu = len([f for f in npz.files if f.startswith('arr_')]) - 1
    pct_arrs = [npz[f'arr_{i}'] for i in range(max_mu+1)]
    data = npz['table']
    if verbose:
        print(f"A total of {len(slmlist)} distinct SLMs and {len(slm_ind)} CSMs were loaded from {filename}.")
        print(f"\tmu  {' '.join(f'{mu:5d}' for mu in range(1,max_mu+1))}")
        print(f"\t#CSM{' '.join(f'{arr.shape[0]:5d}' for arr in pct_arrs[1:])}")
    return crystA, crystB, slmlist, slm_ind, pct_arrs, Table(data, header)

def unzip_csm(csm_ind, crystA, crystB, slmlist, slm_ind, pct_arrs):
    slm = slmlist[slm_ind[csm_ind]]
    mu = imt_multiplicity(crystA, crystB, slm)
    ind_arr = (imt_multiplicity(crystA, crystB, slmlist[slm_ind[:csm_ind]]) == mu).sum()
    p, ks = unzip_pct(pct_arrs[mu][ind_arr])
    return slm, p, ks

def save_poscar(
    filename: Union[str, None],
    cryst: Cryst,
    crystname: Union[str, None] = None
) -> None:
    """
    Save the crystal structure to a POSCAR file.

    Parameters
    ----------
    filename : str
        The name of the file to save, must not already exist in current directory. If `filename = None`, a string will be returned instead.
    cryst : cryst
        The crystal structure to be saved, consisting of the lattice vectors, species, and positions.
    crystname : str, optional
        A system description to write to the comment line of the POSCAR file. If `crystname = None`, `filename` will be used.
    """
    species_name, species_counts = _species_to_poscar(cryst[1])
    if crystname is not None: content = crystname
    else: content = ""
    content += "\n1.0\n"
    content += "\n".join(f"{v[0]:.12f}\t{v[1]:.12f}\t{v[2]:.12f}" for v in cryst[0].tolist())
    content += "\n" + " ".join(species_name.tolist())
    content += "\n" + " ".join(str(n) for n in species_counts.tolist())
    content += "\nDirect\n"
    content += "\n".join(f"{p[0]:.12f}\t{p[1]:.12f}\t{p[2]:.12f}" for p in cryst[2].tolist())
    if filename is not None:
        f = open(filename, mode='x')
        f.write(content)
        f.close()
        return
    else:
        return content

def save_xdatcar(
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
        species_name, species_counts = _species_to_poscar(crystA_sup[1])
        content += '\n' + ' '.join(species_name.tolist())
        content += '\n' + ' '.join(str(n) for n in species_counts.tolist())
        content += f'\nDirect configuration= {i+1:.0f}\n'
        content += '\n'.join(f'{p[0]:.12f}\t{p[1]:.12f}\t{p[2]:.12f}' for p in crystlist[i][2].tolist())
        
    f = open(filename, mode='x')
    f.write(content)
    f.close()
    return