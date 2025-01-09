"""
Enumerate SLMs and CSMs.
"""

from .utilities import *
from time import time
from scipy.optimize import brentq, linear_sum_assignment, basinhopping
from scipy.stats.qmc import Sobol
from scipy.spatial.transform import Rotation

np.set_printoptions(suppress=True)
Cryst = Tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]

def cong_slm(slm: Union[SLM, NDArray[np.int32]], gA: NDArray[np.int32], gB: NDArray[np.int32]) -> tuple[SLM, int]:
    """The representative SLM of the congruence class of `slm`.

    Parameters
    ----------
    slm : slm
        `(hA, hB, q)`, representing a SLM.
    gA : (..., 3, 3) array of ints
        The rotation group of the initial crystal structure, whose elements are \
            integer matrices under fractional coordinates.
    gB : (..., 3, 3) array of ints
        The rotation group of the final crystal structure, whose elements are \
            integer matrices under fractional coordinates.

    Returns
    -------
    ss : slm
        The representative of the congruence class of `slm`.
    len_cl : int
        The size of the congruence class of `slm`.
    """
    hA, hB, q = slm
    cl = np.transpose(np.dot((gB @ hB) @ q, la.inv(gA @ hA)), axes=[2,0,1,3]).reshape(-1,9)
    cl, i = np.unique(cl.round(decimals=4), axis=0, return_index=True)
    iA, iB = np.unravel_index(i[0], (gA.shape[0], gB.shape[0]))
    hAA, qA = hnf_int(gA[iA] @ hA)
    hBB, qB = hnf_int(gB[iB] @ hB)
    ss = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return ss, len(cl)

def enumerate_slm(
    crystA: Cryst, crystB: Cryst, mu: int, kappa_max: float, tol: float = 1e-3,
    kappa: Callable[[NDArray[np.float64]], NDArray[np.float64]] = rmss,
    likelihood_ratio: float = 1e2, verbose: int = 0
) -> List[SLM]:
    """Enumerating all SLMs of multiplicity `mu` with `kappa` smaller than `kappa_max`.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    mu : int
        The multiplicity of SLMs to enumerate.
    kappa_max : float
        A positive threshold value of `kappa` to determine the range of singular values to generate.
    tol : float, optional
        The tolerance for determining the pure rotation group of the crystal structures. Default is 1e-3.
    kappa : callable, optional
        A function that quantifies the strain of a matrix according to its singular values. \
            Default is `rmss`, which computes the root-mean-square strain.
    likelihood_ratio : float, optional
        The expected likelihood ratio of the enumeration being complete and incomplete. Default is 1e2.
    verbose : int, optional
        The level of detail of printing. `0` means no print.

    Returns
    -------
    slmlist : list of slm
        Contains triplets of integer matrices like `(hA, hB, q)`, representing incongruent SLMs.
    """
    assert mu >= 1
    cA = crystA[0].T
    cB = crystB[0].T
    zA = crystA[2].shape[0]
    zB = crystB[2].shape[0]
    hA = hnf_list(np.lcm(zA,zB) // zA * mu)
    hB = hnf_list(np.lcm(zA,zB) // zB * mu)
    hA_inv = la.inv(hA)
    hB_inv = la.inv(hB)
    gA = get_pure_rotation(crystA, tol=tol)
    gB = get_pure_rotation(crystB, tol=tol)
    max_prob_ratio = 20
    # Compute the sampling domain of `s0`s.
    det_s = zA / zB * la.det(cB) / la.det(cA)
    if kappa.__name__ == 'rmss':
        max_strain = kappa_max * 3**0.5
    else:
        print(f"\nWarning: kappa function {kappa.__name__} is not RMSS. I hope you know what you are doing.")
        diff_kappa = lambda x: kappa(np.array([x, (det_s / x) ** 0.5, (det_s / x) ** 0.5])) - kappa_max
        a = brentq(diff_kappa, 1e-2, det_s**(1/3))
        b = brentq(diff_kappa, det_s**(1/3), 1e2)
        max_strain = max(abs(a-1), abs(b-1))
    # Enumerate SLMs.
    print(f"Enumerating SLMs (mu = {mu:d}, {kappa.__name__} <= {kappa_max:.4f}) ...")
    if verbose >= 1:
        print(f"\tprototype sampling domain: SO(3) (±{max_strain:.2f} for each matrix element)")
        print(f"\tassumed maximum probability ratio among classes: {max_prob_ratio:.1f}")
        print(f"\texpected likelihood ratio: {likelihood_ratio:.1f}")
        print("\tnum_s0\tm\tm*\telapsed_time(s)")
    slmlist = []
    iter = 0
    num_s0 = 0
    m = 0
    t = time()
    sobol_seq = Sobol(12)
    while m <= (1 + len(slmlist) * max_prob_ratio) * np.log(likelihood_ratio):    # or True for Debug!
        # Sampling `s0`s around SO(3).
        rand_num = sobol_seq.random(2)
        q0 = np.sqrt(1 - rand_num[:,0]) * np.sin(2 * np.pi * rand_num[:,1])
        q1 = np.sqrt(1 - rand_num[:,0]) * np.cos(2 * np.pi * rand_num[:,1])
        q2 = np.sqrt(rand_num[:,0]) * np.sin(2 * np.pi * rand_num[:,2])
        q3 = np.sqrt(rand_num[:,0]) * np.cos(2 * np.pi * rand_num[:,2])
        s0 = Rotation.from_quat(np.array([q0,q1,q2,q3]).T).as_matrix() + (2 * rand_num[:,3:].reshape(-1,3,3) - 1) * max_strain
        # Round `s0`s to nearest integer matrix triplets.
        q = (np.transpose(np.dot(np.dot(hB_inv, la.inv(cB) @ s0 @ cA), hA), axes=[2,3,0,1,4])).round().astype(int)
        index1 = la.det(q.reshape(-1,3,3)).round().astype(int) == 1
        # Check determinants of `q`s.
        index1 = np.nonzero(index1)[0]
        index1 = np.unravel_index(index1, (s0.shape[0], hA.shape[0], hB.shape[0]))
        # Check strains of `slm`s.
        index2 = np.nonzero(kappa(la.svd(cB @ hB[index1[2]] @ q[index1] @ hA_inv[index1[1]] @ la.inv(cA),
                                        compute_uv=False)) < kappa_max)[0]                                          
        for i in index2:
            slm = (hA[index1[1][i]], hB[index1[2][i]], q[index1[0][i], index1[1][i], index1[2][i]])
            slm, _ = cong_slm(slm, gA, gB)
            repeated = False
            for j in range(len(slmlist)):
                if (slm[0] == slmlist[j][0]).all() and (slm[1] == slmlist[j][1]).all() and (slm[2] == slmlist[j][2]).all():
                    # `slm` is repeated.
                    m = m + 1
                    repeated = True
                    break
            if not repeated:
                # `slm` is new.
                slmlist.append(slm)
                m = 0
        num_s0 += s0.shape[0]
        if verbose >= 2:
            print(f"\t{num_s0}\t{m}\t{(1 + len(slmlist) * max_prob_ratio) * np.log(likelihood_ratio):d}\t{time()-t:.2f}")
        iter = iter + 1
        if iter > 50 and len(slmlist) == 0: break
    print(f'\tFound {len(slmlist):d} incongruent SLMs (in {time()-t:.2f} seconds).')
    return slmlist

def minimize_rmsd_t(
    c: NDArray[np.float64],
    species: NDArray[np.str_],
    pA: NDArray[np.float64],
    pB: NDArray[np.float64]
) -> tuple[float, NDArray[np.int32], NDArray[np.int32]]:
    """Minimize the RMSD with fixed SLM, overall translation, variable permutation and lattice-vector translations.

    Parameters
    ----------
    c : (3, 3) array
        The matrix whose columns are cell vectors of both crystal structures.
    species : (N,) array of strs
        The array that specifies the type of each ion.
    pA, pB : (3, N) array
        The matrices whose columns are fractional coordinates of atoms.
        
    Returns
    -------
    rmsd : float
        The minimized RMSD.
    p : (Z, ) array of ints
        The permutation of the atoms in `pB` that minimizes the RMSD.
    ks : (3, Z) array of ints
        The lattice-vector translations of the atoms in `pB[:,p]` that minimizes the RMSD.
    """
    assert pA.shape[1] == pB.shape[1]
    # Determine the best lattice-vector translations for all possible atomic correspondences.
    shift = np.mgrid[-1:2, -1:2, -1:2].reshape(3,1,1,-1)
    relative_position = np.tensordot(c, (pA.reshape(3,-1,1,1) - pB.reshape(3,1,-1,1) - shift), axes=[1,0])
    d2_tensor = np.sum(relative_position ** 2, axis=0)
    index_shift = np.argmin(d2_tensor, axis=2)
    # Determine the best atomic correspondence.
    p = np.arange(pB.shape[1], dtype=int)
    for s in np.unique(species):
        is_s = species == s
        p[is_s] = p[is_s][linear_sum_assignment(np.take_along_axis(d2_tensor, index_shift[:,:,None], axis=2)[:,:,0][is_s][:,is_s])[1]]
    assert (np.sort(p) == np.arange(pB.shape[1])).all()
    ks = shift[:,0,0,index_shift[np.arange(pA.shape[1]),p]]
    rmsd = la.norm(c @ (pB[:,p] + ks - pA)) / np.sqrt(pA.shape[1])
    return rmsd, p, ks

def minimize_rmsd_tp(
    c: NDArray[np.float64],
    pA: NDArray[np.float64],
    pB: NDArray[np.float64]
) -> Tuple[float, NDArray[np.int32]]:
    """Minimize the RMSD with fixed SLM, overall translation, permutation, and variable lattice-vector translations.
    
    Parameters
    ----------
    c : (3, 3) array
        The matrix whose columns are cell vectors of both crystal structures.
    pA, pB : (3, N) array
        The matrices whose columns are fractional coordinates of atoms.
        
    Returns
    -------
    rmsd : float
        The minimized RMSD.
    ks : (3, Z) array of ints
        The lattice-vector translations (fractional coordinates) of the atoms in `pB` that minimizes the RMSD.
    """
    assert pA.shape[1] == pB.shape[1]
    shift = np.mgrid[-1:2, -1:2, -1:2].reshape(3,1,-1)
    relative_position = np.tensordot(c, (pA.reshape(3,-1,1) - pB.reshape(3,-1,1) - shift), axes=[1,0])
    d2_matrix = np.sum(relative_position ** 2, axis=0)
    ks = shift[:,0,np.argmin(d2_matrix, axis=1)]
    rmsd = la.norm(c @ (pB + ks - pA)) / np.sqrt(pA.shape[1])
    return rmsd, ks

def minimize_rmsd(
    crystA: Cryst, crystB: Cryst, slm: Union[SLM, NDArray[np.int32]], accurate: bool = True
) -> tuple[float, NDArray[np.int32], NDArray[np.int32], NDArray[np.float64]]:
    """Minimize the RMSD with fixed SLM, variable permutation, overall and lattice-vector translations.

    Parameters
    ----------
    crystA : cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : cryst
        The final crystal structure, usually obtained by `load_poscar`.
    slm : slm
        Triplets of integer matrices like `(hA, hB, q)`, representing a SLM.
    accurate: bool
        Whether to use the accurate method to minimize the RMSD. If `False`, the result may be suboptimal. Default is `True`.

    Returns
    -------
    rmsd : float
        The minimum atomic displacement (RMSD) between $(S^{\\text{T}} S)^{1/4}\mathcal{A}$ and $(S^{\\text{T}} S)^{-1/4}\mathcal{B}$.
    p : (Z, ) array of ints
        The permutation of the atoms in `pB_sup` that minimizes the RMSD.
    ks : (3, Z) array of ints
        The lattice-vector translations (fractional coordinates) of the atoms in `pB_sup[:,p]` that minimizes the RMSD.
    t0 : (3, ) array of floats
        The overall translation on `pB` that minimizes the RMSD.
    """
    crystA_sup, crystB_sup, c_sup_half, f_translate = create_common_supercell(crystA, crystB, slm)
    species = crystA_sup[1]
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    t0 = basinhopping(lambda z: minimize_rmsd_t(c_sup_half, species, pA_sup, pB_sup + f_translate @ z.reshape(3,1))[0],
            [0.5, 0.5, 0.5], T=0.05, niter=75, niter_success = None if accurate else 15, minimizer_kwargs={"method": "BFGS"}).x
    _, p, ks = minimize_rmsd_t(c_sup_half, species, pA_sup, pB_sup + f_translate @ t0.reshape(3,1))
    t0 = (pA_sup - pB_sup[:,p] - ks).mean(axis=1, keepdims=True)
    rmsd = la.norm(c_sup_half @ (pA_sup - pB_sup[:,p] - ks - t0)) / np.sqrt(pA_sup.shape[1])
    return rmsd, p, ks, t0