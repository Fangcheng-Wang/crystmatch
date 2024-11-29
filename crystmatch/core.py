"""
The core of crystmatch, including enumeration and optimization algorithms.
"""

from .io import *
from copy import deepcopy
from time import time
from spglib import get_symmetry
from scipy.optimize import brentq, linear_sum_assignment, brute
from scipy.stats.qmc import Sobol
from scipy.spatial.transform import Rotation

np.set_printoptions(suppress=True)
Cryst = Tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]

def get_pure_rotation(cryst: Cryst, symprec: float = 1e-5) -> NDArray[np.int32]:
    """Find all pure rotations appeared in the space group of `cryst`.

    Parameters
    ----------
    cryst : 3-tuple
        `(lattice, species, positions)`, representing the crystal structure, usually obtained by `load_poscar`.
    
    Returns
    -------
    g : (..., 3, 3) array of ints
        A point group of the first kind, containing all pure rotations appeared in the space group of `cryst`, \
            elements of which are integer matrices (under fractional coordinates).

    Examples
    --------
        >>> cryst1 = load_poscar('graphite.txt')
        >>> g1 = get_pure_rotation(cryst1)
    """
    species = cryst[1]
    n = len(species)
    numbers = np.zeros(n, dtype=int)
    temp_s = ''
    temp_n = 0
    for i in range(n):
        if species[i] != temp_s:
            temp_n = temp_n + 1
            temp_s = species[i]
        numbers[i] = temp_n
    g = get_symmetry((cryst[0],cryst[2],numbers), symprec=symprec)['rotations']
    g = g[la.det(g).round(decimals=4)==1,:,:]
    g = np.unique(g, axis=0)
    return g

def check_chem_comp(speciesA, speciesB):
    spA, ctA = np.unique(speciesA, return_counts=True)
    spB, ctB = np.unique(speciesB, return_counts=True)
    assert (spA == spB).all()
    assert np.dot(ctA, ctA) * np.dot(ctB, ctB) == np.dot(ctA, ctB) ** 2
    return

def root_mean_square_strain(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Root-mean-square strain of given singular values.

    Parameters
    ----------
    x : (..., 3) array
        The singular values of 3*3 matrices.
    
    Returns
    -------
    rms_strain : (...) array
        Root-mean-square of `x - 1`.

    Examples
    --------
        >>> x = np.array([[1.1, 0.9, 0.9], [1.2, 1.1, 0.9]])
        >>> root_mean_square_strain(x)
        [0.1        0.14142136]
    """
    return np.sqrt(np.mean((x - 1) ** 2, axis=-1))

def int_fact(n: int) -> List[Tuple[int, int]]:
    """Factorize positive integer `n` into products of two integers.

    Parameters
    ----------
    n : int
        The integer to be factorized.
    
    Returns
    -------
    l : list of 2-tuples of ints
        Contains all `(a, b)` such that a*b=n.
    
    Examples
    --------
        >>> int_fact(6)
        [(1, 6), (2, 3), (3, 2), (6, 1)]
    """
    l = []
    for a in range(1,n+1):
        if n % a == 0: l.append((a, n//a))
    return l

def hnf_list(det: int) -> NDArray[np.int32]:
    """Enumerate all 3*3 column Hermite normal forms (HNFs) with given determinant.

    Parameters
    ----------
    det : int
        The determinant of HNFs.
    
    Returns
    -------
    l : (..., 3, 3) array of ints
        Contains all HNFs with determinant `det`.
    
    Examples
    --------
        >>> hnf_list(2).shape
        (7, 3, 3)
    """
    # Enumerate 3-factorizations of `det`.
    diag_list = []
    for a, aa in int_fact(det):
        for b, c in int_fact(aa):
            diag_list.append((a, b, c))
    # Enumerate HNFs.
    l = []
    for diag in diag_list:
        for h21 in range(diag[1]):
            for h31 in range(diag[2]):
                for h32 in range(diag[2]):
                    h = np.diag(diag)
                    h[1,0] = h21
                    h[2,0] = h31
                    h[2,1] = h32
                    l.append(h)
    l = np.array(l, dtype=int)
    return l

def hnf_decomposition(m: NDArray[np.int32]) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Decompose square integer matrix `m` into product of HNF matrix `h` and unimodular matrix `q`.

    Parameters
    ----------
    m : (N, N) array of ints
        The integer matrix to decompose, with positive determinant.
    
    Returns
    -------
    h : (N, N) array of ints
        The column-style Hermite normal form of `m`.
    q : (N, N) array of ints
        The unimodular matrix satisfying `m` = `h @ q`.
    
    Examples
    --------
        ?
    """
    assert m.dtype == np.int32 and la.det(m) > 0
    N = m.shape[0]
    h = deepcopy(m)
    for i in range(N):
        while not (h[i,i+1:] == 0).all():
            col_nonzero = i + np.nonzero(h[i,i:])[0]
            i0 = col_nonzero[np.argpartition(np.abs(h[i,col_nonzero]), kth=0)[0]]
            h[:,[i,i0]] = h[:,[i0,i]]
            if h[i,i] < 0: h[:,i] = - h[:,i]
            h[:,i+1:] = h[:,i+1:] - np.outer(h[:,i], (h[i,i+1:] / h[i,i]).round().astype(int))
        if h[i,i] < 0: h[:,i] = - h[:,i]
        h[:,:i] = h[:,:i] - np.outer(h[:,i], h[i,:i] // h[i,i])
    q = (la.inv(h) @ m).round().astype(int)
    return h, q

def hnf_rational(m: ArrayLike, max_divisor = 10000) -> NDArray[np.float64]:
    """The Hermite normal form (HNF) of full-row-rank rational matrix `m` (not necessarily square or integer).
    
    Parameters
    ----------
    m : (M, N) array_like, M <= N
        The full-row-rank rational matrix to reduce.
    max_divisor : int
        A positive integer. The least common multiple of all divisors in `m` should not be greater than `max_divisor`.
    
    Returns
    -------
    h : (M, N) array
        The HNF of `m` obtained via elementary column operations over integers.
    
    Examples
    --------
        ?
    """
    for divisor in range(1, max_divisor+1):
        if (np.absolute(np.rint(m * divisor) - m * divisor) <= 1 / max(10000,max_divisor)).all(): break
        elif divisor == max_divisor: print('Warning: input matrix should be rational!')
    h = (m * divisor).round().astype(int)
    M, N = h.shape
    assert M <= N and la.matrix_rank(h, tol=1/max(10000,max_divisor)) == M
    for i in range(M):
        while not (h[i,i+1:] == 0).all():
            col_nonzero = i + np.nonzero(h[i,i:])[0]
            i0 = col_nonzero[np.argpartition(np.abs(h[i,col_nonzero]), kth=0)[0]]
            h[:,[i,i0]] = h[:,[i0,i]]
            if h[i,i] < 0: h[:,i] = - h[:,i]
            h[:,i+1:] = h[:,i+1:] - np.outer(h[:,i], (h[i,i+1:] / h[i,i]).round().astype(int))
        if h[i,i] < 0: h[:,i] = - h[:,i]
        h[:,:i] = h[:,:i] - np.outer(h[:,i], h[i,:i] // h[i,i])
    return h / divisor

def matrix_gcd(m1: ArrayLike, m2: ArrayLike, max_divisor = 10000) -> NDArray[np.float64]:
    """Return a greatest common divisor of rational matrices `m1` and `m2`.
    
    Parameters
    ----------
    m1, m2 : (3, 3) array_like
        Nonsingular rational matrices.
    max_divisor : int
        A positive integer. The least common multiple of all divisors in `m` should not be greater than `max_divisor`.
    
    Returns
    -------
    d : (3, 3) array
        The greatest common divisor of `m1` and `m2` in Hermite normal form.
    """
    assert (la.det([m1, m2]) != 0).all()
    d = hnf_rational(np.hstack((m1, m2)))[:,:3]
    if m1.dtype == np.int32 and m2.dtype == np.int32: return d.round().astype(int)
    return d

def matrix_lcm(m1: ArrayLike, m2: ArrayLike) -> NDArray[np.int32]:
    """Return a least common multiple of integer matrices `m1` and `m2`.
    
    Parameters
    ----------
    m1, m2 : (3, 3) array_like
        Nonsingular integer matrices.
    
    Returns
    -------
    m : (3, 3) array
        The least common multiple of `m1` and `m2` in Hermite normal form.
    """
    assert m1.dtype == np.int32 and m2.dtype == np.int32
    assert (la.det([m1, m2]) != 0).all()
    h = hnf_rational(np.hstack((la.inv(m1.T), la.inv(m2.T))))[:,:3]
    m = la.inv(h.T).round().astype(int)
    return m

def equiv_class_representative(s: Union[SLM, NDArray[np.int32]], gA: NDArray[np.int32], gB: NDArray[np.int32]) -> tuple[SLM, int]:
    """The representative of the equivalence class of `s`.

    Parameters
    ----------
    s : 3-tuple of (3, 3) arrays of ints
        `(hA, hB, q)`, representing a SLM.
    gA : (..., 3, 3) array of ints
        The rotation group of the initial crystal structure, whose elements are \
            integer matrices under fractional coordinates.
    gB : (..., 3, 3) array of ints
        The rotation group of the final crystal structure, whose elements are \
            integer matrices under fractional coordinates.

    Returns
    -------
    ss : 3-tuple of (3, 3) arrays of ints
        The representative of the equivalence class of `s`.
    len_cl : int
        The size of the equivalence class of `s`.
    """
    hA, hB, q = s
    cl = np.transpose(np.dot((gB @ hB) @ q, la.inv(gA @ hA)), axes=[2,0,1,3]).reshape(-1,9)
    cl, i = np.unique(cl.round(decimals=4), axis=0, return_index=True)
    iA, iB = np.unravel_index(i[0], (gA.shape[0], gB.shape[0]))
    hAA, qA = hnf_decomposition(gA[iA] @ hA)
    hBB, qB = hnf_decomposition(gB[iB] @ hB)
    ss = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return ss, len(cl)

def enumerate_slm(
    crystA: Cryst, crystB: Cryst, mu: int, kappa_max: float,
    kappa: Callable[[NDArray[np.float64]], NDArray[np.float64]] = root_mean_square_strain,
    likelihood_ratio: float = 1e2, max_power: int = 1, print_detail: int = 0
) -> List[SLM]:
    """Enumerating all SLMs of multiplicity `mu` with `kappa` smaller than `kappa_max`.

    Parameters
    ----------
    crystA : 3-tuple
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : 3-tuple
        The final crystal structure, usually obtained by `load_poscar`.
    kappa_max : float
        A positive threshold value of `kappa` to determine the range of singular values to generate.
    kappa : callable, optional
        A function that quantifies the strain of a matrix according to its singular values. \
            By default, kappa([x1, x2, x3]) = sqrt(((x1-1)^2 + (x2-1)^2 + (x3-1)^2) / 3).
    likelihood_ratio : float, optional
        The expected likelihood ratio of the enumeration being complete and incomplete. Default is 1e8.
    max_power : int, optional
        Use no more than 2^`max_power` prototypes each time to prevent out of memory. Default is 1.
    print_detail : int, optional
        The level of detail of printing. 0 means no print, 1 means basic settings and 2 means iteration details. \
            Default is 0.


    Returns
    -------
    slist : list of 3-tuples of (3, 3) arrays of ints
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.
    
    Examples
    --------
        >>> slist = enumerate_slm(austenite, martensite, 6, 0.15)
        >>> slist[0]
        ?
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
    gA = get_pure_rotation(crystA)
    gB = get_pure_rotation(crystB)
    max_prob_ratio = 20
    # Compute the sampling domain of `s0`s.
    det_s = zA / zB * la.det(cB) / la.det(cA)
    diff_kappa = lambda x: kappa(np.array([x, (det_s / x) ** 0.5, (det_s / x) ** 0.5])) - kappa_max
    a = brentq(diff_kappa, 0.1, det_s**(1/3))
    b = brentq(diff_kappa, det_s**(1/3), 10)
    max_strain = max(abs(a-1), abs(b-1))
    # Enumerate SLMs.
    print(f"\nEnumerating SLMs (mu = {mu:.0f}, {kappa.__name__} <= {kappa_max:.4f}) ...")
    if print_detail >= 1:
        print(f"\tprototype sampling domain: SO(3) (±{max_strain:.2f} for each matrix element)")
        print(f"\tassumed maximum probability ratio among classes: {max_prob_ratio:.1f}")
        print(f"\texpected likelihood ratio: {likelihood_ratio:.1f}")
        print("\tnum_s0\tm\tm*\telapsed_time(s)")
    slist = []
    iter = 0
    num_s0 = 0
    m = 0
    t = time()
    sobol_seq = Sobol(12)
    while m <= (1 + len(slist) * max_prob_ratio) * np.log(likelihood_ratio):    # or True for Debug!
        # Sampling `s0`s around SO(3).
        if iter == 0: rand_num = sobol_seq.random_base2(iter + 1)
        elif iter <= max_power: rand_num = sobol_seq.random_base2(iter)
        else: rand_num = sobol_seq.random(2 ** max_power)
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
        # Check strains of `s`s.
        index2 = np.nonzero(kappa(la.svd(cB @ hB[index1[2]] @ q[index1] @ hA_inv[index1[1]] @ la.inv(cA), compute_uv=False)) < kappa_max)[0]                                          
        for i in index2:
            s = (hA[index1[1][i]], hB[index1[2][i]], q[index1[0][i], index1[1][i], index1[2][i]])
            s, _ = equiv_class_representative(s, gA, gB)
            repeated = False
            for j in range(len(slist)):
                if (s[0] == slist[j][0]).all() and (s[1] == slist[j][1]).all() and (s[2] == slist[j][2]).all():
                    # `s` is repeated.
                    m = m + 1
                    repeated = True
                    break
            if not repeated:
                # `s` is new.
                slist.append(s)
                m = 0
        num_s0 += s0.shape[0]
        if print_detail >= 2:
            print(f"\t{num_s0}\t{m}\t{(1 + len(slist) * max_prob_ratio) * np.log(likelihood_ratio):.0f}\t{time()-t:.2f}")
        iter = iter + 1
        if iter > 40 and len(slist) == 0: break
    print(f'{len(slist):.0f} inequivalent SLMs (in {time()-t:.2f} seconds)')
    return slist

def complete_slm_list(
    crystA: Cryst, crystB: Cryst, mu_max: int, kappa_max: float,
    kappa: Callable[[NDArray[np.float64]], NDArray[np.float64]] = root_mean_square_strain
) -> NDArray[np.int32]:
    """Enumerating all SLMs with multiplicity not bigger than `mu` and `kappa` not bigger than `kappa_max`.

    Parameters
    ----------
    crystA : 3-tuple
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : 3-tuple
        The final crystal structure, usually obtained by `load_poscar`.
    mu_max : int
        The maximum value of `mu`.
    kappa_max : float
        A positive threshold value of `kappa` to determine the range of singular values to generate.
    kappa : callable, optional
        A function that quantifies the strain of a matrix according to its singular values. \
            By default, kappa([x1, x2, x3]) = sqrt(((x1-1)^2 + (x2-1)^2 + (x3-1)^2) / 3).

    Returns
    -------
    slist : (..., 3, 3, 3) array of ints
        Contains triplets of integer matrices like `(hA, hB, q)`, representing inequivalent SLMs.
    
    Examples
    --------
        >>> slist = complete_slm_list(austenite, martensite, 6, 0.15)
        >>> slist[-1]
        ?
    """
    slist = []
    for mu in range(1,mu_max+1):
        slist = slist + enumerate_slm(crystA, crystB, mu, kappa_max, kappa=kappa)
    slist = np.array(slist)
    hA = slist[:,0,:,:]
    hB = slist[:,1,:,:]
    q = slist[:,2,:,:]
    corelist = (hB @ q @ la.inv(hA)).round(decimals=4)
    _, ind_unique = np.unique(corelist, axis=0, return_index=True)
    print(f"\nA total of {slist.shape[0]:.0f} SLMs ({ind_unique.shape[0]:.0f} deformation gradients) are enumerated.")
    return slist[ind_unique]

def vector_reduce(v: NDArray, divisors: NDArray) -> NDArray:
    """Minimizing the norm of `v` by adding and subtracting columns of `divisors`.

    Parameters
    ----------
    v : (N,) array
        The vector to be reduced.
    divisors : (N, ...) array
        The vectors used to translate `v`.
    
    Returns
    -------
    vv : (N,) array
        The reduced `v` with minimum Euclidean norm.

    Examples
    --------
        ?
    """
    vv = deepcopy(v)
    converged = False
    while not converged:
        converged = True
        for i in range(divisors.shape[1]):
            v0 = divisors[:,i]
            while la.norm(vv + v0) < la.norm(vv):
                converged = False
                vv = vv + v0
            while la.norm(vv - v0) < la.norm(vv):
                converged = False
                vv = vv - v0
    return vv

def niggli_cell(c: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """Reduce cell `c` to its Niggli cell.

    Parameters
    ----------
    c : (3, 3) array
        The cell to be reduced, whose columns are cell vectors.
    
    Returns
    -------
    cc : (3, 3) array
        The Niggli cell, with shortest right-handed cell vectors.
    q : (3, 3) array of ints
        The unimodular matrix satisfying `cc = c @ q`.

    Examples
    --------
        ?
    """
    c0 = np.zeros((3,3))
    cc = deepcopy(c)
    while (cc != c0).any():
        c0 = deepcopy(cc)
        cc = cc[:,np.argsort(la.norm(cc, axis=0))]
        cc[:,2] = vector_reduce(cc[:,2], cc[:,0:2])
    if la.det(cc) < 0: cc = -cc
    q = (la.inv(c) @ cc).round().astype(int)
    return cc, q

def int_vectors_inside(c: NDArray[np.int32]) -> NDArray[np.int32]:
    """Integer vectors inside the cell `c` whose elements are integers.

    Parameters
    ----------
    c : (3, 3) array of ints
        A matrix whose columns are integer cell vectors.
    
    Returns
    -------
    v_int : (3, ...) array of ints
        Its columns are vectors satisfying `v = c @ k`, where 0 <= `k[0]`, `k[1]`, `k[2]` < 1.
    
    Examples
    --------
        ?
    """
    assert c.dtype == np.int32
    vertices = c @ np.mgrid[0:2,0:2,0:2].reshape(3,-1)
    candidates = np.mgrid[np.amin(vertices[0,:]):np.amax(vertices[0,:])+1, np.amin(vertices[1,:]):np.amax(vertices[1,:])+1, \
        np.amin(vertices[2,:]):np.amax(vertices[2,:])+1].reshape(3,-1)
    fractional = (la.inv(c) @ candidates).round(decimals=7)
    is_inside = (np.prod(fractional < 1, axis=0) * np.prod(fractional >= 0, axis=0)).astype(bool)
    assert np.sum(is_inside) == la.det(c).round().astype(int)
    return candidates[:,is_inside]

def local_minimum_rmsd(
    c: NDArray[np.float64], species: NDArray[np.str_], pA: NDArray[np.float64], pB: NDArray[np.float64]
) -> tuple[float, NDArray[np.float64]]:
    """
    The local minimum root-mean-square distance (RMSD) between two crystal structures with same cells.

    Parameters
    ----------
    c : (3, 3) array
        The matrix whose columns are cell vectors of both crystal structures.
    species : (N,) array of strs
        The array that specifies the type of each ion.
    pA, pB : (3, N) array
        The matrices whose columns are fractional coordinates of the ions in the crystal structures.
        
    Returns
    -------
    rmsd : float
        The minimized RMSD (calculated under periodic boundary condition).
    pB_m : (3, N) array
        Modified `pB` with ions permuted, collectively translated, and respectively shifted by lattice vectors.
    
    Examples
    --------
        ?
    """
    assert pA.shape[1] == pB.shape[1]
    # Determine best lattice-vector shift.
    shift = np.mgrid[-1:2, -1:2, -1:2].reshape(3,1,1,-1)
    relative_position = np.tensordot(c, (pA.reshape(3,-1,1,1) - pB.reshape(3,1,-1,1) - shift), axes=[1,0])
    d2_tensor = np.sum(relative_position ** 2, axis=0)
    index_shift = np.argmin(d2_tensor, axis=2)
    d2_matrix = np.amin(d2_tensor, axis=2)
    # Determine best atomic correspondence.
    index_B = np.arange(pB.shape[1], dtype=int)
    for s in np.unique(species):
        is_s = species == s
        index_B[is_s] = index_B[is_s][linear_sum_assignment(d2_matrix[is_s][:,is_s])[1]]
    assert (np.sort(index_B) == np.arange(pB.shape[1])).all()
    pB_m = pB[:,index_B] + shift[:,0,0,index_shift[np.arange(pA.shape[1]),index_B]]
    # Determine translation vector for minimum RMSD under above shift and correspondence
    pB_m = pB_m - np.mean(pB_m - pA, axis=1).reshape(3,1)
    rmsd = la.norm(c @ (pB_m - pA)) / np.sqrt(pA.shape[1])
    return rmsd, pB_m

def minimize_rmsd(
    crystA: Cryst, crystB: Cryst, s: Union[SLM, NDArray[np.int32]], n_grid: int = 5
) -> tuple[float, Cryst, Cryst, Cryst]:
    """
    Minimize the RMSD by optimizing the translation and atomic correspondence compatible with SLM `s`.

    Parameters
    ----------
    crystA : 3-tuple
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : 3-tuple
        The final crystal structure, usually obtained by `load_poscar`.
    s : (3, 3, 3) array_like of ints
        Triplets of integer matrices like `(hA, hB, q)`, representing a SLM.
    n_grid : int, optional
        The number of grid points for translation. Default is 5.

    Returns
    -------
    rmsd : float
        The minimum atomic displacement (RMSD) between distorted `crystA_sup` and `crystB_sup`.
    crystA_sup : 3-tuple
        The crystal structure identical with `crystA`, but described by supercell `cA @ hA`.
    crystB_sup : 3-tuple
        The crystal structure identical with `crystB`, transformed from `crystA` along the path of minimum `distance`. One may \
            generate a linear transition path using `cryst_interpolation(crystA_sup, crystB_sup, num)`.
    cryst_com : 3-tuple
        ?
    
    Examples
    --------
        ?
    """
    # Create distorted `crystA_sup` and `crystB_sup` that share a common periodicity.
    cA = crystA[0].T
    cB = crystB[0].T
    speciesA = crystA[1]
    speciesB = crystB[1]
    pA = crystA[2].T
    pB = crystB[2].T
    check_chem_comp(speciesA, speciesB)
    hA, hB, q = s
    s_mat = cB @ hB @ q @ la.inv(cA @ hA)
    u, sigma, vT = la.svd(s_mat)
    c_com, q_com = niggli_cell(vT.T @ np.diag(sigma ** 0.5) @ vT @ cA @ hA)         # The common supercell, half-distorted.
    mA = hA @ q_com
    mB = hB @ q @ q_com
    cA_sup = cA @ mA
    cB_sup = (u @ vT).T @ cB @ mB           # The natural orientation of `crystB`, such that cB_sup = (v @ sigma @ v.T) @ cA_sup.
    speciesA_sup = np.tile(speciesA, la.det(mA).round().astype(int))
    speciesB_sup = np.tile(speciesB, la.det(mB).round().astype(int))
    pA_sup = la.inv(mA) @ (pA.reshape(3,1,-1) + int_vectors_inside(mA).reshape(3,-1,1)).reshape(3,-1)
    pB_sup = la.inv(mB) @ (pB.reshape(3,1,-1) + int_vectors_inside(mB).reshape(3,-1,1)).reshape(3,-1)
    # Collect each species of ions.
    argsortA = np.argsort(speciesA_sup)
    argsortB = np.argsort(speciesB_sup)
    species_com = speciesA_sup[argsortA]
    assert (species_com == speciesB_sup[argsortB]).all()
    pA_sup = pA_sup[:,argsortA]
    pB_sup = pB_sup[:,argsortB]
    # Calculate the irreducible translation domain.
    tA = la.inv(mA) @ int_vectors_inside(mA)        # Translation elements of `crystA` (in fractional coordinates)
    tB = la.inv(mB) @ int_vectors_inside(mB)        # Translation elements of `crystB` (in fractional coordinates)
    t_com = np.hstack(((tA.reshape(3,-1,1) - tB.reshape(3,1,-1)).reshape(3,-1), np.eye(3)))         # Common translation elements
    t_com = t_com[:,np.argsort(la.norm(t_com, axis=0))]
    t_cell = []
    for j in range(t_com.shape[1]):
        tt = t_com[:,j]
        if len(t_cell) == 0:
            if (tt.round(decimals=4) != 0).any(): t_cell.append(tt)
        elif len(t_cell) == 1:
            if (np.cross(tt, t_cell[0]).round(decimals=4) != 0).any(): t_cell.append(tt)
        elif len(t_cell) == 2:
            if la.det(np.array(t_cell + [tt])).round(decimals=4) != 0:
                t_cell.append(tt)
                break
    t_cell, _ = niggli_cell(np.array(t_cell).T)      # The irreducible translation domain cell.
    # Optimize the translation vector and the ion order in `crystB`.
    z = brute(lambda z: local_minimum_rmsd(c_com, species_com, pA_sup, pB_sup + t_cell @ z.reshape(3,1))[0], ((0,1),(0,1),(0,1)), Ns=n_grid, finish=False)
    rmsd, pB_sup_m = local_minimum_rmsd(c_com, species_com, pA_sup, pB_sup + t_cell @ z.reshape(3,1))
    crystA_sup = (cA_sup.T, species_com, pA_sup.T)
    crystB_sup = (cB_sup.T, species_com, pB_sup_m.T)
    n = species_com.shape[0]
    cryst_com = (c_com.T, np.concatenate((np.char.add(species_com, ['A'] * n), np.char.add(species_com, ['B'] * n))), np.hstack((pA_sup,pB_sup_m)).T)
    return rmsd, crystA_sup, crystB_sup, cryst_com