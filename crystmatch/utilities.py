"""
Load/save crystal structures and CSMs from/to files.
"""

import numpy as np
import numpy.linalg as la
from copy import deepcopy
from spglib import get_spacegroup, find_primitive, get_symmetry
from typing import Union, Tuple, List, Callable
from numpy.typing import NDArray, ArrayLike

np.set_printoptions(suppress=True)
Cryst = Tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]

def load_poscar(filename: str, to_primitive: bool = True, symprec: float = 1e-5) -> Cryst:
    """Load the crystal structure from a POSCAR file.

    Parameters
    ----------
    filename : str
        The name of the POSCAR file to be read.
    to_primitive : bool, optional
        Using the primitive cell instead of the cell given by the POSCAR file. Default is True.

    Returns
    -------
    cryst : cryst
        The loaded crystal structure, consisting of the lattice vectors, species, and positions.
    lattice : (3, 3) array
        In accord with the POSCAR format and the scaling factor is multiplied.
    species : (N,) array of strs
        Respective species of the ions.
    positions : (N, 3) array
        In accord with the POSCAR format ('Direct' mode).
    """
    f = open(filename, mode='r')
    print(f"Loading crystal structure '{f.readline()[:-1]}' from file '{filename}'...")
    a = np.array(f.readline()[:-1], dtype=float)
    lattice = np.zeros((3,3), dtype=float)
    for i in range(3):
        lattice[i,:] = np.array(f.readline().split(), dtype=float)
    if la.det(lattice) < 0: lattice[2,:] = - lattice[2,:]
    lattice = a * lattice
    species_name = f.readline().split()
    species_counts = np.array(f.readline().split(), dtype=int)
    species = []
    for i in range(len(species_counts)):
        species = species + [species_name[i]] * species_counts[i]
    species = np.array(species, dtype=str)
    unit = ''
    while not unit in ['D','d','C','c','K','k']:
        unit = f.readline()[0]
    N = species_counts.sum()
    positions = np.zeros((N,3), dtype=float)
    if unit in ['D','d']:
        for i in range(N):
            positions[i,:] = np.array(f.readline().split()[:3], dtype=float)
    elif unit in ['C','c','K','k']:
        for i in range(N):
            positions[i,:] = np.dot(la.inv(lattice.transpose()), np.array(f.readline().split()[:3], dtype=float))
    f.close()
    _, numbers = np.unique(species, return_inverse=True)
    print(f"\tSpace group: {get_spacegroup((lattice, positions, numbers), symprec=symprec)}")
    if to_primitive:
        lattice, positions, numbers = find_primitive((lattice, positions, numbers))
        if len(numbers) != len(species):
            print(f"\tCell in POSCAR file is not primitive! Using primitive cell (Z = {len(numbers):d}) now.")
        else:
            print(f"\tCell in POSCAR file is already primitive (Z = {len(numbers):d}).")
        species = np.array(species_name)[numbers]
    else: print(f"\tUsing cell in POSCAR file (Z = {len(numbers):d}).")
    cryst = (lattice, species, positions)
    return cryst

def save_poscar(filename: str, cryst: Cryst, crystname: Union[str, None] = None) -> None:
    """
    Save the crystal structure to a POSCAR file.

    Parameters
    ----------
    filename : str
        The name of the file to save, must not already exist in current directory.
    cryst : cryst
        The crystal structure to be saved, consisting of the lattice vectors, species, and positions.
    crystname : str, optional
        A system description to write to the comment line of the POSCAR file. If `crystname = None`, `filename` will be used.
    """
    f = open(filename, mode='x')
    if crystname: f.write(crystname)
    else: f.write(filename.split(sep='.')[0])
    lattice = cryst[0]
    species = cryst[1]
    positions = cryst[2]
    f.write('\n1.0\n')
    f.write('\n'.join(f'{v[0]:.12f}\t{v[1]:.12f}\t{v[2]:.12f}' for v in lattice.tolist()))
    species_unique, species_counts = np.unique(species, return_counts=True)
    f.write('\n' + ' '.join(species_unique.tolist()))
    f.write('\n' + ' '.join(str(n) for n in species_counts.tolist()))
    f.write('\nDirect\n')
    f.write('\n'.join(f'{p[0]:.12f}\t{p[1]:.12f}\t{p[2]:.12f}' for p in positions.tolist()))
    f.close()
    return

def check_chem_comp(speciesA, speciesB):
    spA, ctA = np.unique(speciesA, return_counts=True)
    spB, ctB = np.unique(speciesB, return_counts=True)
    assert (spA == spB).all()
    assert np.dot(ctA, ctA) * np.dot(ctB, ctB) == np.dot(ctA, ctB) ** 2
    return

def create_common_supercell(crystA: Cryst, crystB: Cryst, slm: SLM) -> Tuple[Cryst, Cryst, NDArray[np.float64], NDArray[np.float64]]:
    """Create supercell structures representing A, (S^T S)^(1/2)S^(-1)B. Also return the half-distorted supercell and translation cell.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slm : slm
        The SLM of the CSM.
    
    Returns
    -------
    crystA_sup : cryst
        The supercell of the initial structure.
    crystB_sup_final : cryst
        The supercell of the final structure, but deformed to have the same cell vectors as (S^T S)^(1/2) `crystA_sup` (rotation-free).
    c_sup_half : (3, 3) array of floats
        The half-distorted supercell, i.e., (S^T S)^(1/4) `latticeA_sup` = (S^T S)^(-1/4) `latticeB_sup_final`.
    f_translate : (3, 3) array of floats
        The translation cell (fractional coordinates) of the shuffle, i.e., \
            gcd( (S^T S)^(1/4) `latticeA_sup` , (S^T S)^(-1/4) `latticeB_sup_final` ).
    """
    # Unpacking crystal structures.
    cA = crystA[0].T
    cB = crystB[0].T
    speciesA = crystA[1]
    speciesB = crystB[1]
    pA = crystA[2].T
    pB = crystB[2].T
    check_chem_comp(speciesA, speciesB)
    
    # Determining the supercell geometries from the SLM.
    hA, hB, q = slm
    deform = cB @ hB @ q @ la.inv(cA @ hA)
    u, sigma, vT = la.svd(deform)
    c_sup_half, q_sup = niggli_cell(vT.T @ np.diag(sigma ** 0.5) @ vT @ cA @ hA)         # The half-distorted supercell.
    mA = hA @ q_sup
    mB = hB @ q @ q_sup
    cA_sup = cA @ mA
    cB_sup_final = (u @ vT).T @ cB @ mB                     # The rotation-free orientation of `crystB_sup`.
    
    # Sorting supercell species and positions.
    speciesA_sup = np.tile(speciesA, la.det(mA).round().astype(int))
    speciesB_sup = np.tile(speciesB, la.det(mB).round().astype(int))
    pA_sup = la.inv(mA) @ (pA.reshape(3,1,-1) + int_vec_inside(mA).reshape(3,-1,1)).reshape(3,-1)
    pB_sup = la.inv(mB) @ (pB.reshape(3,1,-1) + int_vec_inside(mB).reshape(3,-1,1)).reshape(3,-1)
    argsortA = np.argsort(speciesA_sup)
    argsortB = np.argsort(speciesB_sup)
    assert (speciesA_sup[argsortA] == speciesB_sup[argsortB]).all(), "Error: Please report this bug to wfc@pku.edu.cn if you see this message."
    species_sup = speciesA_sup[argsortA]
    pA_sup = pA_sup[:,argsortA]
    pB_sup = pB_sup[:,argsortB]
    
    # Computing output.
    crystA_sup = (cA_sup.T, species_sup, pA_sup.T)
    crystB_sup_final = (cB_sup_final.T, species_sup, pB_sup.T)
    f_translate = niggli_cell(matrix_gcd(la.inv(mA), la.inv(mB)))[0]
    return crystA_sup, crystB_sup_final, c_sup_half, f_translate

def int_arrays_to_pair(crystA: Cryst, crystB: Cryst, slm: SLM,
    p: NDArray[np.int32], ks: NDArray[np.int32], centered: bool = True
) -> Tuple[Cryst, Cryst]:
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
    centered : bool, optional
        Whether to make the centers of `crystA_sup` and `crystB_sup_final` coincide. Default is True.
        
    Returns
    -------
    crystA_sup : cryst
        The initial structure, whose lattice vectors and atoms are matched to `crystB_sup` according to the CSM.
    crystB_sup_final : cryst
        The final structure, whose lattice vectors and atoms are matched to `crystA_sup` according to the CSM, with rotation-free orientation.
    """
    crystA_sup, crystB_sup, _, _ = create_common_supercell(crystA, crystB, slm)
    positionsB_sup = crystB_sup[2][p,:] + ks.T
    if centered:
        positionsB_sup = positionsB_sup - np.mean(positionsB_sup - crystA_sup[2], axis=0)
    return crystA_sup, (crystB_sup[0], crystB_sup[1][p], positionsB_sup)

def rms_minus1(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Root-mean-square strain of given singular values.

    Parameters
    ----------
    x : (..., 3) array
        The singular values of 3*3 matrices.
    
    Returns
    -------
    rms_strain : (...) array
        Root-mean-square of `x - 1`.
    """
    return np.sqrt(np.mean((x - 1) ** 2, axis=-1))

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

def int_vec_inside(c: NDArray[np.int32]) -> NDArray[np.int32]:
    """Integer vectors inside the cell `c` whose elements are integers.

    Parameters
    ----------
    c : (3, 3) array of ints
        A matrix whose columns are integer cell vectors.
    
    Returns
    -------
    v_int : (3, ...) array of ints
        Its columns are vectors satisfying `v = c @ k`, where 0 <= `k[0]`, `k[1]`, `k[2]` < 1.
    """
    assert c.dtype == int
    vertices = c @ np.mgrid[0:2,0:2,0:2].reshape(3,-1)
    candidates = np.mgrid[np.amin(vertices[0,:]):np.amax(vertices[0,:])+1, np.amin(vertices[1,:]):np.amax(vertices[1,:])+1, \
        np.amin(vertices[2,:]):np.amax(vertices[2,:])+1].reshape(3,-1)
    fractional = (la.inv(c) @ candidates).round(decimals=7)
    is_inside = (np.prod(fractional < 1, axis=0) * np.prod(fractional >= 0, axis=0)).astype(bool)
    assert np.sum(is_inside) == la.det(c).round().astype(int)
    return candidates[:,is_inside]

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
    assert la.det(m1) != 0 and la.det(m2) != 0
    d = hnf_rational(np.hstack((m1, m2)))[:,:3]
    if m1.dtype == int and m2.dtype == int: d = d.round().astype(int)
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
    assert m1.dtype == int and m2.dtype == int
    assert (la.det([m1, m2]) != 0).all()
    h = hnf_rational(np.hstack((la.inv(m1.T), la.inv(m2.T))))[:,:3]
    m = la.inv(h.T).round().astype(int)
    return m

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

def hnf_int(m: NDArray[np.int32]) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
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
    """
    assert m.dtype == int and la.det(m) > 0
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