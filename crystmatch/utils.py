"""
Load/save crystal structures and CSMs from/to files.
"""

import numpy as np
import numpy.linalg as la
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
    print("\nCrystal structure '" + f.readline()[:-1] + "' loaded")
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
    numbers = np.zeros(N, dtype=int)
    for i in range(1, len(species_counts)):
        numbers[np.sum(species_counts[:i]):] = numbers[np.sum(species_counts[:i]):] + 1
    print('Space group: ' + get_spacegroup((lattice, positions, numbers), symprec=symprec))
    if to_primitive:
        lattice, positions, numbers = find_primitive((lattice, positions, numbers))
        if len(numbers) != len(species):
            print(f"Cell given by the POSCAR file is not primitive. Using primitive cell (Z = {len(numbers):.0f}) now.")
        else:
            print(f"Cell given by POSCAR file is primitive (Z = {len(numbers):.0f}).")
        species = np.array(species_name)[numbers]
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
    """
    assert c.dtype == np.int32
    vertices = c @ np.mgrid[0:2,0:2,0:2].reshape(3,-1)
    candidates = np.mgrid[np.amin(vertices[0,:]):np.amax(vertices[0,:])+1, np.amin(vertices[1,:]):np.amax(vertices[1,:])+1, \
        np.amin(vertices[2,:]):np.amax(vertices[2,:])+1].reshape(3,-1)
    fractional = (la.inv(c) @ candidates).round(decimals=7)
    is_inside = (np.prod(fractional < 1, axis=0) * np.prod(fractional >= 0, axis=0)).astype(bool)
    assert np.sum(is_inside) == la.det(c).round().astype(int)
    return candidates[:,is_inside]

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