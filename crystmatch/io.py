"""
Loading and saving files for crystmatch.
"""

import numpy as np
import numpy.linalg as la
from spglib import get_spacegroup, find_primitive
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
    cryst : 3-tuple
        `(lattice, species, positions)`, representing a crystal structure.
    lattice : (3, 3) array
        In accord with the POSCAR format and the scaling factor is multiplied.
    species : (N,) array of strs
        Respective species of the ions.
    positions : (N, 3) array
        In accord with the POSCAR format ('Direct' mode).

    Examples
    --------
        >>> cryst1 = load_poscar('graphite.txt')
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
    cryst : 3-tuple
        `(lattice, species, positions)`, representing a crystal structure, usually obtained by `load_poscar` or `minimize_rmsd`.
    crystname : str, optional
        A system description to write to the comment line of the POSCAR file. If `crystname = None`, `filename` will be used.
    
    Examples
    --------
        >>> save_poscar('mycryst.txt', mycryst)
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

def load_from_dat(filename: str, index: int) -> Cryst:
    """
    Load the initial structure, final structure, and CSM from a DAT file.

    Parameters
    ----------
    filename : str
        The name of the DAT file to be read.
    index : int
        The index of the CSM to be loaded.

    Returns
    -------
    crystA : 3-tuple
        `(lattice, species, positions)`, representing a crystal structure, usually obtained by `load_poscar` or `minimize_rmsd`.
    crystB : 3-tuple
        `(lattice, species, positions)`, representing a crystal structure, usually obtained by `load_poscar` or `minimize_rmsd`.
    csm : 3-tuple
        `(slm, permutation, translation)`, representing a CSM between `crystA` and `crystB`.
    """
    pass

def csm_to_str(filename: str, crystA: Cryst, crystB: Cryst, csm: SLM) -> str:
    pass

def cryst_interpolation(crystA_sup: Cryst, crystB_sup: Cryst, images: int) -> List[Cryst]:
    """
    Linear interpolation between `crystA` and `crystB`.

    Parameters
    ----------
    crystA_sup, crystB_sup : 3-tuple
        The initial and final crystal structures with specified atomic correspondence, usually obtained by `minimize_rmsd`.
    images : int
        Number of interpolations to generate, must be non-negative.
    
    Returns
    -------
    crystlist : list of 3-tuples
        Contains crystal structures generated from interpolation.
    """
    species = crystA_sup[1]
    assert (species == crystB_sup[1]).all()
    cA = crystA_sup[0].T
    cB = crystB_sup[0].T
    pA = crystA_sup[2].T
    pB = crystB_sup[2].T
    s = cB @ la.inv(cA)
    if ((s.T - s).round(decimals=4) != 0).any(): print('\nWarning: Extra rotation detected, which is removed in output.')
    _, sigma, vT = la.svd(s)
    crystlist = []
    tlist = np.linspace(0, 1, images)
    for t in tlist:
        c = vT.T @ np.diag(sigma ** t) @ vT @ cA
        p = pA * (1-t) + pB * t
        crystlist.append((c.T, species, p.T))
    return crystlist

def save_trajectory(filename: str, crystA_sup: Cryst, crystB_sup: Cryst, images: int = 50, crystname: Union[str, None] = None) -> None:
    """
    Save the linear interpolation between `crystA` and `crystB` to a single XDATCAR file.
    
    Parameters
    ----------
    filename : str
        The name of the file to save, must not already exist in current directory.
    crystA_sup, crystB_sup : 3-tuple
        The initial and final crystal structures with specified atomic correspondence, usually obtained by `minimize_rmsd`.
    images : int, optional
        Number of interpolations to generate, must be non-negative. Default is 50.
    crystname : str, optional
        A system description to write to the comment line of the POSCAR file. If `crystname = None`, `filename` will be used.
    
    Examples
    --------
        >>> save_trajectory('mytrajectory.txt', mycrystA, mycrystB)
    """
    crystlist = cryst_interpolation(crystA_sup, crystB_sup, images)
    f = open(filename, mode='x')
    for i in range(images):
        if i > 0: f.write('\n')
        if crystname: f.write(crystname)
        else: f.write(filename.split(sep='.')[0])
        lattice = crystlist[i][0]
        species = crystlist[i][1]
        positions = crystlist[i][2]
        f.write('\n1.0\n')
        f.write('\n'.join(f'{v[0]:.12f}\t{v[1]:.12f}\t{v[2]:.12f}' for v in lattice.tolist()))
        species_unique, species_counts = np.unique(species, return_counts=True)
        f.write('\n' + ' '.join(species_unique.tolist()))
        f.write('\n' + ' '.join(str(n) for n in species_counts.tolist()))
        f.write(f'\nDirect configuration= {i+1:.0f}\n')
        f.write('\n'.join(f'{p[0]:.12f}\t{p[1]:.12f}\t{p[2]:.12f}' for p in positions.tolist()))
    f.close()
    return