# Copyright (C) 2024 Fangcheng Wang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Common functions used in crystmatch."""

from __future__ import annotations
from os import makedirs, environ
from os.path import sep, exists, splitext
from collections import namedtuple
from scipy.spatial.transform import Rotation
from scipy.optimize import linear_sum_assignment, minimize, basinhopping
from scipy.stats.qmc import Sobol
from spglib import get_spacegroup, get_symmetry, standardize_cell, refine_cell, get_symmetry_dataset
from tqdm import tqdm
from typing import Union, Dict, Callable, Literal, Optional
from numpy.typing import NDArray, ArrayLike
from matplotlib import rcParams, colors
import numpy as np
import numba as nb
import numpy.linalg as la
import matplotlib.pyplot as plt
rcParams.update({
    'font.family': 'serif',
    'pgf.rcfonts': False,
    'figure.dpi': 150
})

np.set_printoptions(suppress=True)
Cryst = tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]
Table = namedtuple("Table", ["data", "header"])
NPZ_ARR_COMMENT = "arr_mu.npy saves the PCTs of those CSMs with multiplicity mu"

def load_poscar(
    filename: str,
    to_primitive: bool = True,
    tol: float = 1e-3,
    verbose: bool = True
) -> Cryst:
    """Load the crystal structure from a POSCAR file.

    Parameters
    ----------
    filename : str
        The name of the POSCAR file to be read.
    to_primitive : bool, optional
        If True, reduced the crystal structure to its primitive structure.
    tol : float, optional
        The tolerance for `spglib` symmetry detection.
    verbose : bool, optional
        If True, print verbose output during loading.

    Returns
    -------
    cryst : cryst
        The loaded crystal structure, represented as a tuple of `(lattice, species, positions)`.
    """
    with open(filename, mode='r') as f:
        f.readline()
        if verbose: print(f"Loading crystal structure from file '{filename}' ...")
        a = np.array(f.readline()[:-1], dtype=float)
        lattice = np.zeros((3,3), dtype=float)
        for i in range(3):
            lattice[i,:] = np.array(f.readline().split(), dtype=float)
        if la.det(lattice) < 0: lattice[2,:] = - lattice[2,:]
        lattice = a * lattice
        sp_name = f.readline().split()
        sp_counts = np.array(f.readline().split(), dtype=int)
        species = []
        for i in range(len(sp_name)):
            species = species + [sp_name[i]] * sp_counts[i]
        species = np.array(species, dtype=str)
        unit = ''
        while not unit.startswith(('D','d','C','c','K','k')):
            unit = f.readline().strip()
        z = sp_counts.sum()
        positions = np.zeros((z,3), dtype=float)
        for i in range(z):
            if unit.startswith(('D','d')):
                positions[i,:] = np.array(f.readline().split()[:3], dtype=float)
            elif unit.startswith(('C','c','K','k')):
                positions[i,:] = np.dot(la.inv(lattice.transpose()), np.array(f.readline().split()[:3], dtype=float))
    indsort = np.argsort(species)
    cryst = (lattice, species[indsort], positions[indsort])

    if verbose: print(f"\tSpace group: {get_spacegroup(cryst_to_spglib(cryst), symprec=tol)}.")
    if to_primitive:
        cryst = primitive_cryst(cryst, tol=tol)
        if verbose:
            if len(cryst[1]) != len(species): print(f"\tCell in POSCAR file is not primitive! Using primitive cell (Z = {len(cryst[1]):d}) now.")
            else: print(f"\tCell in POSCAR file is already primitive (Z = {len(cryst[1]):d}).")
    elif verbose: print(f"\tUsing cell in POSCAR file (Z = {len(species):d}).")
    return cryst

def load_csmcar(
    filename: str,
    verbose: bool = True
) -> tuple[NDArray[np.float64], NDArray[np.float64], Dict[str, float], NDArray[np.float64]]:
    """Load the CSMCAR file, which contains `crystmatch` parameters.

    Parameters
    ----------
    filename : str
        The name of the POSCAR file to be read.
    verbose : bool, optional
        If True, print verbose output during loading.

    Returns
    -------
    voigtA, voigtB : (6, 6) array
        The loaded elastic tensor for the initial and final structure, in Voigt notation (ordered as XX, YY, ZZ, YZ, ZX, XY).
    weight_func : dict
        The loaded weight function for the shuffle distance.
    ori_rel : (2, 2, 3) array
        The two loaded parallelisms, representing the orientation relationship between the initial and final structure.
    """
    with open(filename, mode='r') as f:
        if verbose: print(f"Loading crystmatch parameters from file '{filename}' ...")
        voigtA = None
        voigtB = None
        weight_func = None
        ori_rel = None
        standard_axes = ['XX', 'YY', 'ZZ', 'YZ', 'ZX', 'XY']
        nl = '\n'
        tab = '\t'
        info = f.readline()
        while info:
            info = info.split('#')[0].strip()
            if info.startswith(('I', 'i')):
                if voigtA: print("Warning: Initial elastic tensor is defined multiple times! The last definition will be used.")
                voigtA = np.zeros((6, 6))
                axes = []
                for i in range(6):
                    l = f.readline().strip().split()
                    axes.append(l[0])
                    voigtA[i,:] = [float(x) for x in l[1:7]]
                axis_yz = 'YZ' if 'YZ' in axes else 'ZY'
                axis_zx = 'ZX' if 'ZX' in axes else 'XZ'
                axis_xy = 'XY' if 'XY' in axes else 'YX'
                ind = [axes.index('XX'), axes.index('YY'), axes.index('ZZ'), axes.index(axis_yz), axes.index(axis_zx), axes.index(axis_xy)]
                voigtA = voigtA[np.ix_(ind, ind)]
                if verbose: print(f"initial elastic tensor:\t{tab.join(standard_axes)}\n"
                                    + f"{nl.join([tab + tab + standard_axes[i] + tab + tab.join(row.round(2).astype(str)) for i, row in enumerate(voigtA)])}")
            elif info.startswith(('F', 'f')):
                if voigtB: print("Warning: Final elastic tensor is defined multiple times! The last definition will be used.")
                voigtB = np.zeros((6, 6))
                axes = []
                for i in range(6):
                    l = f.readline().strip().split()
                    axes.append(l[0])
                    voigtB[i,:] = [float(x) for x in l[1:7]]
                axis_yz = 'YZ' if 'YZ' in axes else 'ZY'
                axis_zx = 'ZX' if 'ZX' in axes else 'XZ'
                axis_xy = 'XY' if 'XY' in axes else 'YX'
                ind = [axes.index('XX'), axes.index('YY'), axes.index('ZZ'), axes.index(axis_yz), axes.index(axis_zx), axes.index(axis_xy)]
                voigtB = voigtB[np.ix_(ind, ind)]
                if verbose: print(f"final elastic tensor:\t{tab.join(standard_axes)}\n"
                                    + f"{nl.join([tab + tab + standard_axes[i] + tab + tab.join(row.round(2).astype(str)) for i, row in enumerate(voigtB)])}")
            elif info.startswith(('A', 'a', 'W', 'w', 'S', 's', 'D', 'd')):
                species = f.readline().strip().split()
                weights = np.array(f.readline().strip().split(), dtype=float)
                weight_func = {s: w for s, w in zip(species, weights)}
                if verbose: print(f"atomic weights used to define the shuffle distance:\n{nl.join([tab + key + tab + str(value) for key, value in weight_func.items()])}")
            elif info.startswith(('O', 'o')):
                para1 = f.readline().split("||")
                h_i1, k_i1, l_i1 = [float(x) for x in para1[0].strip().split()]
                h_f1, k_f1, l_f1 = [float(x) for x in para1[1].strip().split()]
                para2 = f.readline().split("||")
                h_i2, k_i2, l_i2 = [float(x) for x in para2[0].strip().split()]
                h_f2, k_f2, l_f2 = [float(x) for x in para2[1].strip().split()]
                ori_rel = (((h_i1, k_i1, l_i1), (h_f1, k_f1, l_f1)), ((h_i2, k_i2, l_i2), (h_f2, k_f2, l_f2)))
                if verbose:
                    print("orientation relationship:")
                    print(f"\t({h_i1:.3f}, {k_i1:.3f}, {l_i1:.3f}) || ({h_f1:.3f}, {k_f1:.3f}, {l_f1:.3f})")
                    print(f"\t({h_i2:.3f}, {k_i2:.3f}, {l_i2:.3f}) || ({h_f2:.3f}, {k_f2:.3f}, {l_f2:.3f})")
            info = f.readline()
    return voigtA, voigtB, weight_func, ori_rel

def unique_filename(
    message: Optional[str],
    filename: str,
    ext: bool = True
) -> str:
    """Generate a unique filename by appending a counter to the base name if the file already exists.
    
    Parameters
    ----------
    message : str, optional
        A message to print when generating the unique filename.
    filename : str
        The original filename.
    ext : bool, optional
        If True, assume that `filename` has an extension. This should be set to False if `filename` does not have an extension.

    Returns
    -------
    new_filename : str
        The unique filename, which is the original filename with a counter appended if necessary.
    """
    if ext: base, ext = splitext(filename)
    counter = 1
    new_filename = filename
    while exists(new_filename):
        new_filename = f"{base}-{counter}{ext}" if ext else f"{filename}-{counter}"
        counter += 1
    if message != None: print(f"{message} '{new_filename}' ...")
    return new_filename

def cryst_to_spglib(
    cryst: Cryst,
    return_dict: bool = False
) -> Union[tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int32]],
tuple[tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int32]], Dict[str, int]]]:
    """ Convert a crystal structure to a spglib-compatible format.

    Parameters
    ----------
    cryst : cryst
        The crystal structure, represented as a tuple of `(lattice, species, positions)`.
    return_dict : bool, optional
        If True, return a dictionary mapping species to unique numbers.

    Returns
    -------
    spglib_cell : 3-tuple
        The spglib-compatible format.
    species_dict : dict, optional
        A dictionary mapping species to unique numbers. Only returned if `return_dict` is True.
    """
    species_dict, numbers = np.unique(cryst[1], return_inverse=True)
    spglib_cell = (cryst[0], cryst[2], numbers)
    if return_dict: return spglib_cell, species_dict
    else: return spglib_cell

def spglib_to_cryst(
    spglib_cell: tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int32]],
    species_dict: Dict[str, int]
) -> Cryst:
    """Convert a spglib-compatible cell to a crystal structure.
    
    Parameters
    ----------
    spglib_cell : 3-tuple
        The spglib-compatible cell, represented as `(lattice, positions, species_numbers)`.
    species_dict : dict
        A dictionary mapping species to unique numbers, usually obtained from `cryst_to_spglib()`.
    
    Returns
    -------
    cryst : cryst
        The crystal structure, represented as a tuple of `(lattice, species, positions)`.
    """
    return (spglib_cell[0], species_dict[spglib_cell[2]], spglib_cell[1])

def primitive_cryst(
    cryst_sup: Cryst,
    tol: float = 1e-3
) -> Cryst:
    """Reduce a crystal structure to its primitive structure using spglib.
    
    Parameters
    ----------
    cryst_sup : cryst
        The crystal structure to be reduced, represented as a tuple of `(lattice, species, positions)`.
    tol : float, optional
        The tolerance for symmetry detection.
        
    Returns
    -------
    cryst : cryst
        The primitive crystal structure, represented as a tuple of `(lattice, species, positions)`.
    """
    cell_sup, species_dict = cryst_to_spglib(cryst_sup, return_dict=True)
    cell = standardize_cell(cell_sup, to_primitive=True, no_idealize=True, symprec=tol)
    return spglib_to_cryst(cell, species_dict)

def check_stoichiometry(
    speciesA: NDArray[np.int32],
    speciesB: NDArray[np.int32]
) -> None:
    """Check if the stoichiometry of two species lists is the same.
    
    Parameters
    ----------
    speciesA, speciesB : array-like
        The species lists to be checked, which should be 1D arrays of atomic species.
    """
    spA, ctA = np.unique(speciesA, return_counts=True)
    spB, ctB = np.unique(speciesB, return_counts=True)
    if not (spA == spB).all(): raise ValueError("Atomic species are not the same!")
    if not (ctB.sum() * ctA == ctA.sum() * ctB).all(): raise ValueError("Stoichiometric ratios are not the same!")
    return

def create_common_supercell(
    crystA: Cryst,
    crystB: Cryst,
    slm: SLM
) -> tuple[Cryst, Cryst, NDArray[np.float64], NDArray[np.int32], NDArray[np.int32]]:
    """Return the initial, final, and half-distorted supercell, as well as the transformation matrices.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slm : slm
        The SLM of the CSM.
    
    Returns
    -------
    crystA_sup : cryst
        The supercell structure of `crystA`.
    crystB_sup : cryst
        The supercell structure of `crystB`.
    c_sup_half : (3, 3) array of floats
        The half-distorted supercell.
    mA : (3, 3) array of ints
        The matrix that transforms `crystA` to `crystA_sup`.
    mB : (3, 3) array of ints
        The matrix that transforms `crystB` to `crystB_sup`.
    """
    # Unpacking crystal structures.
    cA = crystA[0].T
    cB = crystB[0].T
    speciesA = crystA[1]
    speciesB = crystB[1]
    pA = crystA[2].T
    pB = crystB[2].T
    check_stoichiometry(speciesA, speciesB)
    
    # Determining the supercell geometries from the SLM.
    hA, hB, q = slm
    u, sigma, vT = la.svd(deformation_gradient(crystA, crystB, slm))
    c_sup_half, q_sup = lll_reduce(vT.T @ np.diag(sigma ** 0.5) @ vT @ cA @ hA)         # The half-distorted supercell.
    mA = hA @ q_sup
    mB = hB @ q @ q_sup
    cA_sup = cA @ mA
    cB_sup = cB @ mB
    
    # Sorting supercell species and positions.
    speciesA_sup = np.tile(speciesA, la.det(mA).round().astype(int))
    speciesB_sup = np.tile(speciesB, la.det(mB).round().astype(int))
    pA_sup = (la.inv(mA) @ (pA.reshape(3,1,-1) + int_vec_inside(mA).reshape(3,-1,1)).reshape(3,-1))
    pB_sup = (la.inv(mB) @ (pB.reshape(3,1,-1) + int_vec_inside(mB).reshape(3,-1,1)).reshape(3,-1))
    argsortA = np.argsort(speciesA_sup)
    argsortB = np.argsort(speciesB_sup)
    if not (speciesA_sup[argsortA] == speciesB_sup[argsortB]).all():
        raise AssertionError("Species array is ill-sorted. Please report this bug to wfc@pku.edu.cn.")
    species_sup = speciesA_sup[argsortA]
    pA_sup = pA_sup[:,argsortA]
    pB_sup = pB_sup[:,argsortB]
    
    # Computing output.
    crystA_sup = (cA_sup.T, species_sup, pA_sup.T)
    crystB_sup = (cB_sup.T, species_sup, pB_sup.T)
    return crystA_sup, crystB_sup, c_sup_half, mA, mB

def frac_cell(
    mA: NDArray[np.int32],
    mB: NDArray[np.int32]
) -> NDArray[np.float64]:
    """The primitive cell of the lattice generated by mA^{-1} and mB^{-1}.
    
    Parameters
    ----------
    mA : (3, 3) array of ints
        The transformation matrix from the initial crystal structure to the supercell.
    mB : (3, 3) array of ints
        The transformation matrix from the final crystal structure to the supercell.
        
    Returns
    -------
    fcell : (3, 3) array of floats
        The primitive cell of the lattice generated by mA^{-1} and mB^{-1}.
    """
    lcm = np.lcm(la.det(mA).round().astype(int), la.det(mB).round().astype(int))
    dual = (np.hstack((la.inv(mA), la.inv(mB))) * lcm).round().astype(int)
    return lll_reduce(_hnf(dual)[:,:3] / lcm)[0]

def imt_multiplicity(
    crystA: Cryst,
    crystB: Cryst,
    slmlist: NDArray[np.int32]
) -> int | NDArray[np.int32]:
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
    if slmlist.ndim == 3: return la.det(slmlist[0]).round().astype(int) // dA
    else: return la.det(slmlist[:,0,:,:]).round().astype(int) // dA

def deformation_gradient(
    crystA: Cryst,
    crystB: Cryst,
    slmlist: NDArray[np.int32]
) -> NDArray[np.float64]:
    """Compute the deformation gradient matrices of given IMTs.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slmlist : list of slm
        A list of uuLMs, each represented by a triplet of integer matrices like `(hA, hB, q)`.
    
    Returns
    -------
    slist : (..., 3, 3) array
        A list of deformation gradient matrices.
    """
    cA = crystA[0].T
    cB = crystB[0].T
    slms = np.array(slmlist)
    hA = slms[...,0,:,:]
    hB = slms[...,1,:,:]
    q = slms[...,2,:,:]
    return cB @ hB @ q @ la.inv(cA @ hA)

def rmss(
    slist: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Root-mean-square strains of given deformation gradient matrices.

    Parameters
    ----------
    slist : (..., 3, 3) array
        A list of deformation gradient matrices.
    
    Returns
    -------
    rmss : (...) array
        Root-mean-square strains.
    """
    return np.sqrt(np.mean((la.svd(slist, compute_uv=False) - 1) ** 2, axis=-1))

def zip_pct(
    p: NDArray[np.int32],
    ks: NDArray[np.int32]
) -> NDArray[np.int32]:
    """Zip the permutation `p` and the class-wise translations `ks` into a single PCT array.
    
    Parameters
    ----------
    p : (N,) array of ints
        The permutation part of the PCT.
    ks : (3, N) array of ints
        The class-wise translations of the PCT.
    
    Returns
    -------
    pct : (N, 4) array of ints
        The zipped PCT array.
    """
    return np.hstack((p.reshape(-1,1), ks.T), dtype=int)

def unzip_pct(
    pct: NDArray[np.int32]
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Unzip the PCT array into the permutation and class-wise translations.
    
    Parameters
    ----------
    pct : (N, 4) array of ints
        The PCT array, where the first column is the permutation and the rest are class-wise translations.
    
    Returns
    -------
    p : (N,) array of ints
        The permutation part of the PCT.
    ks : (3, N) array of ints
        The class-wise translations of the PCT.
    """
    return pct[:,0], pct[:,1:].T

def get_pure_rotation(
    cryst: Cryst,
    tol: float = 1e-3
) -> NDArray[np.int32]:
    """Find all pure rotations appeared in the space group of `cryst`.

    Parameters
    ----------
    cryst : cryst
        The crystal structure, usually obtained by `load_poscar`.
    tol : float, optional
        The tolerance for `spglib` symmetry detection.
    
    Returns
    -------
    g : (..., 3, 3) array of ints
        A point group of the first kind, containing all pure rotations appeared in the space group of `cryst`, \
            elements of which are integer matrices (using the columns of `cryst[1]` as basis vectors).
    """
    g = get_symmetry(cryst_to_spglib(cryst), symprec=tol)['rotations']
    g = g[la.det(g).round(decimals=4)==1,:,:]
    g = np.unique(g, axis=0)
    return g

@nb.njit
def _mul_xor_hash(arr, init=65537, k=37):
    # This function is provided by @norok2 on StackOverflow: https://stackoverflow.com/a/66674679.
    result = init
    for x in arr.view(np.uint64):
        result = (result * k) ^ x
    return result

@nb.njit
def setdiff2d(
    arr1: NDArray[np.int32],
    arr2: NDArray[np.int32]
) -> NDArray[np.int32]:
    """Find the set difference of two 2D integer arrays, `arr1` and `arr2`.
    This function is provided by @norok2 on StackOverflow: https://stackoverflow.com/a/66674679.
    
    Parameters
    ----------
    arr1 : (N, M) array of ints
        The first array, from which elements will be removed.
    arr2 : (P, M) array of ints
        The second array, containing elements to be removed from `arr1`.
        
    Returns
    -------
    result : (Q, M) array of ints
        The set difference of `arr1` and `arr2`.
    """
    delta = {_mul_xor_hash(arr2[0])}
    for i in range(1, arr2.shape[0]):
        delta.add(_mul_xor_hash(arr2[i]))
    n = 0
    for i in range(arr1.shape[0]):
        if _mul_xor_hash(arr1[i]) not in delta:
            n += 1
    result = np.empty((n, arr1.shape[-1]), dtype=arr1.dtype)
    j = 0
    for i in range(arr1.shape[0]):
        if _mul_xor_hash(arr1[i]) not in delta:
            result[j] = arr1[i]
            j += 1
    return result

def int_vec_inside(c: NDArray[np.int32]) -> NDArray[np.int64]:
    """Integer vectors inside the cell `c @ [0, 1)^3` whose elements are integers.

    Parameters
    ----------
    c : (3, 3) array of ints
        A matrix whose columns are integer cell vectors.
    
    Returns
    -------
    v_int : (3, ...) array of ints
        Its columns are vectors satisfying `v = c @ k`, where `k[0]`, `k[1]`, `k[2]` are all in `[0, 1)`.
    """
    if not c.dtype == int: raise TypeError(f"Input matrix must be integer.")
    vertices = c @ np.mgrid[0:2,0:2,0:2].reshape(3,-1)
    candidates = np.mgrid[np.amin(vertices[0,:]):np.amax(vertices[0,:])+1, np.amin(vertices[1,:]):np.amax(vertices[1,:])+1, \
        np.amin(vertices[2,:]):np.amax(vertices[2,:])+1].reshape(3,-1)
    fractional = (la.inv(c) @ candidates).round(decimals=7)
    is_inside = (np.prod(fractional < 1, axis=0) * np.prod(fractional >= 0, axis=0)).astype(bool)
    assert np.sum(is_inside) == np.abs(la.det(c)).round().astype(int)
    return candidates[:,is_inside]

def _int_fact(n: int) -> list[tuple[int, int]]:
    l = []
    for a in range(1,n+1):
        if n % a == 0: l.append((a, n//a))
    return l

def all_hnfs(det: int) -> NDArray[np.int32]:
    """Enumerate all (3, 3) Hermite normal forms (HNFs) with given determinant.

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
    for a, aa in _int_fact(det):
        for b, c in _int_fact(aa):
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

def _hnf(m: NDArray[np.int32]) -> NDArray[np.int32]:
    h = m.copy()
    n_row = h.shape[0]
    for i in range(n_row):
        while (h[i,i+1:] != 0).any():
            col_nonzero = i + np.nonzero(h[i,i:])[0]
            i0 = col_nonzero[np.argpartition(np.abs(h[i,col_nonzero]), kth=0)[0]]
            h[:,[i,i0]] = h[:,[i0,i]]
            if h[i,i] < 0: h[:,i] = - h[:,i]
            h[:,i+1:] = h[:,i+1:] - np.outer(h[:,i], h[i,i+1:] // h[i,i])
        if h[i,i] < 0: h[:,i] = - h[:,i]
        h[:,:i] = h[:,:i] - np.outer(h[:,i], h[i,:i] // h[i,i])
    return h

def hnf(
    m: NDArray[np.int32],
    return_q: bool = False
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Compute the Hermite normal form of integer matrix `m`.

    Parameters
    ----------
    m : (M, N) array of ints
        The integer matrix to decompose, with positive determinant and M <= N.
    return_q : bool, optional
        Whether to return the unimodular matrix `q`.
    
    Returns
    -------
    h : (M, N) array of ints
        The column-style Hermite normal form of `m`.
    q : (N, N) array of ints
        The unimodular matrix satisfying `m` = `h @ q`. Only returned if `return_q` is True.
    """
    if not m.dtype == int: raise TypeError(f"Input matrix must be integer:\n{m}")
    if not la.matrix_rank(m, tol=1e-6) == m.shape[0]: raise ValueError(f"Input matrix must be full-row-rank:\n{m}")
    h = _hnf(m)
    if return_q:
        if not m.shape[1] == m.shape[0]: raise ValueError(f"Input matrix must be square:\n{m}")
        return h, (la.inv(h) @ m).round().astype(int)
    else:
        return h

def standardize_imt(
    slm: NDArray[np.int32],
    gA: NDArray[np.int32],
    gB: NDArray[np.int32]
) -> SLM:
    """The standard SLM of the congruence class of `slm`.

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
    slm0 : slm
        The standardized SLM.
    """
    hA, hB, q = slm
    cl = np.dot((gB @ hB) @ q, la.inv(gA @ hA)).transpose((2,0,1,3)).reshape(-1,9).round(4)
    ind = np.arange(cl.shape[0])
    for j in range(9):
        ind = ind[cl[ind,j] == np.min(cl[ind,j])]
        if len(ind) == 1: break
    iA, iB = np.unravel_index(ind, (gA.shape[0], gB.shape[0]))
    i = np.lexsort(np.array([gA[iA], gB[iB]]).transpose((0,2,3,1)).reshape(18,-1))[0]
    hAA, qA = hnf(gA[iA[i]] @ hA, return_q=True)
    hBB, qB = hnf(gB[iB[i]] @ hB, return_q=True)
    slm0 = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return slm0

def _gram_schmidt(c: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = c.shape[1]
    cc = np.zeros_like(c, dtype=float)
    mu = np.zeros((n, n), dtype=float)
    norm2 = np.zeros(n, dtype=float)

    for i in range(n):
        v = c[:,i].copy()
        for j in range(i):
            mu[i,j] = np.dot(c[:,i], cc[:,j]) / norm2[j]
            v -= mu[i,j] * cc[:,j]
        cc[:,i] = v
        norm2[i] = np.dot(v, v)
    return cc, mu, norm2

def lll_reduce(
    c: NDArray[np.float64],
    delta: float = 0.75
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """Return the LLL-reduced cell `cc` and the unimodular matrix `q` such that `cc = c @ q`

    Parameters
    ----------
    c : (3, 3) array
        The cell to be reduced, whose columns are cell vectors.
    delta : float, optional
        The parameter for the LLL algorithm.
    
    Returns
    -------
    cc : (3, 3) array
        The LLL-reduced cell.
    q : (3, 3) array of ints
        The unimodular matrix satisfying `cc = c @ q`.
    """
    cc = c.copy()
    n = cc.shape[1]
    k = 1
    _, mu, norm2 = _gram_schmidt(cc)

    while k < n:
        for j in range(k-1, -1, -1):
            q = round(mu[k,j])
            if abs(q) > 0:
                cc[:,k] -= q * cc[:,j]
        _, mu, norm2 = _gram_schmidt(cc)

        if norm2[k] >= (delta - mu[k,k-1] ** 2) * norm2[k-1]:
            k += 1
        else:
            cc[:,[k,k-1]] = cc[:,[k-1,k]]
            _, mu, norm2 = _gram_schmidt(cc)
            k = max(k-1, 1)
    
    if la.det(cc) < 0: cc[:,-1] = -cc[:,-1]
    q = (la.inv(c) @ cc).round().astype(int)
    return cc, q

def voigt_to_tensor(
    voigt_matrix: NDArray[np.float64],
    cryst: Cryst = None,
    tol: float = 1e-3,
    verbose: bool = True
) -> NDArray[np.float64]:
    """Convert a Voigt-notation tensor to a rank-4 tensor, and symmetrize it according to the symmetry of `cryst` (if provided).
    
    Parameters
    ----------
    voigt_matrix : (6, 6) array
        The elastic tensor, in Voigt notation (ordered as XX, YY, ZZ, YZ, ZX, XY).
    cryst : Cryst, optional
        The crystal structure, whose symmetry is used to symmetrize the elastic tensor.
    tol : float, optional
        The tolerance for symmetry finding.
    verbose : bool, optional
        Whether to print information about the symmetrization.

    Returns
    -------
    tensor : (3, 3, 3, 3) array
        The rank-4 elastic tensor.
    """
    tensor = np.ones((9,9)) * np.inf                # XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
    voigt_ind = [0,4,8,5,6,1]                       # XX,YY,ZZ,YZ,ZX,XY
    tensor[np.ix_(voigt_ind, voigt_ind)] = voigt_matrix
    tensor = tensor.reshape((3,3,3,3))
    tensor = np.min([tensor, tensor.transpose((1,0,2,3))], axis=0)
    tensor = np.min([tensor, tensor.transpose((2,3,0,1))], axis=0)
    tensor = np.min([tensor, tensor.transpose((0,1,3,2))], axis=0)
    if cryst is not None:
        spglib_cryst = (cryst[0],cryst[2],np.unique(cryst[1], return_inverse=True)[1])
        g = get_symmetry(spglib_cryst, symprec=tol)['rotations']
        r = cryst[0].T @ g @ la.inv(cryst[0].T)
        tensor_sym = np.mean(np.einsum('ijkl,qim,qjn,qko,qlp->qmnop', tensor, r, r, r, r), axis=0)
        spg = get_spacegroup(spglib_cryst, symprec=tol)
        dev = np.max(np.abs(tensor_sym - tensor)) / np.max(np.abs(tensor))
        if verbose:
            print(f"Symmetrizing elastic tensor using space group {spg} ...")
            print(f"\tmax |Y_ijkl| = {np.max(np.abs(tensor_sym)):.3f}")
            print(f"\tmax |ΔY_ijkl| = {np.max(np.abs(tensor_sym - tensor)):.3f}")
            print(f"\tdeviation = {100 * dev:.2f}%")
        if dev > 0.2: print(f"\nWarning: Elastic tensor does not have the expected symmetry ({spg})! Check if the input POSCAR and elastic tensor are consistent.\n")
        return tensor_sym
    else: return tensor

def pct_distance(
    c: NDArray[np.float64],
    pA: NDArray[np.float64],
    pB: NDArray[np.float64],
    p: NDArray[np.int32],
    ks: NDArray[np.int32],
    weights: Optional[NDArray[np.float64]] = None,
    l: float = 2.0,
    min_t0: bool = True,
    return_t0: bool = False
) -> float | tuple[float, NDArray[np.float64]]:
    """Return the shuffle distance of a PCT.
    
    Parameters
    ----------
    c : (3, 3) array
        The lattice vectors of the crystal structure.
    pA, pB : (3, Z) array
        The fractional coordinates of the atoms in the initial and final structures, respectively.
    p : (Z, ) array of ints
        The permutation of the atoms.
    ks : (3, Z) array of ints
        The class-wise translations (fractional coordinates) of the atoms in `pB`.
    weights : (Z, ) array of floats, optional
        The weights of each atom. If None, all atoms have the same weight.
    l : float, optional
        The l-norm to be used for distance calculation, must not be less than 1.
    min_t0 : bool, optional
        Set to True to minimize the shuffle distance by translating the final structure.
    return_t0 : bool, optional
        Whether to return the best overall translation if `min_t0` is True.
    
    Returns
    -------
    distance : float
        The shuffle distance.
    t0 : (3, 1) array
        The best overall translation, reshaped as a 3x1 matrix. Only returned if `return_t0` is True.
    """
    if ks.shape != (3,len(p)): raise ValueError("'p' and 'ks' must have the same number of atoms.")
    if not min_t0: return np.average(((c @ (pB[:,p] + ks - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l)
    if np.allclose(l, 2.0, atol=1e-6):
        t0 = -np.average(pB[:,p] + ks - pA, axis=1, weights=weights, keepdims=True)
        d = np.average(((c @ (pB[:,p] + ks + t0 - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l)
    else:
        res = minimize(lambda t: np.average(la.norm(c @ (pB[:,p] + ks + t.reshape(3,1) - pA), axis=0)**l, weights=weights),
                        -np.average(pB[:,p] + ks - pA, axis=1, weights=weights), method='SLSQP', options={'disp': False})
        d = res.fun ** (1/l)
        t0 = res.x.reshape(3,1)
    if return_t0: return d, t0
    else: return d