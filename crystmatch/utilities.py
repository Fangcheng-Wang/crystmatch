"""
Load/save crystal structures and CSMs from/to files.
"""

from os import makedirs, environ
from os.path import sep, exists, splitext
import numpy as np
import numba as nb
import numpy.linalg as la
from collections import namedtuple
from spglib import get_spacegroup, get_symmetry, standardize_cell, refine_cell, get_symmetry_dataset
from tqdm import tqdm
from typing import Union, Tuple, List, Dict, Callable
from numpy.typing import NDArray, ArrayLike
from scipy.spatial.transform import Rotation

np.set_printoptions(suppress=True)
Cryst = Tuple[NDArray[np.float64], NDArray[np.str_], NDArray[np.float64]]
SLM = Tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int32]]

def load_poscar(filename: str, to_primitive: bool = True, tol: float = 1e-3, verbose: bool = True) -> Cryst:
    """Load the crystal structure from a POSCAR file.

    Parameters
    ----------
    filename : str
        The name of the POSCAR file to be read.
    to_primitive : bool, optional
        Using the primitive cell instead of the cell given by the POSCAR file. Default is True.
    tol : float, optional
        The tolerance for `spglib` symmetry detection; default is 1e-3.

    Returns
    -------
    cryst : cryst
        The loaded crystal structure, consisting of the lattice vectors, species, and positions.
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
    cryst = (lattice, species, positions)

    if verbose: print(f"\tSpace group: {get_spacegroup(cryst_to_spglib(cryst), symprec=tol)}.")
    if to_primitive:
        cryst = primitive_cryst(cryst, tol=tol)
        if verbose:
            if len(cryst[1]) != len(species): print(f"\tCell in POSCAR file is not primitive! Using primitive cell (Z = {len(cryst[1]):d}) now.")
            else: print(f"\tCell in POSCAR file is already primitive (Z = {len(cryst[1]):d}).")
    elif verbose: print(f"\tUsing cell in POSCAR file (Z = {len(species):d}).")
    return cryst

def load_csmcar(filename: str, verbose: bool = True):
    """Load the CSMCAR file, which contains `crystmatch` parameters.

    Parameters
    ----------
    filename : str
        The name of the POSCAR file to be read.

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
                if verbose: print(f"Initial elastic tensor:\t{tab.join(standard_axes)}\n{nl.join([tab + standard_axes[i] + tab + tab.join(row.astype(str)) for i, row in enumerate(voigtA)])}")
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
                if verbose: print(f"Final elastic tensor:\t{tab.join(standard_axes)}\n{nl.join([tab + standard_axes[i] + tab + tab.join(row.astype(str)) for i, row in enumerate(voigtB)])}")
            elif info.startswith(('D', 'd')):
                species = f.readline().strip().split()
                weights = np.array(f.readline().strip().split(), dtype=float)
                weight_func = {s: w for s, w in zip(species, weights)}
                if verbose: print(f"Distance weight function:\n{nl.join([tab + key + tab + str(value) for key, value in weight_func.items()])}")
            elif info.startswith(('O', 'o')):
                vs = f.readline().split("||")
                vix, viy, viz = [float(x) for x in vs[0].strip().split()]
                vfx, vfy, vfz = [float(x) for x in vs[1].strip().split()]
                ws = f.readline().split("||")
                wix, wiy, wiz = [float(x) for x in ws[0].strip().split()]
                wfx, wfy, wfz = [float(x) for x in ws[1].strip().split()]
                ori_rel = np.array([[[vix, viy, viz], [vfx, vfy, vfz]], [[wix, wiy, wiz], [wfx, wfy, wfz]]])
                if verbose:
                    print(f"Orientation relationship:\n\t({vix:.3f}, {viy:.3f}, {viz:.3f}) || ({vfx:.3f}, {vfy:.3f}, {vfz:.3f})")
                    print(f"\t({wix:.3f}, {wiy:.3f}, {wiz:.3f}) || ({wfx:.3f}, {wfy:.3f}, {wfz:.3f})")
            info = f.readline()
    return voigtA, voigtB, weight_func, ori_rel

def unique_filename(message: Union[str, None], filename: str) -> str:
    """Get a unique filename by appending a number to the end of the given filename.

    Parameters
    ----------
    filename : str
        The filename to be modified.
    message : str, optional
        A message to print before the filename.

    Returns
    -------
    new_filename : str
        The modified filename with a unique number appended.
    """
    base, ext = splitext(filename)
    counter = 1
    new_filename = filename
    while exists(new_filename):
        new_filename = f"{base}-{counter}{ext}"
        counter += 1
    if message != None: print(f"{message} '{new_filename}' ...")
    return new_filename

def species_poscar_format(species: NDArray[np.str_]) -> Tuple[NDArray[np.str_], NDArray[np.int32]]:
    """
    Examine whether a species array is sorted. If so, convert it to the POSCAR format.
    
    Parameters
    ----------
    species : (N,) array of str
        The species array.
    
    Returns
    -------
    species_unique : (M,) array of str
        The unique species in the order of their first appearance in `species`.
    species_counts : (M,) array of int
        The number of occurrences of each unique species in `species`.
    """
    _, sp_idx, sp_inv, sp_counts = np.unique(species, return_index=True, return_inverse=True, return_counts=True)
    if np.sum(np.diff(sp_inv) != 0) != sp_idx.shape[0] - 1:
        raise ValueError("Species array is ill-sorted. Please report this bug to wfc@pku.edu.cn if you see this message.")
    return species[np.sort(sp_idx)], sp_counts

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
    species_name, species_counts = species_poscar_format(cryst[1])
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

def cryst_to_spglib(cryst, return_dict=False):
    species_dict, numbers = np.unique(cryst[1], return_inverse=True)
    spglib_cell = (cryst[0], cryst[2], numbers)
    if return_dict: return spglib_cell, species_dict
    else: return spglib_cell

def spglib_to_cryst(spglib_cell, species_dict):
    return (spglib_cell[0], species_dict[spglib_cell[2]], spglib_cell[1])

def primitive_cryst(cryst_sup, tol=1e-3):
    cell_sup, species_dict = cryst_to_spglib(cryst_sup, return_dict=True)
    cell = standardize_cell(cell_sup, to_primitive=True, no_idealize=True, symprec=tol)
    return spglib_to_cryst(cell, species_dict)

def supercell_decomposition(cryst_sup, cryst=None, tol=1e-3):
    """Compute the rotation `r` and the integer transformation matrix `m` such that `c_sup = r @ c @ m`
    """
    if cryst is None:
        q = np.eye(3, dtype=int)
        r0 = np.eye(3)
    else:
        sym0 = get_symmetry_dataset(cryst_to_spglib(cryst), symprec=tol)
        q = sym0.transformation_matrix.round().astype(int)
        r0 = sym0.std_rotation_matrix
        if not np.abs(la.det(q)).round(8) == 1: raise ValueError("'cryst' must be primitive or left as None.")
    sym = get_symmetry_dataset(cryst_to_spglib(cryst_sup), symprec=tol)
    m = (la.inv(q) @ sym.transformation_matrix).round().astype(int)
    r = r0 @ sym.std_rotation_matrix.T
    return r, m

def triangularize_cryst(cryst_sup, return_primitive=False, tol=1e-3):
    """Rotate the crystal structure and change its lattice basis such that `c` is lower triangular.
    """
    c_sup = cryst_sup[0].T
    p_sup = cryst_sup[2].T
    r, m = supercell_decomposition(cryst_sup, tol=tol)
    h, q = hnf_int(m, return_q=True)
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

def check_chem_comp(speciesA, speciesB):
    spA, ctA = np.unique(speciesA, return_counts=True)
    spB, ctB = np.unique(speciesB, return_counts=True)
    if not (spA == spB).all(): raise ValueError("Atomic species are not the same!")
    if not (ctB.sum() * ctA == ctA.sum() * ctB).all(): raise ValueError("Stoichiometric ratios are not the same!")
    return

def create_common_supercell(crystA: Cryst, crystB: Cryst, slm: SLM) -> Tuple[Cryst, Cryst, NDArray[np.float64], NDArray[np.float64]]:
    """Create supercell structures representing $\mathcal{A}$ (initial structure)and $\sqrt{S^{\\text{T}} S}S^{-1}\mathcal{B}$ \
        (rotation-free final structure). Also return the half-distorted supercell and translation cell.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slm : slm
        The SLM of the CSM.
    
    Returns
    -------
    crystA_sup : cryst
        The supercell of $\mathcal{A}$.
    crystB_sup_final : cryst
        The supercell of $\sqrt{S^{\\text{T}} S}\mathcal{A}$, which equals the supercell of $(S^{\\text{T}} S)^{-1/2}\mathcal{B}$.
    c_sup_half : (3, 3) array of floats
        The half-distorted supercell, which equals the supercel of $(S^{\\text{T}} S)^{1/4}\mathcal{A}$ and that of $(S^{\\text{T}} S)^{-1/4}\mathcal{B}$.
    mA, mB : (3, 3) array of ints
        The matrices that transform `crystA` to `crystA_sup` and `crystB` to `crystB_sup`, respectively.
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
    c_sup_half, q_sup = cell_reduce(vT.T @ np.diag(sigma ** 0.5) @ vT @ cA @ hA)         # The half-distorted supercell.
    mA = hA @ q_sup
    mB = hB @ q @ q_sup
    cA_sup = cA @ mA
    cB_sup_final = (u @ vT).T @ cB @ mB                     # The rotation-free orientation of `crystB_sup`.
    
    # Sorting supercell species and positions.
    speciesA_sup = np.tile(speciesA, la.det(mA).round().astype(int))
    speciesB_sup = np.tile(speciesB, la.det(mB).round().astype(int))
    pA_sup = (la.inv(mA) @ (pA.reshape(3,1,-1) + int_vec_inside(mA).reshape(3,-1,1)).reshape(3,-1))
    pB_sup = (la.inv(mB) @ (pB.reshape(3,1,-1) + int_vec_inside(mB).reshape(3,-1,1)).reshape(3,-1))
    argsortA = np.argsort(speciesA_sup)
    argsortB = np.argsort(speciesB_sup)
    if not (speciesA_sup[argsortA] == speciesB_sup[argsortB]).all():
        raise ValueError("Species array is ill-sorted. Please report this bug to wfc@pku.edu.cn if you see this message.")
    species_sup = speciesA_sup[argsortA]
    pA_sup = pA_sup[:,argsortA]
    pB_sup = pB_sup[:,argsortB]
    
    # Computing output.
    crystA_sup = (cA_sup.T, species_sup, pA_sup.T)
    crystB_sup_final = (cB_sup_final.T, species_sup, pB_sup.T)
    return crystA_sup, crystB_sup_final, c_sup_half, mA, mB

def f_translate(mA, mB):
    """The translation cell (fractional coordinates) generated by mA^{-1} and mB^{-1}.
    """
    return cell_reduce(matrix_gcd(la.inv(mA), la.inv(mB)))[0]

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
    crystA_sup, crystB_sup, _, _, _ = create_common_supercell(crystA, crystB, slm)
    if not (crystA_sup[1] == crystB_sup[1][p]).all():
        raise ValueError("Species array is ill-sorted. Please report this bug to wfc@pku.edu.cn if you see this message.")
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T[:,p] + ks
    if centered:
        pB_sup = pB_sup - np.mean(pB_sup - pA_sup, axis=1, keepdims=True)
    return crystA_sup, (crystB_sup[0], crystB_sup[1][p], pB_sup.T)

def imt_multiplicity(crystA: Cryst, crystB: Cryst, slmlist: Union[SLM, List[SLM], NDArray[np.int32]]) -> Union[int, NDArray[np.int32]]:
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

def cube_to_so3(vec):
    """Map `vec` in [0,1)^3 to SO(3), preserving the uniformity of distribution.
    """
    q0 = np.sqrt(1 - vec[0]) * np.sin(2 * np.pi * vec[1])
    q1 = np.sqrt(1 - vec[0]) * np.cos(2 * np.pi * vec[1])
    q2 = np.sqrt(vec[0]) * np.sin(2 * np.pi * vec[2])
    q3 = np.sqrt(vec[0]) * np.cos(2 * np.pi * vec[2])
    return Rotation.from_quat([q0,q1,q2,q3]).as_matrix()

def deformation_gradient(crystA: Cryst, crystB: Cryst, slmlist: List[SLM]) -> NDArray[np.float64]:
    """Compute the deformation gradient matrices of given IMTs.
    
    Parameters
    ----------
    crystA, crystB : cryst
        The initial and final structures.
    slmlist : list of slm
        The IMTs.
    
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

def rmss(slist: NDArray[np.float64]) -> NDArray[np.float64]:
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

def zip_pct(p, ks):
    return np.hstack((p.reshape(-1,1), ks.T), dtype=int)

def unzip_pct(pct):
    return pct[:,0], pct[:,1:].T

def get_pure_rotation(cryst: Cryst, tol: float = 1e-3) -> NDArray[np.int32]:
    """Find all pure rotations appeared in the space group of `cryst`.

    Parameters
    ----------
    cryst : 3-tuple
        `(lattice, species, positions)`, representing the crystal structure, usually obtained by `load_poscar`.
    tol : float, optional
        The tolerance for `spglib` symmetry detection; default is 1e-3.
    
    Returns
    -------
    g : (..., 3, 3) array of ints
        A point group of the first kind, containing all pure rotations appeared in the space group of `cryst`, \
            elements of which are integer matrices (under fractional coordinates).
    """
    g = get_symmetry(cryst_to_spglib(cryst), symprec=tol)['rotations']
    g = g[la.det(g).round(decimals=4)==1,:,:]
    g = np.unique(g, axis=0)
    return g


@nb.njit
def mul_xor_hash(arr, init=65537, k=37):
    """This function is provided by @norok2 on StackOverflow: https://stackoverflow.com/a/66674679.
    """
    result = init
    for x in arr.view(np.uint64):
        result = (result * k) ^ x
    return result

@nb.njit
def setdiff2d(arr1, arr2):
    """This function is provided by @norok2 on StackOverflow: https://stackoverflow.com/a/66674679.
    """
    delta = {mul_xor_hash(arr2[0])}
    for i in range(1, arr2.shape[0]):
        delta.add(mul_xor_hash(arr2[i]))
    n = 0
    for i in range(arr1.shape[0]):
        if mul_xor_hash(arr1[i]) not in delta:
            n += 1
    result = np.empty((n, arr1.shape[-1]), dtype=arr1.dtype)
    j = 0
    for i in range(arr1.shape[0]):
        if mul_xor_hash(arr1[i]) not in delta:
            result[j] = arr1[i]
            j += 1
    return result

def int_vec_inside(c: NDArray[np.int32]) -> NDArray[np.int32]:
    """Integer vectors inside the cell `c` whose elements are integers.

    Parameters
    ----------
    c : (3, 3) array of ints
        A matrix whose columns are integer cell vectors.
    
    Returns
    -------
    v_int : (3, ...) array of ints
        Its columns are vectors satisfying `v = c @ k`, where `k[0]`, `k[1]`, `k[2]` $\in [0, 1)$.
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
    d = hnf_rational(np.hstack((m1, m2)), max_divisor=max_divisor)[:,:3]
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

def all_hnf(det: int) -> NDArray[np.int32]:
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

def hnf(m: NDArray[np.int32]) -> NDArray[np.int32]:
    """Return the Hermite normal form (HNF) of square integer matrix `m`.

    Parameters
    ----------
    m : (M, N) array of ints
        The integer matrix to reduce, with positive determinant and M <= N (the function will not check this).
    
    Returns
    -------
    h : (M, N) array of ints
        The column-style Hermite normal form of `m`.
    """
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

def hnf_int(m: NDArray[np.int32], return_q: bool = True) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Decompose square integer matrix `m` into product of HNF matrix `h` and unimodular matrix `q`.

    Parameters
    ----------
    m : (M, N) array of ints
        The integer matrix to decompose, with positive determinant and M <= N.
    return_q : bool, optional
        Whether to return the unimodular matrix `q`. Default is True.
    
    Returns
    -------
    h : (M, N) array of ints
        The column-style Hermite normal form of `m`.
    q : (N, N) array of ints
        The unimodular matrix satisfying `m` = `h @ q`.
    """
    if not m.dtype == int: raise TypeError(f"Input matrix must be integer:\n{m}")
    if not la.matrix_rank(m, tol=1e-6) == m.shape[0]: raise ValueError(f"Input matrix must be full-row-rank:\n{m}")
    h = hnf(m)
    if not return_q or m.shape[1] != m.shape[0]: return h
    else: return h, (la.inv(h) @ m).round().astype(int)

def hnf_rational(m: ArrayLike, max_divisor = 100, tol=1e-3) -> NDArray[np.float64]:
    """The Hermite normal form (HNF) of full-row-rank rational matrix `m` (not necessarily square or integer).
    
    Parameters
    ----------
    m : (M, N) array_like, M <= N
        The full-row-rank rational matrix to reduce.
    max_divisor : int
        A positive integer. The least common multiple of all divisors in `m` should not be greater than `max_divisor`.
    tol : float
        The tolerance for rational approximation.
    
    Returns
    -------
    h : (M, N) array
        The HNF of `m` obtained via elementary column operations over integers.
    """
    for divisor in range(1, max_divisor+1):
        if (np.abs((m * divisor) % 1.0) <= tol * divisor).all(): break
        elif divisor == max_divisor: raise ValueError(f"Input matrix must be rational with divisor LCM <= {max_divisor}:\n{m}")
    h = (m * divisor).round().astype(int)
    n_row, n_col = h.shape
    if not (n_row <= n_col and la.matrix_rank(h, tol=tol) == n_row): raise ValueError(f"Input matrix must be full-row-rank:\n{m}")
    return hnf_int(h, return_q=False) / divisor

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
    vv = v.copy()
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

def cell_reduce(c: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
    """Shorten the cell vectors of `c`.

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
    cc = c.copy()
    while (cc != c0).any():
        c0 = cc.copy()
        cc = cc[:,np.argsort(la.norm(cc, axis=0))]
        cc[:,2] = vector_reduce(cc[:,2], cc[:,0:2])
    if la.det(cc) < 0: cc = -cc
    q = (la.inv(c) @ cc).round().astype(int)
    return cc, q

def voigt_to_tensor(voigt_matrix, cryst=None, tol=1e-3, verbose=True):
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
    tensor = np.ones((9,9)) * np.inf       # XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
    voigt_ind = [0,4,8,5,6,1]              # XX,YY,ZZ,YZ,ZX,XY
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
