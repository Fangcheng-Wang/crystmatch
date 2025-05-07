"""
Enumerate SLMs and CSMs.
"""

from .utilities import *
from scipy.optimize import brentq, linear_sum_assignment, minimize, basinhopping
from scipy.stats.qmc import Sobol

np.set_printoptions(suppress=True)
Constraint = namedtuple("Constraint", ["enforce", "prevent"])
IMPOSSIBLE_DISTANCE = 1e10
DELTA_K = np.mgrid[-1:2, -1:2, -1:2].reshape(3,1,1,-1)

def sigma_weight(sigma):
    """Return the weight of a given (3,) array of singular values, which is proportional to the standard measure on R^9 
    """
    return np.abs(np.diff((np.array(sigma)[...,[[0,1,2],[1,2,0]]])**2, axis=-2).squeeze(axis=-2)).prod(axis=-1)

def w_function(elasticA, elasticB):
    """Return the strain energy density function.
    """
    if elasticA.shape != (3,3,3,3) or elasticB.shape != (3,3,3,3):
        raise ValueError("elasticA and elasticB must be (3,3,3,3) arrays")
    def w(s, return_quadratic_form=False, u=None, vt=None):
        if return_quadratic_form:
            if u is None or vt is None: raise ValueError("u and v must be provided if return_quadratic_form is True")
            return 1/8 * (np.einsum('ijkl,mi,mj,nk,nl->mn', elasticA, vt, vt, vt, vt)
                        + np.einsum('ijkl,im,jm,kn,ln->mn', elasticB, u, u, u, u))
        else:
            s = np.array(s)
            u, sigma, vt = la.svd(s, compute_uv=True)
            return 1/2 * (np.einsum('ijkl,...mi,...mj,...nk,...nl,...m,...n->...', elasticA, vt, vt, vt, vt, sigma**0.5 - 1, sigma**0.5 - 1)
                        + np.einsum('ijkl,...im,...jm,...kn,...ln,...m,...n->...', elasticB, u, u, u, u, sigma**-0.5 - 1, sigma**-0.5 - 1))
    return w

def standardize_imt(slm: Union[SLM, NDArray[np.int32]], gA: NDArray[np.int32], gB: NDArray[np.int32]) -> tuple[SLM, int]:
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
    slm0 : slm
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
    slm0 = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return slm0, len(cl)

def enumerate_imt(
    crystA: Cryst, crystB: Cryst, mu: int, w_max: float,
    w = rmss, tol: float = 1e-3,
    max_iter: int = 2000, verbose: int = 1
) -> List[SLM]:
    """Enumerating all IMTs of multiplicity `mu` with `w` smaller than `w_max`.

    Parameters
    ----------
    crystA : Cryst
        The initial crystal structure, usually obtained by `load_poscar`.
    crystB : Cryst
        The final crystal structure, usually obtained by `load_poscar`.
    mu : int
        The multiplicity of SLMs to enumerate.
    w_max : float
        The maximum strain energy density, with the same units as `w`.
    w : callable, optional
        The strain energy density function. Default is `rmss`.
    tol : float, optional
        The tolerance for determining the pure rotation group of the crystal structures. Default is 1e-3.
    max_iter : int, optional
        The maximum number of consequtive iterations without finding any new SLMs. Default is 100.
    verbose : int, optional
        The level of verbosity. Default is 1.
    
    Returns
    -------
    
    """
    cA = crystA[0].T
    cB = crystB[0].T
    zA = crystA[2].shape[0]
    zB = crystB[2].shape[0]
    hA = all_hnf(np.lcm(zA,zB) // zA * mu)
    hB = all_hnf(np.lcm(zA,zB) // zB * mu)
    hA_inv = la.inv(hA)
    hB_inv = la.inv(hB)
    gA = get_pure_rotation(crystA, tol=tol)
    gB = get_pure_rotation(crystB, tol=tol)
    
    sobol6 = Sobol(6)
    if verbose >= 1:
        progress_bar = tqdm(total=max_iter, position=0, desc=f"\r\tmu={mu:d}", ncols=60, mininterval=0.5,
                            bar_format=f'{{desc}}: {{percentage:3.0f}}% |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]')
    slmlist = []
    iter = 0
    
    while iter < max_iter:
        # generate a random `s0`
        rand6 = sobol6.random(1)
        u = cube_to_so3(rand6[0,0:3])
        vt = cube_to_so3(rand6[0,3:6])
        if w == rmss: quad_form = np.eye(3) / (3 * w_max**2)
        else: quad_form = w(None, return_quadratic_form=True, u=u, vt=vt) / w_max
        # sampling `sigma - 1` in {x|quad_form(x,x)<1}
        size = 1000
        vec = np.random.normal(size=(size,3))
        mesh_ball = (np.random.uniform(low=0, high=1, size=size)**(1/3) / la.norm(vec, axis=1)).reshape(-1,1) * vec
        l = la.cholesky(quad_form)
        mesh_sigma = 1 + np.einsum('ij,kj->ki', la.inv(l.T), mesh_ball)
        max_weight = sigma_weight(mesh_sigma).max()
        sigma0 = None
        while sigma0 is None:
            temp = mesh_sigma[np.random.randint(mesh_sigma.shape[0])]
            if np.random.rand() > sigma_weight(temp) / max_weight: sigma0 = temp
        s0 = np.einsum('ij,j,jk->ik', u, sigma0, vt)
        
        # Round `s0`s to nearest integer matrix triplets
        q = np.transpose(np.dot(np.dot(hB_inv, la.inv(cB) @ s0 @ cA), hA), axes=(2,0,1,3)).round().astype(int)
        # Check determinants of `q`
        index1 = la.det(q.reshape(-1,3,3)).round().astype(int) == 1
        index1 = np.nonzero(index1)[0]
        index1 = np.unravel_index(index1, (hA.shape[0], hB.shape[0]))
        # Check strains
        index2 = np.nonzero(w(cB @ hB[index1[1]] @ q[index1] @ hA_inv[index1[0]] @ la.inv(cA)) <= w_max)[0]
        
        for i in index2:
            slm = (hA[index1[0][i]], hB[index1[1][i]], q[index1[0][i], index1[1][i]])
            slm, _ = standardize_imt(slm, gA, gB)
            repeated = False
            for j in range(len(slmlist)):
                if (slm[0] == slmlist[j][0]).all() and (slm[1] == slmlist[j][1]).all() and (slm[2] == slmlist[j][2]).all():
                    # `slm` is repeated.
                    repeated = True
                    break
            if not repeated:
                # `slm` is new.
                slmlist.append(slm)
                if iter >= max_iter / 10: max_iter *= 2          # current max_iter is not safe.
                iter = 0
                if verbose >= 1:
                    progress_bar.n = 0
                    progress_bar.total = max_iter
                    progress_bar.refresh()
        iter = iter + 1
        if verbose >= 1: progress_bar.update(1)
    if verbose >= 1: progress_bar.close()
    return np.array(slmlist, dtype=int)

def pct_distance(c, pA, pB, p, ks, weights=None, l=2):
    """
    c : (3, 3) array
        The lattice vectors of the crystal structure.
    pA, pB : (3, Z) array
        The fractional coordinates of the atoms in the initial and final structures, respectively.
    p : (Z, ) array of ints
        The permutation of the atoms.
    ks : (3, Z) array of floats
        The class-wise translations (fractional coordinates) of the atoms in `pB`.
    weights : (Z, ) array of floats, optional
        The weights of each atom. Default is None, which means all atoms have the same weight.
    l : float, optional
        The l-norm to be used for distance calculation, must not be less than 1. Default is 2.
    
    Returns
    -------
    distance : float
        The shuffle distance.
    """
    if l == 2:
        t = np.average(pA - pB[:,p] - ks, axis=1, weights=weights, keepdims=True)
        return np.average(la.norm(c @ (pA - pB[:,p] - ks - t), axis=0)**l, weights=weights) ** (1/l), t
    else:
        res = minimize(lambda t: np.average(la.norm(c @ (pA - pB[:,p] - ks - t), axis=0)**l, weights=weights),
                        np.average(pA - pB[:,p] - ks, axis=1, weights=weights, keepdims=True), method='BFGS', options={'disp': False})
        return res.fun ** (1/l), res.x
        

def minimize_rmsd_tp(
    c: NDArray[np.float64],
    pA: NDArray[np.float64],
    pB: NDArray[np.float64]
) -> Tuple[float, NDArray[np.int32]]:
    """Minimize the RMSD with fixed SLM, overall translation, permutation, and variable class-wise translations.
    
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
        The class-wise translations (fractional coordinates) of the atoms in `pB` that minimizes the RMSD.
    """
    raise NotImplementedError("This function is not implemented yet.")
    assert pA.shape[1] == pB.shape[1]
    shift = np.mgrid[-1:2, -1:2, -1:2].reshape(3,1,-1)
    relative_position = np.tensordot(c, (pA.reshape(3,-1,1) - pB.reshape(3,-1,1) - shift), axes=[1,0])
    d2_matrix = np.sum(relative_position ** 2, axis=0)
    ks = shift[:,0,np.argmin(d2_matrix, axis=1)]
    rmsd = la.norm(c @ (pB + ks - pA)) / np.sqrt(pA.shape[1])
    return rmsd, ks

def punishment_matrix(z, prevent, l=2):
    m = np.zeros((z,z))
    if not prevent: return m
    prevent_arr = np.array(list(prevent), dtype=int).reshape(-1,2)
    m[prevent_arr[:,0],prevent_arr[:,1]] = IMPOSSIBLE_DISTANCE**l
    return m

def shuffle_distance(c, species, pA, pB, constraint=Constraint(set(),set()), weight_function=None, l=2):
    """
    c : (3, 3) array
        The base matrix of the shuffle lattice.
    species : (Z, ) array of str
        The species of each atom.
    pA, pB : (3, Z) array
        The fractional coordinates of the atoms in the initial and final structures, respectively.
    constraint : Constraint, optional
        The constraint on the permutation. Default is an empty constraint.
    weight_function : dict, optional
        The weight for each atomic species. Default is None, which means all atoms have the same weight.
    l : float, optional
        The l-norm to be used for distance calculation, must not be less than 1. Default is 2.
    
    Returns
    -------
    d_hat : float
        The least shuffle distance.
    p : (Z, ) array of ints
        The permutation part of the PCT with the least shuffle distance.
    ks : (3, Z) array of floats
        The class-wise translation part of the PCT with the least shuffle distance.
    """
    if pA.shape[1] == pB.shape[1]: z = pA.shape[1]
    else: ValueError("Atom numbers do not match. Please report this bug to wfc@pku.edu.cn if you see this message.")

    p = np.ones(z, dtype=int) * -1
    enforce = np.array(list(constraint.enforce), dtype=int).reshape(-1,2)
    if enforce.shape[0] > 0: p[enforce[:,0]] = enforce[:,1]

    relative_position = (c @ pA).reshape(3,-1,1,1) - (c @ pB).reshape(3,1,-1,1) - np.tensordot(c, DELTA_K, axes=(1,0))
    cost_tensor = la.norm(relative_position, axis=0) ** l
    where_kij = np.argmin(cost_tensor, axis=2)
    cost_matrix = np.min(cost_tensor, axis=2) + punishment_matrix(z, constraint.prevent)
        
    for s in np.unique(species):
        row_ind = (species == s) * (np.arange(z)[:,None] != enforce[:,0][None,:]).all(axis=1)
        if not row_ind.any(): continue      # all assignments of species s are enforced
        col_ind = (species == s) * (np.arange(z)[:,None] != enforce[:,1][None,:]).all(axis=1)
        result = linear_sum_assignment(cost_matrix[row_ind,:][:,col_ind])[1]
        if cost_matrix[row_ind,:][:,col_ind][np.arange(len(result)),result].sum() >= IMPOSSIBLE_DISTANCE**2:         # all assignments of species s are prevented
            return IMPOSSIBLE_DISTANCE, None, None
        p[row_ind] = np.nonzero(col_ind)[0][result]
        
    d_hat = np.average(cost_matrix[np.arange(z),p], weights=weight_function[species]) ** (1/l)
    ks = DELTA_K[:,0,0,where_kij[np.arange(z),p]]
    return d_hat, p, ks

def optimize_pct(crystA, crystB, slm, constraint=Constraint(set(),set()), weight_function=None, l=2):
    crystA_sup, crystB_sup, c_sup_half, f_translate = create_common_supercell(crystA, crystB, slm)
    species = crystA_sup[1]
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    z = pA_sup.shape[1]
    def d_hat(t):
        return shuffle_distance(c_sup_half, species, pA_sup, pB_sup + f_translate @ t.reshape(3,1),
                                constraint=constraint, weight_function=weight_function, l=l)[0]
    t0 = basinhopping(d_hat, [0.5, 0.5, 0.5], T=0.05, niter=50+5*int(z**0.5), minimizer_kwargs={"method": "BFGS"}).x
    _, p, ks = shuffle_distance(c_sup_half, species, pA_sup, pB_sup + f_translate @ t0.reshape(3,1),
                                constraint=constraint, weight_function=weight_function, l=l)
    if (punishment_matrix(z, constraint.prevent)[np.arange(z),p] != 0).any():
        return IMPOSSIBLE_DISTANCE, None, None
    d, t0 = pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weight_function[species], l=l)
    return d, p, ks, t0

def pct_fill(crystA, crystB, slm, p, ks0, d_max, weight_function=None, l=2, warning_threshold=1000):
    crystA_sup, crystB_sup, c_sup_half, _ = create_common_supercell(crystA, crystB, slm)
    if not (crystA_sup[1] == crystB_sup[1][p]).all():
        ValueError("Permuation does not preserve atom species.")
    z = p.shape[0]
    weights = weight_function[crystA_sup[1]]
    vs = (crystB_sup[2].T[:,p] - crystA_sup[2].T).reshape(1,3,z)
    def dl(ks):
        if l == 2:
            return np.average(la.norm(c_sup_half @ (vs + ks - np.average(vs + ks, axis=2, weights=weights, keepdims=True)), axis=1) ** l, axis=1, weights=weights)
        else:
            t0 = minimize(lambda t: np.average(la.norm(c_sup_half @ (vs + ks - t), axis=1) ** l, axis=1, weights=weights).sum(),
                            np.average(vs + ks, axis=2, weights=weights, keepdims=True), method='BFGS', options={'disp': False}).x
            return np.average(la.norm(c_sup_half @ (vs + ks - t0), axis=1) ** l, axis=1, weights=weights)
    if not (dl(ks0.reshape(1,3,z)) <= d_max ** l).all():
        ValueError("The initial PCT already exceeds the maximum RMSD.")

    # flood-filling ks
    dks = np.eye(3*z, dtype=int).reshape(3*z,z,3).transpose(0,2,1)
    dks = np.concatenate([dks[3:], -dks[3:]], axis=0)
    valid = np.array([ks0], dtype=int)
    visited = np.array([ks0], dtype=int)
    queue = ks0.reshape(1,3,z) + dks
    while queue.shape[0] > 0:
        visited = np.concatenate([visited, queue], axis=0)
        is_valid = dl(queue) <= d_max ** l
        valid = np.concatenate([valid, queue[is_valid,:,:]], axis=0)
        if len(valid) > warning_threshold:
            print(f"\nWarning: Already filled {len(valid)} PCTs, which exceeds 'alarm_threshold'. We suggest you to use a smaller 'd_max'.")
            warning_threshold *= 2
        queue = np.unique((queue[is_valid,:,:].reshape(-1,1,3,z) + dks.reshape(1,-1,3,z)).reshape(-1,3,z), axis=0)
        is_nonvisited = (1 - (queue.reshape(-1,1,3,z) == visited.reshape(1,-1,3,z)).all(axis=(2,3)).any(axis=1)).astype(bool)
        queue = queue[is_nonvisited]
    return valid

def murty_split(node, p):
    new_nodes = []
    tuples = [(i,int(p[i])) for i in range(len(p))]
    for i in range(len(p) - 1):
        enforce = node.enforce | set(tuples[:i])
        prevent = node.prevent | {tuples[i]}
        if len({e[0] for e in enforce}) != len(enforce): continue       # conflict enforces in crystA
        if len({e[1] for e in enforce}) != len(enforce): continue       # conflict enforces in crystB
        if enforce & prevent: continue                                  # conflict enforce and prevent
        new_nodes.append(Constraint(enforce, prevent))
    return new_nodes

def enumerate_pct(crystA, crystB, slm, d_max, weight_function=None, l=2, verbose=1):
    mypcts = []
    stack = [Constraint(set(),set())]
    p_count = 0
    while stack:
        node = stack.pop()
        rmsd, p, ks0, _ = optimize_pct(crystA, crystB, slm, constraint=node, weight_function=weight_function, l=l)
        if rmsd <= d_max:
            p_count += 1
            for ks in pct_fill(crystA, crystB, slm, p, ks0, d_max, weight_function=weight_function, l=l):
                mypcts.append(zip_pct(p, ks))
            for new_node in murty_split(node, p):
                stack.append(new_node)
    if verbose >= 1:
        print(f"Found {len(mypcts)} PCTs (with {p_count} unique permutations) with RMSD <= {d_max}.")
    return np.array(mypcts, dtype=int)