"""
Enumerate SLMs and CSMs.
"""

from .utilities import *
from scipy.optimize import linear_sum_assignment, minimize, basinhopping
from scipy.stats.qmc import Sobol

np.set_printoptions(suppress=True)
Constraint = namedtuple("Constraint", ["enforce", "prevent"])
IMPOSSIBLE_DISTANCE = 1e10
DELTA_K = np.mgrid[-1:1,-1:1,-1:1].reshape(3,1,1,-1)
NEIGHBOR_K = np.mgrid[-1:2,-1:2,-1:2].reshape(3,-1).T[np.arange(27)!=13]

def sigma_weight(sigma):
    """Return the weight of a given (3,) array of singular values, which is proportional to the standard measure on R^9 
    """
    return np.abs(np.diff((np.array(sigma)[...,[[0,1,2],[1,2,0]]])**2, axis=-2).squeeze(axis=-2)).prod(axis=-1)

def strain_energy_function(elasticA, elasticB):
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

def pct_distance(c, pA, pB, p, ks, weights=None, l=2, min_t0=True, return_t0=False):
    """Return the shuffle distance of a PCT.
    
    Parameters
    ----------
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
    min_t0 : bool, optional
        Set to True to minimize the shuffle distance by translating the final structure. Default is True.
    return_t0 : bool, optional
        Whether to return the best overall translation if `min_t0` is True. Default is False.
    
    Returns
    -------
    distance : float
        The shuffle distance.
    t0 : (3, 1) array
        The best overall translation, reshaped as a 3x1 matrix.
    """
    if ks.shape != (3,len(p)): raise ValueError("'p' and 'ks' must have the same number of atoms.")
    if not min_t0: return np.average(((c @ (pB[:,p] + ks - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l)
    if l == 2:
        t0 = -np.average(pB[:,p] + ks - pA, axis=1, weights=weights, keepdims=True)
        d = np.average(((c @ (pB[:,p] + ks + t0 - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l)
    else:
        res = minimize(lambda t: np.average(la.norm(c @ (pB[:,p] + ks + t.reshape(3,1) - pA), axis=0)**l, weights=weights),
                        -np.average(pB[:,p] + ks - pA, axis=1, weights=weights), method='SLSQP', options={'disp': False})
        d = res.fun ** (1/l)
        t0 = res.x
    if return_t0: return d, t0
    else: return d

def csm_distance(crystA, crystB, slm, p, ks, weight_func=None, l=2, min_t0=True, return_t0=False):
    """Return the shuffle distance of a CSM.
    """
    crystA_sup, crystB_sup, c_sup_half, _, _ = create_common_supercell(crystA, crystB, slm)
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    species = crystA_sup[1]
    weights = [weight_func[s] for s in species] if weight_func else None
    return pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weights, l=l, min_t0=min_t0, return_t0=return_t0)

def standardize_imt(slm: Union[SLM, NDArray[np.int32]], gA: NDArray[np.int32], gB: NDArray[np.int32]) -> SLM:
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
    """
    hA, hB, q = slm
    cl = np.dot((gB @ hB) @ q, la.inv(gA @ hA)).transpose((2,0,1,3)).reshape(-1,9).round(4)
    ind = np.arange(cl.shape[0])
    for j in range(9):
        ind = ind[cl[ind,j] == np.min(cl[ind,j])]
        if len(ind) == 1: break
    iA, iB = np.unravel_index(ind, (gA.shape[0], gB.shape[0]))
    i = np.lexsort(np.array([gA[iA], gB[iB]]).transpose((0,2,3,1)).reshape(18,-1))[0]
    hAA, qA = hnf_square(gA[iA[i]] @ hA)
    hBB, qB = hnf_square(gB[iB[i]] @ hB)
    slm0 = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return slm0

def enumerate_imt(
    crystA: Cryst, crystB: Cryst, mu: int, w_max: float,
    w: Callable = rmss, tol: float = 1e-3,
    max_iter: int = 1000, verbose: int = 1
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
    hA = all_hnfs(np.lcm(zA,zB) // zA * mu)
    hB = all_hnfs(np.lcm(zA,zB) // zB * mu)
    hA_inv = la.inv(hA)
    hB_inv = la.inv(hB)
    gA = get_pure_rotation(crystA, tol=tol)
    gB = get_pure_rotation(crystB, tol=tol)
    
    sobol6 = Sobol(6)
    if verbose >= 1:
        progress_bar = tqdm(total=max_iter, position=0, desc=f"\tmu={mu:-2d}", ncols=60, mininterval=0.5,
                            bar_format=f'{{desc}}: {{percentage:3.0f}}% |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]')
    slmlist = np.array([], dtype=int).reshape(-1,3,3,3)
    iter = 0
    
    while iter < max_iter:
        # generate a random `s0`
        rand6 = sobol6.random(1)
        u = cube_to_so3(rand6[0,0:3])
        vt = cube_to_so3(rand6[0,3:6])
        if w == rmss: quad_form = np.eye(3) / (3 * w_max**2)
        else: quad_form = w(None, return_quadratic_form=True, u=u, vt=vt) / (w_max * 1.0)
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
            slm = np.array(standardize_imt(slm, gA, gB), dtype=int)
            if not (slm == slmlist).all(axis=(1,2,3)).any():
                slmlist = np.concatenate((slmlist, slm.reshape(-1,3,3,3)))
                if iter >= max_iter / 10: max_iter *= 2          # current max_iter is not safe.
                iter = 0
                if verbose >= 1:
                    progress_bar.total = max_iter
                    progress_bar.n = 0
                    progress_bar.last_print_n = 0
                    progress_bar.refresh()
        iter = iter + 1
        if verbose >= 1: progress_bar.update(1)
    if verbose >= 1: progress_bar.close()
    return np.array(slmlist, dtype=int)

def punishment_matrix(z, prevent, l=2):
    m = np.zeros((z,z))
    if not prevent: return m
    prevent_arr = np.array(list(prevent), dtype=int).reshape(-1,2)
    m[prevent_arr[:,0],prevent_arr[:,1]] = IMPOSSIBLE_DISTANCE**l
    return m

def optimize_pct_fixed(c, species, pA, pB, constraint=Constraint(set(),set()), weight_func=None, l=2):
    """
    c : (3, 3) array
        The base matrix of the shuffle lattice.
    species : (Z, ) array of str
        The species of each atom.
    pA, pB : (3, Z) array
        The fractional coordinates of the atoms in the initial and final structures, respectively.
    constraint : Constraint, optional
        The constraint on the permutation. Default is an empty constraint.
    weight_func : dict, optional
        The weight function, with keys as atomic species. Default is None, which means all atoms have the same weight.
    l : float, optional
        The l-norm to be used for distance calculation, must not be less than 1. Default is 2.
    
    Returns
    -------
    d_hat : float
        The least shuffle distance.
    p : (Z, ) array of ints
        The permutation part of the PCT with the least shuffle distance.
    ks : (3, Z) array of ints
        The class-wise translation part of the PCT with the least shuffle distance.
    """
    z = len(species)
    if not (pA.shape[1] == z and pB.shape[1] == z): raise ValueError("Atom numbers do not match.")
    weights = [weight_func[s] for s in species] if weight_func else None
    
    p = np.ones(z, dtype=int) * z
    where_kij = np.ones((z,z), dtype=int) * z
    enforce = np.array(list(constraint.enforce), dtype=int).reshape(-1,2)
    if enforce.shape[0] > 0:
        p[enforce[:,0]] = enforce[:,1]
        diff_frac = (pB[:,enforce[:,1]] - pA[:,enforce[:,0]]).reshape(3,-1,1) % 1.0 + DELTA_K.reshape(3,1,8)
        norm_squared = np.sum(np.tensordot(c, diff_frac, axes=(1,0))**2, axis=0)
        where_kij[enforce[:,0], enforce[:,1]] = np.argmin(norm_squared, axis=1)
    punish = punishment_matrix(z, constraint.prevent, l=l)
    
    for s in np.unique(species):
        row_ind = species == s
        row_ind[enforce[:,0]] = False
        if not row_ind.any(): continue      # all assignments of species s are enforced
        zz = np.sum(row_ind)
        col_ind = species == s
        col_ind[enforce[:,1]] = False
        diff_frac = (pB[:,col_ind].reshape(3,1,-1,1) - pA[:,row_ind].reshape(3,-1,1,1)) % 1.0 + DELTA_K
        norm_squared = np.sum(np.tensordot(c, diff_frac, axes=(1,0))**2, axis=0)
        where_kij[np.ix_(row_ind, col_ind)] = np.argmin(norm_squared, axis=2)
        cost_matrix = np.min(norm_squared, axis=2) ** (l/2) + punish[np.ix_(row_ind, col_ind)]
        result = linear_sum_assignment(cost_matrix)[1]
        if cost_matrix[np.arange(zz),result].sum() >= IMPOSSIBLE_DISTANCE**l:         # all assignments of species s are prevented
            return IMPOSSIBLE_DISTANCE, None, None
        p[row_ind] = np.nonzero(col_ind)[0][result]
        
    ks = (-np.floor(pB[:,p] - pA) + DELTA_K[:,0,0,where_kij[np.arange(z),p]]).round().astype(int)
    d = pct_distance(c, pA, pB, p, ks, weights=weights, l=l, min_t0=False)
    return d, p, ks

def optimize_pct_local(c, species, pA, pB, t, constraint=Constraint(set(),set()), weight_func=None, l=2):
    t0 = t.reshape(3,1)
    n_iter = 0
    weights = [weight_func[s] for s in species] if weight_func else None
    while True:
        n_iter += 1
        dh, p, ks = optimize_pct_fixed(c, species, pA, pB+t0, constraint=constraint, weight_func=weight_func, l=l)
        if dh >= IMPOSSIBLE_DISTANCE: return IMPOSSIBLE_DISTANCE, None, None, None
        d, t0 = pct_distance(c, pA, pB, p, ks, weights=weights, l=l, return_t0=True)
        if dh - d < 1e-4: break
        if n_iter > 100: raise RecursionError("PCT optimization failed to converge. Please report this bug to wfc@pku.edu.cn if you see this message.")
    return d, p, ks, t0

def optimize_pct(crystA, crystB, slm, constraint=Constraint(set(),set()), weight_func=None, l=2, t_grid=64):
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = frac_cell(mA, mB)
    species = crystA_sup[1]
    z = species.shape[0]
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    reslist = [optimize_pct_local(c_sup_half, species, pA_sup, pB_sup, t, constraint=constraint, weight_func=weight_func, l=l)
                for t in Sobol(3).random(t_grid) @ f.T]
    ind = np.argmin([res[0] for res in reslist])
    d, p, ks, t0 = reslist[ind]
    return d, p, ks, t0

def pct_fill(crystA, crystB, slm, d_max, p, ks0=None, weight_func=None, l=2, warning_threshold=1000):
    """Returns all class-wise translations with permutation p that has d <= d_max, including ks0.
    """
    z = p.shape[0]
    if z == 1: return np.array([ks0], dtype=int)
    
    crystA_sup, crystB_sup, c_sup_half, _, _ = create_common_supercell(crystA, crystB, slm)
    species = crystA_sup[1]
    if not (crystB_sup[1][p] == species).all(): raise ValueError("Permuation does not preserve atom species.")
    vs = (crystB_sup[2].T[:,p] - crystA_sup[2].T).reshape(1,3,z)
    weights = [weight_func[s] for s in species] if weight_func else None
    if l == 2:
        def dl(ks):
            return np.average(la.norm(c_sup_half @ (vs + ks - np.average(vs + ks, axis=2, weights=weights, keepdims=True)), axis=1) ** l, axis=1, weights=weights)
    else:
        def dl(ks):
            t0 = minimize(lambda t: np.average(la.norm(c_sup_half @ (vs + ks - t.reshape(-1,3,1)), axis=1) ** l, axis=1, weights=weights).sum(),
                            np.average(vs + ks, axis=2, weights=weights).reshape(-1), method='SLSQP', options={'disp': False}).x
            return np.average(la.norm(c_sup_half @ (vs + ks - t0), axis=1) ** l, axis=1, weights=weights)
    if ks0 is None:
        d0, ks0, _ = optimize_ct(crystA, crystB, slm, p, weight_func=weight_func, l=l)
    else: d0 = dl(ks0.reshape(1,3,z)) ** (1/l)
    if d0 > d_max: raise ValueError("The shuffle distance of the given PCT already exceeds d_max.")

    # flood-filling ks
    dks = []
    for i in range(z):
        if z == 2 and i == 0: continue
        dk = np.zeros((3,z), dtype=int)
        for nk in NEIGHBOR_K:
            if i == 0: dk[:,1:] = -nk.reshape(3,1)
            else: dk[:,i] = nk
            dks.append(dk.copy())
    dks = np.array(dks, dtype=int)

    valid = np.array([ks0], dtype=int)
    visited = np.array([ks0], dtype=int)
    queue = ks0.reshape(1,3,z) + dks
    while queue.shape[0] > 0:
        visited = np.concatenate([visited, queue], axis=0)
        is_valid = dl(queue) <= d_max ** l
        new_valid = queue[is_valid,:,:]
        valid = np.concatenate([valid, new_valid], axis=0)
        if len(valid) > warning_threshold:
            print(f"\nWarning: Already filled {len(valid)} PCTs, which exceeds 'alarm_threshold'. We suggest you to use a smaller 'd_max'.")
            warning_threshold *= 2
        queue = np.unique((new_valid.reshape(-1,1,3,z) + dks.reshape(1,-1,3,z)).reshape(-1,3,z), axis=0)
        queue = setdiff2d(queue.reshape(-1,3*z), visited.reshape(-1,3*z)).reshape(-1,3,z)
    return valid

def is_compatible(p, constraint):
    enforce = np.array(list(constraint.enforce), dtype=int).reshape(-1,2)
    prevent = np.array(list(constraint.prevent), dtype=int).reshape(-1,2)
    return (p[enforce[:,0]] == enforce[:,1]).all() and (p[prevent[:,0]] != prevent[:,1]).all()

def murty_split(node, p):
    if not is_compatible(p, node): raise ValueError("The given permutation is not compatible with the given node.")
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

def cong_permutations(p, crystA, crystB, slm):
    crystA_sup, crystB_sup, _, mA, mB = create_common_supercell(crystA, crystB, slm)
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    tA = la.inv(mA) @ int_vec_inside(mA)
    tB = la.inv(mB) @ int_vec_inside(mB)
    lA = [np.nonzero((((pA_sup.reshape(3,1,-1) + tA[:,i].reshape(3,1,1) - pA_sup.reshape(3,-1,1))
                        .round(7) % 1.0) == 0).all(axis=0))[1] for i in range(tA.shape[1])]
    lB = [np.nonzero((((pB_sup.reshape(3,1,-1) + tB[:,j].reshape(3,1,1) - pB_sup.reshape(3,-1,1))
                        .round(7) % 1.0) == 0).all(axis=0))[1] for j in range(tB.shape[1])]
    return np.unique([kB[p[kA]] for kA in lA for kB in lB], axis=0)

def enumerate_pct(crystA, crystB, slm, d_max, weight_func=None, l=2, t_grid=16, incong=True, verbose=1):
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = frac_cell(mA, mB)
    species = crystA_sup[1]
    weights = [weight_func[s] for s in species] if weight_func else None
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T

    if verbose >= 1:
        print(f"Enumerating bottom PCTs with d <= {d_max:.4f} ...")
        prob_solved = 0
        prob_total = 1
        progress_bar = tqdm(total=prob_total, position=prob_solved, desc="\tnodes", ncols=50, mininterval=0.5,
                            bar_format=f'{{desc}}: {{n_fmt:>2s}}/{{total_fmt:>2s}} |{{bar}}| [elapsed {{elapsed:5s}}]')
    bottom_pcts = []
    node_stack = [Constraint(set(),set())]
    split_stack = []
    sobol3 = Sobol(3)
    while node_stack:
        node = node_stack.pop(0)
        split = False
        for i, p in enumerate(split_stack):
            if is_compatible(p, node):
                split = True
                for new_node in murty_split(node, p):
                    node_stack.append(new_node)
                    if verbose >= 1: prob_total += 1
                split_stack.pop(i)
                break
        if not split:
            for t in sobol3.random(t_grid) @ f.T:
                d, p, ks, _ = optimize_pct_local(c_sup_half, species, pA_sup, pB_sup, t, constraint=node, weight_func=weight_func, l=l)
                if d <= d_max:
                    bottom_pcts.append(zip_pct(p, ks))
                    for new_node in murty_split(node, p):
                        node_stack.append(new_node)
                        if verbose >= 1: prob_total += 1
                    if incong:
                        for p_cong in cong_permutations(p, crystA, crystB, slm):
                            split_stack.append(p_cong)
                if d <= d_max or d >= IMPOSSIBLE_DISTANCE: break
        if verbose >= 1:
            prob_solved += 1
            progress_bar.total = prob_total
            progress_bar.n = prob_solved
            progress_bar.last_print_n = prob_solved
            progress_bar.refresh()
    if verbose >= 1:
        progress_bar.close()
        print(f"\tFound {len(bottom_pcts)} {'incongruent bottom PCTs' if incong else 'bottom PCTs (not incongruent)'} with d <= {d_max}.")
        print(f"Filling non-bottom PCTs with d <= {d_max:.4f} ...")

    pctlist = []
    dlist = []
    if verbose >= 1:
        progress_bar = tqdm(total=len(bottom_pcts), position=0, desc="\tbottom PCTs", ncols=70, mininterval=0.5,
                        bar_format=f'{{desc}}: {{n_fmt:>2s}}/{{total_fmt:>2s}} |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]')
    for pct in bottom_pcts:
        for ks in pct_fill(crystA, crystB, slm, d_max, pct[:,0], ks0=pct[:,1:].T, weight_func=weight_func, l=l):
            pctlist.append(zip_pct(pct[:,0], ks))
            dlist.append(pct_distance(c_sup_half, pA_sup, pB_sup, pct[:,0], ks, weights=weights, l=l))
        if verbose >= 1: progress_bar.update(1)
    if verbose >= 1:
        progress_bar.close()
        print(f"\tFound {len(pctlist)} {'incongruent PCTs' if incong else 'PCTs (not incongruent)'} with d <= {d_max}.")
    return np.array(pctlist, dtype=int), np.array(dlist)

def optimize_ct_fixed(c, pA, pB, p, weights=None, l=2):
    """
    """
    k_grid = (-np.floor(pB[:,p] - pA).reshape(3,-1,1) + DELTA_K.reshape(3,1,-1)).round().astype(int)
    norm_squared = (np.tensordot(c, (pB[:,p] - pA).reshape(3,-1,1) + k_grid, axes=(1,0))**2).sum(axis=0)
    ks = k_grid[:,np.arange(len(p)),np.argmin(norm_squared, axis=1)]
    return pct_distance(c, pA, pB, p, ks, weights=weights, l=l, min_t0=False), ks

def optimize_ct_local(c, pA, pB, p, t, weights=None, l=2):
    t0 = t.reshape(3,1)
    n_iter = 0
    while True:
        n_iter += 1
        dh, ks = optimize_ct_fixed(c, pA, pB+t0, weights=weights, l=l)
        d, t0 = pct_distance(c, pA, pB, p, ks, weights=weights, l=l, return_t0=True)
        if dh - d < 1e-4: break
        if n_iter > 100: raise RecursionError("CT optimization failed to converge. Please report this bug to wfc@pku.edu.cn if you see this message.")
    return d, ks, t0

def optimize_ct(crystA, crystB, slm, p, weight_func=None, l=2, t_grid=64):
    """
    """
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = frac_cell(mA, mB)
    species = crystA_sup[1]
    if not (crystB_sup[1][p] == species).all(): raise ValueError("The given permutation does not preserve atomic species.")
    weights = [weight_func[s] for s in species] if weight_func else None
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    reslist = [optimize_ct_local(c_sup_half, pA_sup, pB_sup, t, weights=weights, l=l) for t in Sobol(3).random(t_grid) @ f.T]
    ind = np.argmin([res[0] for res in reslist])
    d, ks, t0 = reslist[ind]
    return d, ks, t0

def enumerate_rep_csm(crystA: Cryst, crystB: Cryst, mu_max: int, w_max: float,
    w: Callable[[NDArray[np.float64]], NDArray[np.float64]] = rmss, tol: float = 1e-3, max_iter: int = 1000,
    weight_func: Union[Dict[str, float], None] = None, l: int = 2, t_grid: int = 64, verbose: int = 1):
    """Enumerating all IMTs of multiplicity `mu` with `w` smaller than `w_max` and their representative CSMs.
    
    Parameters
    ----------
    crystA : Cryst
        The initial crystal structure.
    crystB : Cryst
        The final crystal structure.
    mu_max : int
        The maximum multiplicity of the IMTs to be enumerated.
    w_max : float
        The maximum value of the function `w`.
    w : Callable, optional
        The strain energy function to be used, usually obtained by `strain_energy_function()`; default is `rmss`.
    tol : float, optional
        The tolerance for `spglib` symmetry detection; default is 1e-3.
    max_iter : int, optional
        The maximum iteration number for IMT generation; default is 1000.
    weight_func : dict, optional
        A dictionary of atomic weights for each species; default is None, which means all atoms have the same weight.
    l : float, optional
        The type of norm to be used for distance calculation; default is 2.
    t_grid : int, optional
        The number of grid points for PCT optimization; default is 64.
    verbose : int, optional
        The level of verbosity; default is 1.
    
    Returns
    -------
    pct_arrs : list of (..., 4) arrays of ints
        The PCTs of representative CSMs, where `...` is the number of CSMs with multiplicity `mu`.
    slmlist : (N, 3, 3, 3) array of ints
        The IMTs enumerated.
    mulist : (N,) array of ints
        The multiplicities of the IMTs.
    wlist : (N,) array of floats
        The values of the function `w` for the IMTs.
    d0list : (N, ) array of floats
        The shuffle distances of the representative CSMs, i.e., the minimal shuffle distance of each SLM.
    """
    
    # enumerate SLMs
    if verbose: print(f"\nEnumerating incongruent SLMs for mu <= {mu_max} and {'rmss' if w == rmss else 'w'} <= {w_max:.4f} ...")
    slmlist = np.array([], dtype=int).reshape(0,3,3,3)
    for mu in range(1,mu_max+1):
        slmlist = np.concatenate([slmlist, enumerate_imt(crystA, crystB, mu, w_max, w=w, tol=tol, max_iter=max_iter, verbose=verbose)], axis=0)
    slmlist = np.array(slmlist)
    if verbose: print(f"A total of {len(slmlist):d} incongruent SLMs are enumerated:")
    if len(slmlist) == 0: raise Warning("No SLM is found. Try larger arguments for '--enumeration' or check if the input POSCARs are correct.")
    mulist = imt_multiplicity(crystA, crystB, slmlist)
    if verbose:
        print(f"\tmu  {' '.join(f'{i:5d}' for i in range(1,mu_max+1))}")
        print(f"\t#SLM{' '.join(f'{s:5d}' for s in np.bincount(mulist, minlength=mu_max+1)[1:])}")
    
    # exclude SLMs with the same deformation gradient
    _, ind = np.unique(deformation_gradient(crystA, crystB, slmlist).round(decimals=4), axis=0, return_index=True)
    slmlist = slmlist[ind]
    mulist = imt_multiplicity(crystA, crystB, slmlist)
    wlist = w(deformation_gradient(crystA, crystB, slmlist))
    
    # sort SLMs by 'mu' and then by 'w'
    ind = np.lexsort((wlist.round(decimals=4), mulist))
    slmlist = slmlist[ind]
    mulist = mulist[ind]
    wlist = wlist[ind]
    if verbose:
        print(f"Among them, a total of {len(slmlist):d} SLMs have distinct deformation gradients:")
        print(f"\tmu  {' '.join(f'{i:5d}' for i in range(1,mu_max+1))}")
        print(f"\t#SLM{' '.join(f'{s:5d}' for s in np.bincount(mulist, minlength=mu_max+1)[1:])}")

    # computing representative CSMs
    pct_arrs = [np.array("arr_mu.npy saves the PCTs of those CSMs with multiplicity mu", dtype=str)]
    d0list = np.zeros(len(slmlist))
    if verbose: print(f"Computing representative CSMs (the one with the lowest d among all CSMs with the same deformation gradient):")
    n_digit = np.floor(np.log10(np.bincount(mulist).max())).astype(int) + 1
    for mu in range(1,mu_max+1):
        n_csm = np.sum(mulist == mu)
        if n_csm == 0:
            pct_arrs.append(np.array(f"There is no CSM with multiplicity {mu}", dtype=str))
            continue
        pctlist = []
        if verbose: 
            progress_bar = tqdm(total=n_csm, position=0, desc=f"\r\tmu={mu:d}", ncols=57+2*n_digit, mininterval=0.5,
                    bar_format=f'{{desc}}: {{n_fmt:>{n_digit}s}}/{{total_fmt:>{n_digit}s}} |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]')
        for i in range(n_csm):
            d, p, ks, _ = optimize_pct(crystA, crystB, slmlist[mulist == mu][i], weight_func=weight_func, l=l, t_grid=t_grid)
            d0list[np.sum(mulist < mu) + i] = d
            pctlist.append(zip_pct(p, ks))
            if verbose: progress_bar.update(1)
        if verbose: progress_bar.close()
        pctlist = np.array(pctlist, dtype=int)
        
        # sort CSMs by 'w' and then by 'd'
        ind = np.lexsort((d0list[mulist == mu].round(decimals=4), wlist[mulist == mu].round(decimals=4)))
        slmlist[mulist == mu] = slmlist[mulist == mu][ind]
        mulist[mulist == mu] = mulist[mulist == mu][ind]
        wlist[mulist == mu] = wlist[mulist == mu][ind]
        d0list[mulist == mu] = d0list[mulist == mu][ind]
        pctlist = pctlist[ind]
        pct_arrs.append(pctlist)
    return pct_arrs, slmlist, mulist, wlist, d0list