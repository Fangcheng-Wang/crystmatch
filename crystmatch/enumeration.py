"""
Enumerate SLMs and CSMs.
"""

from .utilities import *
from scipy.optimize import linear_sum_assignment, minimize, OptimizeResult
from scipy.stats.qmc import Sobol

np.set_printoptions(suppress=True)
Constraint = namedtuple("Constraint", ["enforce", "prevent"])
IMPOSSIBLE_DISTANCE = 1e10
DELTA_K = np.mgrid[-1:1,-1:1,-1:1].reshape(3,1,1,-1)

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
    hAA, qA = hnf_int(gA[iA[i]] @ hA)
    hBB, qB = hnf_int(gB[iB[i]] @ hB)
    slm0 = (hAA, hBB, qB @ q @ la.inv(qA).round().astype(int))
    return slm0

def enumerate_imt(
    crystA: Cryst, crystB: Cryst, mu: int, w_max: float,
    w = rmss, tol: float = 1e-3,
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

def pct_distance(c, pA, pB, p, ks, weights=None, l=2, fixed=False):
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
    fixed : bool, optional
        Set to True to fix and not return the overall translation. Default is False.
    
    Returns
    -------
    distance : float
        The shuffle distance.
    t0 : (3, ) array
        The best overall translation.
    """
    if ks.shape != (3,len(p)): raise ValueError("'p' and 'ks' must have the same number of atoms.")
    if fixed:
        return np.average(((c @ (pB[:,p] + ks - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l)
    if l == 2:
        t0 = -np.average(pB[:,p] + ks - pA, axis=1, weights=weights, keepdims=True)
        return np.average(((c @ (pB[:,p] + ks + t0 - pA))**2).sum(axis=0)**(l/2), weights=weights) ** (1/l), t0
    else:
        res = minimize(lambda t: np.average(la.norm(c @ (pB[:,p] + ks + t.reshape(3,1) - pA), axis=0)**l, weights=weights),
                        -np.average(pB[:,p] + ks - pA, axis=1, weights=weights), method='BFGS', options={'disp': False})
        return res.fun ** (1/l), res.x
        
def punishment_matrix(z, prevent, l=2):
    m = np.zeros((z,z))
    if not prevent: return m
    prevent_arr = np.array(list(prevent), dtype=int).reshape(-1,2)
    m[prevent_arr[:,0],prevent_arr[:,1]] = IMPOSSIBLE_DISTANCE**l
    return m

def optimize_pct_fixed(c, species, pA, pB, constraint=Constraint(set(),set()), weights=None, l=2):
    """
    c : (3, 3) array
        The base matrix of the shuffle lattice.
    species : (Z, ) array of str
        The species of each atom.
    pA, pB : (3, Z) array
        The fractional coordinates of the atoms in the initial and final structures, respectively.
    constraint : Constraint, optional
        The constraint on the permutation. Default is an empty constraint.
    weights : (Z, ) array of floats, optional
        The weights of each atom. Default is None, which means all atoms have the same weight.
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
    if pA.shape[1] == pB.shape[1]: z = pA.shape[1]
    else: ValueError("Atom numbers do not match. Please report this bug to wfc@pku.edu.cn if you see this message.")
    
    p = np.ones(z, dtype=int) * -1
    where_kij = np.ones((z,z), dtype=int) * -1
    enforce = np.array(list(constraint.enforce), dtype=int).reshape(-1,2)
    if enforce.shape[0] > 0:
        p[enforce[:,0]] = enforce[:,1]
        diff_frac = (pB[:,enforce[:,1]] - pA[:,enforce[:,0]]).reshape(3,-1,1) % 1.0 + DELTA_K.reshape(3,1,8)
        norm_squared = np.sum(np.tensordot(c, diff_frac, axes=(1,0))**2, axis=0)
        where_kij[enforce[:,0], enforce[:,1]] = np.argmin(norm_squared, axis=1)
    punish = punishment_matrix(z, constraint.prevent)
    
    for s in np.unique(species):
        row_ind = species == s
        row_ind[enforce[:,0]] = False
        if not row_ind.any(): continue      # all assignments of species s are enforced
        zz = np.sum(row_ind)
        col_ind = species == s
        col_ind[enforce[:,1]] = False
        diff_frac = (pB[:,col_ind].reshape(3,1,-1,1) - pA[:,row_ind].reshape(3,-1,1,1)) % 1.0 + DELTA_K
        temp = np.zeros((3,zz*zz*8))
        np.dot(c, diff_frac.reshape(3,-1), out=temp)
        norm_squared = np.sum(temp**2, axis=0).reshape(zz,zz,8)
        where_kij[np.ix_(row_ind, col_ind)] = np.argmin(norm_squared, axis=2)
        cost_matrix = np.min(norm_squared, axis=2) ** (l/2) + punish[np.ix_(row_ind, col_ind)]
        result = linear_sum_assignment(cost_matrix)[1]
        if cost_matrix[np.arange(zz),result].sum() >= IMPOSSIBLE_DISTANCE**l:         # all assignments of species s are prevented
            return IMPOSSIBLE_DISTANCE, None, None
        p[row_ind] = np.nonzero(col_ind)[0][result]
        
    ks = (- np.floor(pB[:,p] - pA) + DELTA_K[:,0,0,where_kij[np.arange(z),p]]).round().astype(int)
    return pct_distance(c, pA, pB, p, ks, weights=weights, l=l, fixed=True), p, ks

def optimize_pct(crystA, crystB, slm, constraint=Constraint(set(),set()), weights=None, l=2, t_grid=64):
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = f_translate(mA, mB)
    species = crystA_sup[1]
    z = species.shape[0]
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    def d_hat(t):
        return optimize_pct_fixed(c_sup_half, species, pA_sup, pB_sup + f @ t.reshape(3,1),
                                constraint=constraint, weights=weights, l=l)[0]
    temp = np.inf
    for t in Sobol(3).random(t_grid):
        res = minimize(d_hat, t, method='SLSQP', options={'disp': False, 'ftol': 1e-7})
        if res.fun < temp:
            temp = res.fun
            t0 = res.x
    dd, p, ks = optimize_pct_fixed(c_sup_half, species, pA_sup, pB_sup + f @ t0.reshape(3,1),
                                constraint=constraint, weights=weights, l=l)
    if (punishment_matrix(z, constraint.prevent)[np.arange(z),p] != 0).any(): return IMPOSSIBLE_DISTANCE, None, None, None
    d, t0 = pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weights, l=l)
    if np.abs(d - dd) > 1e-5:
        print(f"Grid-minimization ({dd:.5f}) and analytical solution ({d:.5f}) disagree.")
        print("If this continues to happen, even if a larger 't_grid' is used, please report this bug to wfc@pku.edu.cn.")
    return d, p, ks, t0

def pct_fill(crystA, crystB, slm, p, ks0, d_max, weights=None, l=2, warning_threshold=1000):
    """Returns all class-wise translations with permutation p that has d <= d_max, including ks0.
    """
    crystA_sup, crystB_sup, c_sup_half, _, _ = create_common_supercell(crystA, crystB, slm)
    if not (crystA_sup[1] == crystB_sup[1][p]).all():
        ValueError("Permuation does not preserve atom species.")
    z = p.shape[0]
    vs = (crystB_sup[2].T[:,p] - crystA_sup[2].T).reshape(1,3,z)
    if l == 2:
        def dl(ks):
            return np.average(la.norm(c_sup_half @ (vs + ks - np.average(vs + ks, axis=2, weights=weights, keepdims=True)), axis=1) ** l, axis=1, weights=weights)
    else:
        def dl(ks):
            t0 = minimize(lambda t: np.average(la.norm(c_sup_half @ (vs + ks - t.reshape(-1,3,1)), axis=1) ** l, axis=1, weights=weights).sum(),
                            np.average(vs + ks, axis=2, weights=weights).reshape(-1), method='BFGS', options={'disp': False}).x
            return np.average(la.norm(c_sup_half @ (vs + ks - t0), axis=1) ** l, axis=1, weights=weights)
    if not (dl(ks0.reshape(1,3,z)) <= d_max ** l).all():
        ValueError("The shuffle distance of the given PCT already exceeds d_max.")

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

def enumerate_pct(crystA, crystB, slm, d_max, weights=None, l=2, t_grid=16, incong=True, verbose=1):
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = f_translate(mA, mB)
    species = crystA_sup[1]
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T

    if verbose >= 1:
        prob_solved = 0
        prob_total = 1
        progress_bar = tqdm(total=prob_total, position=prob_solved, desc=f"\r\titerations", ncols=60, mininterval=0.5,
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
            def d_hat(t):
                return optimize_pct_fixed(c_sup_half, species, pA_sup, pB_sup + f @ t.reshape(3,1),
                                        constraint=node, weights=weights, l=l)[0]
            for t in sobol3.random(t_grid):
                res = minimize(d_hat, t, method='SLSQP', options={'disp': False})
                if res.fun <= d_max:
                    t0 = res.x
                    _, p, ks = optimize_pct_fixed(c_sup_half, species, pA_sup, pB_sup + f @ t0.reshape(3,1),
                                            constraint=node, weights=weights, l=l)
                    bottom_pcts.append(zip_pct(p, ks))
                    for new_node in murty_split(node, p):
                        node_stack.append(new_node)
                        if verbose >= 1: prob_total += 1
                    if incong:
                        for p_cong in cong_permutations(p, crystA, crystB, slm):
                            split_stack.append(p_cong)
                    break
                elif res.fun >= IMPOSSIBLE_DISTANCE:
                    break
        if verbose >= 1:
            prob_solved += 1
            progress_bar.total = prob_total
            progress_bar.n = prob_solved
            progress_bar.last_print_n = prob_solved
            progress_bar.refresh()
    progress_bar.close()
    if verbose >= 1: print(f"Found {len(bottom_pcts)} {'incongruent permutations' if incong else 'permutations (not incongruent)'} with d <= {d_max}.")

    pctlist = []
    dlist = []
    for pct in bottom_pcts:
        for ks in pct_fill(crystA, crystB, slm, pct[:,0], pct[:,1:].T, d_max, weights=weights, l=l):
            pctlist.append(zip_pct(pct[:,0], ks))
            dlist.append(pct_distance(c_sup_half, pA_sup, pB_sup, pct[:,0], ks, weights=weights, l=l)[0])
    if verbose >= 1: print(f"Found {len(pctlist)} {'incongruent PCTs' if incong else 'PCTs (not incongruent)'} with d <= {d_max}.")
    return np.array(pctlist, dtype=int), np.array(dlist)

def optimize_ct_fixed(c, pA, pB, p, weights=None, l=2):
    """
    """
    k_grid = (- np.floor(pB[:,p] - pA).reshape(3,-1,1) + DELTA_K.reshape(3,1,-1)).round().astype(int)
    norm_squared = (np.tensordot(c, (pB[:,p] - pA).reshape(3,-1,1) + k_grid, axes=(1,0))**2).sum(axis=0)
    ks = k_grid[:,np.arange(len(p)),np.argmin(norm_squared, axis=1)]
    return pct_distance(c, pA, pB, p, ks, weights=weights, l=l, fixed=True), ks

def optimize_ct(crystA, crystB, slm, p, weights=None, l=2, t_grid=64):
    """
    """
    crystA_sup, crystB_sup, c_sup_half, mA, mB = create_common_supercell(crystA, crystB, slm)
    f = f_translate(mA, mB)
    if not (crystA_sup[1] == crystB_sup[1][p]).all():
        raise ValueError("The given permutation does not preserve atomic species.")
    pA_sup = crystA_sup[2].T
    pB_sup = crystB_sup[2].T
    def d_hat(t):
        return optimize_ct_fixed(c_sup_half, pA_sup, pB_sup + f @ t.reshape(3,1), p, weights=weights, l=l)[0]
    temp = np.inf
    for t in Sobol(3).random(t_grid):
        res = minimize(d_hat, t, method='SLSQP', options={'disp': False, 'ftol': 1e-7})
        if res.fun < temp:
            temp = res.fun
            t0 = res.x
    _, ks = optimize_ct_fixed(c_sup_half, pA_sup, pB_sup + f @ t0.reshape(3,1), p, weights=weights, l=l)
    d, t0 = pct_distance(c_sup_half, pA_sup, pB_sup, p, ks, weights=weights, l=l)
    return d, ks, t0


