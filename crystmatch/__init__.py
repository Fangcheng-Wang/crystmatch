from .utilities import *
from .enumeration import *
from .analysis import *
import sys
import argparse

__name__ = "crystmatch"
__version__ = "1.1.4"
__author__ = "Fang-Cheng Wang"
__email__ = "wfc@pku.edu.cn"
__description__ = 'Enumerating and analyzing crystal-structure matches for solid-solid phase transitions.'
__url__ = 'https://fangcheng-wang.github.io/crystmatch/'
__epilog__ = 'The current version (v' + __version__ + ') may contain bugs. To get the latest version, please run:\
\n\n\t$ pip3 install --upgrade crystmatch\n\nWe also recommend you to see the documentation at:\
\n\n\t' + __url__ + '\n\nIf you use crystmatch in your research, please cite the following paper: \n\n\t\
[1] FC Wang, QJ Ye, YC Zhu, and XZ Li, Physical Review Letters 132, 086101 (2024) (https://arxiv.org/abs/2305.05278)\n\n\
You are also welcome to contact me at ' + __email__ + ' for any questions, feedbacks or comments.'

def enum_rep(crystA: Cryst, crystB: Cryst, mu_max: int, rmss_max: float, tol: float, accurate: bool):
    
    # enumerating SLMs
    print(f"\nEnumerating incongruent SLMs for mu <= {mu_max} and rmss <= {rmss_max:.4f}:")
    slmlist = []
    for mu in range(1,mu_max+1): slmlist = slmlist + enumerate_slm(crystA, crystB, mu, rmss_max, tol=tol)
    print(f"A total of {len(slmlist):d} incongruent SLMs are enumerated:")
    if len(slmlist) == 0: raise Warning("No SLM is found. Try larger arguments for '--enumeration' or check if the input POSCARs are correct.")
    slmlist = np.array(slmlist)
    mulist = multiplicity(crystA, crystB, slmlist)
    print(f"\tmu  {' '.join(f'{i:5d}' for i in range(1,mu_max+1))}")
    print(f"\t#SLM{' '.join(f'{s:5d}' for s in np.bincount(mulist, minlength=mu_max+1)[1:])}")
    
    _, ind = np.unique((slmlist[:,1,:,:] @ slmlist[:,2,:,:] @ la.inv(slmlist[:,0,:,:])).round(decimals=4), axis=0, return_index=True)
    slmlist = slmlist[ind]
    mulist = multiplicity(crystA, crystB, slmlist)
    rmsslist = rmss(sing_val(crystA, crystB, slmlist))
    ind = np.lexsort((rmsslist.round(decimals=4), mulist))
    slmlist = slmlist[ind]
    mulist = mulist[ind]
    rmsslist = rmsslist[ind]
    print(f"Among them, a total of {len(slmlist):d} SLMs has distinct deformation gradients:")
    print(f"\tmu  {' '.join(f'{i:5d}' for i in range(1,mu_max+1))}")
    print(f"\t#SLM{' '.join(f'{s:5d}' for s in np.bincount(mulist, minlength=mu_max+1)[1:])}")

    # computing representative CSMs
    csm_bins = [np.array("arr_mu.npy saves the shuffles of those CSMs with multiplicity mu", dtype=str)]
    rmsdlist = np.inf * np.ones(len(slmlist))
    print(f"Minimizing RMSD to obtain the representative CSM of each deformation gradient:")
    zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
    n_digit = np.floor(np.log10(np.bincount(mulist).max())).astype(int) + 1
    for mu in range(1,mu_max+1):
        n_csm = np.sum(mulist == mu)
        if n_csm == 0:
            csm_bins.append(np.array(f"No CSM with multiplicity {mu}", dtype=str))
            continue
        shufflelist = np.zeros((n_csm, mu * zlcm, 4), dtype=int)
        for i in tqdm(range(n_csm), desc=f"\r\tmu={mu:d}", ncols=57+2*n_digit,
                    bar_format=f'{{desc}}: {{n_fmt:>{n_digit}s}}/{{total_fmt:>{n_digit}s}} |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]'):
            d, p, ks, _ = minimize_rmsd(crystA, crystB, slmlist[mulist == mu][i], accurate=accurate)
            rmsdlist[np.sum(mulist < mu) + i] = d
            shufflelist[i,:,0] = p
            shufflelist[i,:,1:] = ks.T
        ind = np.lexsort((rmsdlist[mulist == mu].round(decimals=4), rmsslist[mulist == mu].round(decimals=4)))
        slmlist[mulist == mu] = slmlist[mulist == mu][ind]
        mulist[mulist == mu] = mulist[mulist == mu][ind]
        rmsslist[mulist == mu] = rmsslist[mulist == mu][ind]
        rmsdlist[mulist == mu] = rmsdlist[mulist == mu][ind]
        shufflelist = shufflelist[ind]
        csm_bins.append(shufflelist)
    
    return csm_bins, slmlist, mulist, rmsslist, rmsdlist

def main():
    time0 = time()
    sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)
    sys.stderr = open(sys.stderr.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)

    parser = argparse.ArgumentParser(
        prog=__name__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='%(prog)s [-h] [-initial POSCAR_I] [-final POSCAR_F] [mode options] [output options]',
        description=__description__, epilog=__epilog__)
    parser.add_argument("-v", "--version", action='version', version="%(prog)s " + __version__)
    parser.add_argument("-E", "--enumeration", nargs=2, type=float, metavar=('MAX_MU','MAX_RMSS'),
                        help="use 'enumeration' mode, with MAX_MU and MAX_RMSS as the multiplicity and RMSS upper bounds")
    parser.add_argument("-R", "--read", type=str, metavar='CSM_LIST',
                        help="use 'read' mode, with CSMs loaded from CSM_LIST (an NPZ file generated by crystmatch)")
    parser.add_argument("-S", "--single", action='store_true', default=None,
                        help="use 'single-CSM' mode, with POSCARs specified by -I and -F to determine a single CSM")
    parser.add_argument("-I", "--initial", type=str, metavar='POSCAR_I',
                        help="POSCAR file of the initial crystal structure")
    parser.add_argument("-F", "--final", type=str, metavar='POSCAR_F',
                        help="POSCAR file of the final crystal structure")
    parser.add_argument("-a", "--accurate", action='store_true', help="use more accurate algorithm to optimize RMSD in 'enumeration' mode")
    parser.add_argument("-e", "--export", nargs='+', type=int, metavar=('index1', 'index2'),
                        help="export CSMs from CSM_LIST with the given indices")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-3, metavar='TOL',
                        help="tolerance for determining crystal symmetry; default is 1e-3")
    parser.add_argument("-o", "--orientation", nargs=12, type=float, metavar=('vix','viy','viz','vfx','vfy','vfz','wix','wiy','wiz','wfx','wfy','wfz'),
                        help="benchmark CSMs by the orientation relationship: vi||vf, wi||wf")
    parser.add_argument("-c", "--csv", action='store_true', help="whether to create CSV file")
    parser.add_argument("-p", "--plot", action='store_true', help="whether to create scatter plot")
    parser.add_argument("-u", "--uspfix", action='store_true', help="use USP-fixed manner [1] in '--orientation'")
    parser.add_argument("-A", nargs='*', type=str, metavar='', help="deprecated")
    parser.add_argument("-B", nargs='*', type=str, metavar='', help="deprecated")
    args = parser.parse_args()

    if args.A is not None or args.B is not None:
        raise DeprecationWarning("'-A -B' is deprecated since version 1.1.0, use '-I -F' instead. Run:\n\n\t$ crystmatch --help\n\nfor usage.")
        
    if sum(1 for x in [args.enumeration, args.read, args.single] if x is not None) >= 2:
        raise ValueError("Cannot choose more than one mode.")

    if sum(1 for x in [args.enumeration, args.read, args.single] if x is not None) == 0:
        print("No mode is specified (not recommended, see 'crystmatch -h' for common usage).")
        mode = input("Please enter a mode ('enumeration' / 'read' / 'single'): ")
        if mode in ["'enumeration'", 'enumeration', 'e', 'E']:
            mu_max = int(input("Enter the maximum multiplicity: "))
            rmss_max = float(input("Enter the maximum root-mean-square-strain: "))
            args.enumeration = (mu_max, rmss_max)
        elif mode in ["'read'", 'read', 'r', 'R']:
            args.read = input("Enter the path of the CSM_LIST file: ")
            args.export = np.array(input("Enter the indices of CSMs to be exported (separated by space), or leave blank for no export: ").split(), dtype=int)
            if args.export.shape[0] == 0: args.export = None
        elif mode in ["'single'", 'single', 'single-CSM', 's', 'S']:
            args.single = True
            print("Warning: The CSM will be uniquely determined by the initial and final POSCARs. I hope you know what you are doing.")
        else:
            raise ValueError(f"Invalid mode '{mode}'.")
    
    if args.enumeration is not None or (args.read is not None and args.orientation is not None):
        # in which case CSV and PLOT will always be generated
        args.csv = True
        args.plot = True
    
    # determining CSMs

    if args.enumeration is not None:
        print("\nMode: 'enumeration'")
        mu_max, rmss_max = args.enumeration
        mu_max = np.rint(mu_max).astype(int)
        if mu_max < 1: raise ValueError("Multiplicity should be a positive integer.")
        elif rmss_max <= 0: raise ValueError("Root-mean-square-strain should be a positive float.")
        if mu_max >= 8 or rmss_max > 0.5:
            print(f"Warning: Current MAX_MU = {mu_max:d} and MAX_RMSS = {rmss_max:.2f} may result in a large number of CSMs, which may take "
                    + "a long time to enumerate. If you are new to 'crystmatch', we recommend using MAX_MU <= 4 and MAX_RMSS = 0.4 first, "
                    + "and then increase MAX_MU and decrease MAX_RMSS gradually.")
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        crystA = load_poscar(args.initial, tol=args.tolerance)
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystB = load_poscar(args.final, tol=args.tolerance)
        check_chem_comp(crystA[1], crystB[1])
        job = f"m{mu_max:d}s{rmss_max:.2f}"
        csm_bins, slmlist, mulist, rmsslist, rmsdlist = enum_rep(crystA, crystB, mu_max, rmss_max, args.tolerance, args.accurate)
        table = np.array([np.arange(len(mulist)), np.arange(len(mulist)), mulist, rmsslist, rmsdlist], dtype=float).T

    elif args.read is not None:
        print(f"\nMode: 'read'")
        data = np.load(args.read, allow_pickle=True)
        table = data['table']
        slmlist = data['slmlist']
        job, crystA, crystB, mu_max, rmss_max, symprec = data['metadata']
        print(f"Metadata from '{args.read}':\n\ttotal number of CSMs: {table.shape[0]:d}\n\tmultipliticy <= {mu_max}"
                + f"\n\trms strain <= {rmss_max}\n\tsymmetry tolerance: {symprec}")
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))

    elif args.single is not None:
        print("\nMode: 'single-CSM'")
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        crystA_sup = load_poscar(args.initial, tol=args.tolerance, to_primitive=False)
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystB_sup = load_poscar(args.final, tol=args.tolerance, to_primitive=False)
        if not crystA_sup[1].shape[0] == crystB_sup[1].shape[0]: raise ValueError("\nNumbers of atoms in two crysts are not equal!")
        if not (crystA_sup[1] == crystB_sup[1]).all(): raise ValueError("\nAtom species in two crysts do not match!")
        crystA = load_poscar(args.initial, tol=args.tolerance, to_primitive=True, verbose=False)
        crystB = load_poscar(args.final, tol=args.tolerance, to_primitive=True, verbose=False)
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))
        job = "single"
        
        # single-CSM analysis
        print(f"\nGeometric properties of the CSM:")
        cA = crystA[0].T
        cB = crystB[0].T
        cA_sup = crystA_sup[0].T
        cB_sup = crystB_sup[0].T
        # computing SLM
        mA = (la.inv(cA) @ cA_sup).round().astype(int)
        mB = (la.inv(cB) @ cB_sup).round().astype(int)
        hA, qA = hnf_int(mA)
        hB, qB = hnf_int(mB)
        slmlist = [cong_slm((hA, hB, qB @ la.inv(qA).round().astype(int)),
                            get_pure_rotation(crystA, tol=args.tolerance), get_pure_rotation(crystB, tol=args.tolerance))[0]]
        mulist = multiplicity(crystA, crystB, slmlist)
        # computing principal strains
        _, sigma, vT = la.svd(cB_sup @ la.inv(cA_sup))
        rmsslist = [rmss(sigma)]
        # computing rmsd and optimized final structure
        cB_sup_final = vT.T @ np.diag(sigma) @ vT @ cA_sup
        c_sup_half = vT.T @ np.diag(sigma ** 0.5) @ vT @ cA_sup
        pA_sup = crystA_sup[2].T
        pB_sup = crystB_sup[2].T
        t0 = basinhopping(lambda z: minimize_rmsd_tp(c_sup_half, pA_sup, pB_sup + z.reshape(3,1))[0],
            [0.5, 0.5, 0.5], T=0.05, niter=75, minimizer_kwargs={"method": "BFGS"}).x
        _, ks = minimize_rmsd_tp(c_sup_half, pA_sup, pB_sup + t0.reshape(3,1))
        t0 = (pA_sup - pB_sup - ks).mean(axis=1, keepdims=True)
        pB_sup_final = pB_sup + ks + t0
        rmsdlist = [la.norm(c_sup_half @ (pA_sup - pB_sup_final)) / np.sqrt(pA_sup.shape[1])]
        table = np.array([[0], [0], mulist, rmsslist, rmsdlist], dtype=float).T
        print(f"\tmultiplicity ∈ {{{','.join(str(mu) for mu in range(1,mulist[0]+1) if mulist[0] % mu == 0)}}}\n\tprincipal strains: "
                + ", ".join([f"{100*(s-1):.2f}%" for s in sigma]) + f"\n\trmss = {rmsslist[0]:.4f}\n\trmsd = {rmsdlist[0]:.4f} Å.")
        print(f"\nTo produce a list of CSMs that may contain this CSM, you can run:\n"
                + f"\n\t$ crystmatch -E {mulist[0]} {rmsslist[0]+0.005:.3f} -I '{args.initial}' -F '{args.final}'")

    print("")

    # orientation analysis
    
    if args.orientation != None:
        vix, viy, viz, vfx, vfy, vfz, wix, wiy, wiz, wfx, wfy, wfz = args.orientation
        print(f"Comparing each CSM with the given OR:\n\t[{vix},{viy},{viz}] || [{vfx},{vfy},{vfz}]"
                + f"\n\t[{wix},{wiy},{wiz}] || [{wfx},{wfy},{wfz}]")
        r = orient_matrix([vix,viy,viz], [vfx,vfy,vfz], [wix,wiy,wiz], [wfx,wfy,wfz])
        anglelist = deviation_angle(crystA, crystB, slmlist, r, uspfix=args.uspfix)
        table = np.hstack((table, anglelist.reshape(-1,1)))
        if args.single is not None: print(f"The deviation angle is {anglelist[0]:.6f} = {anglelist[0]*180/np.pi:.4f}°\n")

    # saving NPZ and POSCAR files
    
    if args.enumeration is not None:
        # saving enumeration results in NPZ format
        metadata = np.array([job, crystA, crystB, mu_max, rmss_max, args.tolerance], dtype=object)
        np.savez(unique_filename("Saving CSMs to", f"CSM_LIST-{job}.npz"), *csm_bins, metadata=metadata, slmlist=slmlist, table=table)
    
    if args.single is not None:
        # saving primitive cells
        dirname = unique_filename(f"Saving optimized final structure (rotation-free, translated to minimize RMSD) to", "CSM_single")
        makedirs(dirname)
        save_poscar(f"{dirname}{sep}{args.initial.split(sep)[-1]}", crystA_sup)
        save_poscar(f"{dirname}{sep}{splitext(args.final)[0].split(sep)[-1]}-optimized{splitext(args.final)[1]}", (cB_sup_final.T, crystB_sup[1], pB_sup_final.T))
    
    if args.export is not None:
        if args.read is None:
            print("\nValueWarning: '--export' option is only available in 'read' mode.")
        else:
            print(f"\nExporting CSMs with indices {', '.join(map(str, args.export))} in {args.read} ...")
            for i in args.export:
                if i >= table.shape[0]:
                    print(f"\nIndexWarning: Index {i} is out of range (there are only {table.shape[0]} CSMs).")
                    continue
                slm = data['slmlist'][table[i,1].round().astype(int)]
                mu = table[i,2].round().astype(int)
                ind = np.sum(table[:i,2].round().astype(int) == mu)
                p = data[f'arr_{mu:d}'][ind,:,0]
                ks = data[f'arr_{mu:d}'][ind,:,1:].T
                crystA_sup, crystB_sup = int_arrays_to_pair(crystA, crystB, slm, p, ks)
                dirname = unique_filename(f"Saving the CSM with index {i:d} to", f"CSM_{i:d}")
                makedirs(dirname)
                save_poscar(f"{dirname}{sep}POSCAR_I", crystA_sup, crystname=f"CSM_{i:d} initial")
                save_poscar(f"{dirname}{sep}POSCAR_F", crystB_sup, crystname=f"CSM_{i:d} final")
    
    # saving CSV file and plotting

    if args.csv:
        if args.orientation is None:
            np.savetxt(unique_filename("Saving results to", f"TABLE-{job}.csv"), table * np.array([1,1,1,100,1]), fmt='%8d,%8d,%8d,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å')
        else:
            np.savetxt(unique_filename("Saving results to", f"TABLE-OR-{job}.csv"), table * np.array([1,1,1,100,1,180/np.pi]), fmt='%8d,%8d,%8d,%8.4f,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å, theta/°')
    
    if args.plot:
        save_scatter(unique_filename("Creating scatter plot in", f"PLOT-{job}.pdf"), table[:,3], table[:,4], table[:,2].round().astype(int), cbarlabel="Multiplicity")
        if args.orientation is not None: save_scatter(unique_filename("Creating scatter plot in", f"OR-{job}.pdf"),
                                                        table[:,3], table[:,4], table[:,5], cbarlabel="Deviation / °", cmap=plt.cm.get_cmap('jet'))
    
    print(f"\nTotal time spent: {time()-time0:.2f} seconds")

if __name__ == "__main__":
    main()