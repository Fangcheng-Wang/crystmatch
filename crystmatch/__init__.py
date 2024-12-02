from .utilities import *
from .enumeration import *
from .analysis import *
from os import makedirs
from os.path import sep, exists
import argparse

__name__ = "crystmatch"
__version__ = "1.1.0"
__author__ = "Fang-Cheng Wang"
__email__ = "wfc@pku.edu.cn"
__description__ = 'Enumerating and analyzing crystal-structure matches for solid-solid phase transitions.'
__url__ = 'https://github.com/Fangcheng-Wang/crystmatch'
__epilog__ = 'The current version (v' + __version__ + ') may contain bugs. To get the latest version, please see: \
\n\n\t' + __url__ + '\n\nIf you use crystmatch in your research, please cite the following paper: \n\n\t\
[1] F.-C. Wang et al., Physical Review Letters 132, 086101 (2024) (https://arxiv.org/abs/2305.05278)\n\n\
You are also welcome to contact me at ' + __email__ + ' for any questions or comments.'

def main():
    time0 = time()
    parser = argparse.ArgumentParser(
        prog=__name__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='%(prog)s [-h] [-initial POSCAR_I] [-final POSCAR_F] [mode options] [output options]',
        description=__description__, epilog=__epilog__)
    parser.add_argument("-v", "--version", action='version', version="%(prog)s " + __version__)
    parser.add_argument("-I", "--initial", type=str, metavar='POSCAR_I',
                        help="POSCAR file of the initial crystal structure")
    parser.add_argument("-F", "--final", type=str, metavar='POSCAR_F',
                        help="POSCAR file of the final crystal structure")
    parser.add_argument("-E", "--enumeration", nargs=2, type=float, metavar=('MAX_MU','MAX_RMSS'),
                        help="enumeration mode, using MAX_MU and MAX_RMSS as upper bounds")
    parser.add_argument("-A", "--analysis", nargs='?', const=-1, type=str, metavar='CSM_LIST',
                        help="analysis mode, using the CSM given by POSCAR_I and POSCAR_F (or from CSM_LIST if provided)")
    parser.add_argument("-e", "--export", nargs='+', type=int, metavar=('index1', 'index2'),
                        help="export CSMs from CSM_LIST with the given indices")
    parser.add_argument("-s", "--symprec", type=float, default=1e-5, metavar='TOL',
                        help="tolerance for determining crystal symmetry, default is 1e-5")
    parser.add_argument("--orientation", nargs=12, type=float, metavar=('vix','viy','viz','vfx','vfy','vfz','wix','wiy','wiz','wfx','wfy','wfz'),
                        help="benchmark CSMs by the orientation relationship (OR) vi||vf, wi||wf")
    parser.add_argument("--fixusp", action='store_true', help="use USP-fixed manner [1] when determining OR")
    parser.add_argument("--nocsv", action='store_true', help="not creating CSV file")
    parser.add_argument("--noplot", action='store_true', help="not creating PDF plots")
    parser.add_argument("-B", nargs='*', type=str, metavar='', help="deprecated")
    args = parser.parse_args()

    if args.B != None:
        print("\nError: '-A -B' is deprecated since version 1.1.0, use '-I -F' instead. See 'crystmatch -h' for usage.")
        return
        
    if args.enumeration != None and args.analysis != None:
        # mode: conflict
        print("\nError: Cannot choose both enumeration and analysis modes.")
        return

    if args.enumeration == None and args.analysis == None:
        # mode: unspecified
        mode = input("\nPlease enter a mode ('enumeration' or 'analysis'): ")
        if mode == 'enumeration':
            mu_max = int(input("Enter the maximum multiplicity: "))
            rmss_max = float(input("Enter the maximum root-mean-square-strain: "))
            args.enumeration = (mu_max, rmss_max)
        elif mode == 'analysis':
            args.analysis = input("\nEnter the path of the CSM_LIST file or press Enter to use POSCARs: ")
            if args.analysis == '': args.analysis = -1
        else:
            print("\nError: Invalid mode.")
            return
    
    if (args.analysis == None or args.analysis == -1) and args.export != None:
        print("\nError: Cannot export CSMs without specifying the CSM_LIST file.")
        args.export = None
    
    if args.analysis == None and args.orientation != None:
        print("\nWarning: Cannot perform orientation analysis in enumeration mode.")
        args.orientation = None
    
    # main computation (enumeration or analysis)

    if args.enumeration != None:
        
        # mode: enumeration
        print("\nMode: Enumeration.")
        mu_max, rmss_max = args.enumeration
        mu_max = np.rint(mu_max).astype(int)
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        crystA = load_poscar(args.initial, symprec=args.symprec)
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystB = load_poscar(args.final, symprec=args.symprec)
        check_chem_comp(crystA[1], crystB[1])
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        job = f"{args.initial.split('.')[0]}-{args.final.split('.')[0]}-m{mu_max:d}s{rmss_max:.2f}"
        
        # enumerating SLMs
        print(f"\nEnumerating incongruent SLMs for mu <= {mu_max} and rmss <= {rmss_max:.4f} ...")
        slmlist = []
        for mu in range(1,mu_max+1): slmlist = slmlist + enumerate_slm(crystA, crystB, mu, rmss_max)
        print(f"\nA total of {len(slmlist):d} incongruent SLMs are enumerated.")
        slmlist = np.array(slmlist)
        _, ind = np.unique((slmlist[:,1,:,:] @ slmlist[:,2,:,:] @ la.inv(slmlist[:,0,:,:])).round(decimals=4), axis=0, return_index=True)
        slmlist = slmlist[ind]
        mulist = multiplicity(crystA, crystB, slmlist)
        rmsslist = rms_minus1(sing_val(crystA, crystB, slmlist))
        ind = np.lexsort((rmsslist.round(decimals=4), mulist))
        slmlist = slmlist[ind]
        mulist = mulist[ind]
        rmsslist = rmsslist[ind]
        print(f"They have {len(slmlist):d} distinct deformation gradients, and thus {len(slmlist):d} representative CSMs.")

        # computing representative CSMs
        mu_bins = [np.array("arr_mu.npy saves the shuffles of those CSMs with multiplicity mu", dtype=str)]
        rmsdlist = []
        print(f"\nOptimizing RMSDs for {slmlist.shape[0]} shuffles ...")
        time1 = time()
        for mu in range(1,mu_max+1):
            n_csm = np.sum(mulist == mu)
            if n_csm == 0:
                mu_bins.append(np.array(f"No CSM with multiplicity {mu}", dtype=str))
                continue
            shufflelist = np.zeros((n_csm, mu * zlcm, 4), dtype=int)
            for i in range(n_csm):
                d, p, ks, _ = minimize_rmsd(crystA, crystB, slmlist[mulist == mu][i])
                rmsdlist.append(d)
                shufflelist[i,:,0] = p
                shufflelist[i,:,1:] = ks.T
            mu_bins.append(shufflelist)
        print(f"\tObtained representative CSMs and their RMSDs (in {time()-time1:.2f} seconds).")
        table = np.array([np.arange(len(mulist)), np.arange(len(mulist)), mulist, rmsslist, rmsdlist], dtype=float).T

    elif args.analysis == -1:
        
        # mode: analysis (single CSM given by POSCARs)
        print("\nMode: Analysis (single CSM).")
        print("\nWarning: The CSM is uniquely determined by initial and final POSCAR files. I hope you know what you are doing.")
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        crystA_sup = load_poscar(args.initial, symprec=args.symprec, to_primitive=False)
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystB_sup = load_poscar(args.final, symprec=args.symprec, to_primitive=False)
        if not crystA_sup[1].shape[0] == crystB_sup[1].shape[0]:
            print("\nError: Numbers of atoms in two crysts are not equal!")
            return
        if not (crystA_sup[1] == crystB_sup[1]).all():
            print("\nError: Atom species in two crysts do not match!")
            return
        job = f"{args.initial.split('.')[0]}-{args.final.split('.')[0]}-single"
        if args.orientation != None: job += "-OR" + "".join(map(str, args.orientation))
        
        # computing mu
        _, numbers = np.unique(crystA_sup[1], return_inverse=True)
        _, _, numbersA = find_primitive((crystA_sup[0], crystA_sup[2], numbers))
        _, _, numbersB = find_primitive((crystB_sup[0], crystB_sup[2], numbers))
        mulist = [len(numbers) // np.lcm(len(numbersA), len(numbersB))]
        # computing strains
        cA_sup = crystA_sup[0].T
        cB_sup = crystB_sup[0].T
        _, sigma, vT = la.svd(cB_sup @ la.inv(cA_sup))
        rmsslist = [rms_minus1(sigma)]
        # computing rmsd
        cB_sup_final = vT.T @ np.diag(sigma) @ vT @ cA_sup
        c_sup_half = vT.T @ np.diag(sigma ** 0.5) @ vT @ cA_sup
        pA_sup = crystA_sup[2].T
        pB_sup = crystB_sup[2].T
        t0 = brute(lambda z: minimize_rmsd_tp(c_sup_half, pA_sup, pB_sup + z.reshape(3,1))[0],
            ((0,1),(0,1),(0,1)), Ns=10, finish=None)
        _, ks = minimize_rmsd_tp(c_sup_half, pA_sup, pB_sup + t0.reshape(3,1))
        t0 = (pA_sup - pB_sup - ks).mean(axis=1)
        pB_sup_final = pB_sup + ks + t0
        rmsdlist = [la.norm(c_sup_half @ (pA_sup - pB_sup_final)) / np.sqrt(pA_sup.shape[1])]
        assert args.orientation == None, "Error: Sorry, this function is still in development."
        print(f"\nThe CSM has properties:\n\tmultiplicity = {mulist[0]}\n\tprincipal strains = " + 
                ", ".join([f"{100*(s-1):.2f}%" for s in sigma]) + f"\n\trmss = {rmsslist[0]:.4f}\n\trmsd = {rmsdlist[0]:.4f} Å.")
        table = np.array([[0], [0], mulist, rmsslist, rmsdlist], dtype=float).T
    
    else:

        # mode: analysis (loading CSM_LIST from NPZ file)
        print("\nMode: Analysis (loading CSM_LIST from {args.analysis}).")
        data = np.load(args.analysis)
        job, crystA, crystB, mu_max, rmss_max, symprec = data['metadata']
        print(f"Loaded: {job}\n\tMAX_MU: {mu_max}\n\tMAX_RMSS: {rmss_max}\n\tSYMPREC: {symprec}")
        if args.orientation != None: job += "-OR" + "".join(map(str, args.orientation))
        
        # analyzing the CSMs
        pass
    
    # save and plot results
    
    outdir = job
    if exists(f"{job}"):
        r = 1
        while exists(f"{job}_{r}"): r += 1
        outdir = f"{job}_{r}"
    makedirs(outdir)
    
    if args.enumeration != None:
        # saving enumeration results in NPZ format
        print(f"\nSaving results to {outdir}{sep}CSM_LIST({job}).npz ...")
        metadata = np.array([job, crystA, crystB, mu_max, rmss_max, args.symprec], dtype=object)
        np.savez(f"{outdir}{sep}CSM_LIST({job}).npz", *mu_bins, metadata=metadata, slmlist=slmlist, table=table)
    elif args.analysis == -1:
        # saving final POSCAR (rotation-free and RMSD-minimized)
        print(f"\nSaving final crystal structure (with rotation-free orientation and RMSD-minimized positions) \
            to {outdir}{sep}POSCAR_F ...")
        save_poscar(f"{outdir}{sep}POSCAR_I", crystA_sup, crystname = args.initial.split('.')[0])
        save_poscar(f"{outdir}{sep}POSCAR_F", (cB_sup_final.T, crystB_sup[1], pB_sup_final.T),
                    crystname = args.final.split('.')[0] + "-final")
    
    if args.export != None:
        pass
    
    if not args.nocsv:
        print(f"\nSaving results to {outdir}{sep}TABLE({job}).csv ...")
        if args.orientation == None:
            np.savetxt(f"{outdir}{sep}TABLE({job}).csv", table * np.array([1,1,1,100,1]), fmt='%8d,%8d,%8d,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å')
        else:
            np.savetxt(f"{outdir}{sep}TABLE({job}).csv", table * np.array([1,1,1,100,1,180/np.pi]), fmt='%8d,%8d,%8d,%8.4f,%8.4f,%8.6f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å, theta/°')
    
    if not (args.noplot or args.analysis == -1):
        print(f"\nCreating plots in {outdir}{sep}PLOT({job}).pdf ...")
        scatter_colored(f"{outdir}{sep}rmsd-rmss-mu({job}).pdf", rmsslist, rmsdlist, mulist, cbarlabel="Multiplicity")
        if args.orientation != None:
            pass    # plot OR
    
    print(f"\nTotal time spent: {time()-time0:.2f} seconds")

if __name__ == "__main__":
    main()