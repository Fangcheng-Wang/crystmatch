from .utils import *
from .enumeration import *
from .analysis import *
from os import makedirs
from os.path import sep
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
    t0 = time()
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
        job = f"{args.initial}-{args.final}-m{mu_max:d}s{rmss_max:.2f}"
        
        # enumerating SLMs
        print(f"\nEnumerating incongruent SLMs for mu <= {mu_max} and rmss <= {rmss_max:.4f} ...")
        slmlist = []
        for mu in range(1,mu_max+1): slmlist = slmlist + enumerate_slm(crystA, crystB, mu, rmss_max)
        print(f"\nA total of {len(slmlist):d} incongruent SLMs are enumerated.")
        slmlist = np.array(slmlist)
        _, ind = np.unique((slmlist[:,1,:,:] @ slmlist[:,2,:,:] @ la.inv(slmlist[:,0,:,:])).round(decimals=4), axis=0, return_index=True)
        slmlist = slmlist[ind]
        mulist = multiplicity(crystA, crystB, slmlist)
        rmsslist = rmss(sing_val(crystA, crystB, slmlist))
        print(f"They have {len(slmlist):d} distinct deformation gradients, and thus {len(slmlist):d} representative CSMs.")

        # computing representative CSMs
        mu_bins = []
        rmsdlist = []
        print(f"\nOptimizing RMSDs for {slmlist.shape[0]} shuffles ...")
        t = time()
        for mu in range(1,mu_max+1):
            n_csm = np.sum(mulist == mu)
            if n_csm == 0:
                mu_bins.append(None)
                continue
            shufflelist = np.zeros((n_csm, mu * zlcm, 4), dtype=int)
            for i in range(n_csm):
                d, p, ks = minimize_rmsd_translation(crystA, crystB, slmlist[mulist == mu][i])
                rmsdlist.append(d)
                shufflelist[i,:,0] = p
                shufflelist[i,:,1:] = ks.T
            mu_bins.append(shufflelist)
        rmsdlist = np.array(rmsdlist)
        print(f"\tObtained representative CSMs and their RMSDs (in {time()-t:.2f} seconds).")
        
        # saving results in NPZ format
        print(f"\nSaving results to CSM_LIST({job}).npz ...")
        metadata = np.array([job, crystA, crystB, mu_max, rmss_max, args.symprec], dtype=object)
        ind = np.lexsort((rmsdlist.round(decimals=4), rmsslist.round(decimals=4), mulist))
        slmlist = slmlist[ind]
        mulist = mulist[ind]
        rmsslist = rmsslist[ind]
        rmsdlist = rmsdlist[ind]
        table = np.array([np.arange(len(mulist)), np.arange(len(mulist)), mulist, rmsslist, rmsdlist], dtype=float).T
        np.savez(f".{sep}CSM_LIST({job}).npz", *mu_bins, metadata=metadata, slmlist=slmlist, table=table)

    elif args.analysis == -1:
        
        # mode: analysis (single CSM given by POSCARs)
        print("\nMode: Analysis (single CSM given by POSCARs).")
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        crystA = load_poscar(args.initial, symprec=args.symprec, to_primitive=False)
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystB = load_poscar(args.final, symprec=args.symprec, to_primitive=False)
        check_chem_comp(crystA[1], crystB[1])
        job = f"{args.initial}-{args.final}"
        if args.orientation != None: job += "-OR" + "".join(map(str, args.orientation))
        
        # analyzing the CSM
        pass
    
    else:

        # mode: analysis (loading CSM_LIST from NPZ file)
        print("\nMode: Analysis (loading CSM_LIST from {args.analysis}).")
        data = np.load(args.analysis)
        job, crystA, crystB, mu_max, rmss_max, symprec = data['metadata']
        print(f"Enumeration job: {job}\n\tMAX_MU: {mu_max}\n\tMAX_RMSS: {rmss_max}\n\tSYMPREC: {symprec}")

        # analyzing the CSMs
        pass
    
    # save and plot results
    if not args.nocsv:
        if args.orientation == None:
            np.savetxt(f".{sep}TABLE({job}).csv", table * np.array([1,1,1,100,1]), fmt='%8d,%8d,%8d,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å')
        else:
            np.savetxt(f".{sep}TABLE({job}).csv", table * np.array([1,1,1,100,1,180/np.pi]), fmt='%8d,%8d,%8d,%8.4f,%8.4f,%8.6f',
                        header=' index,  slm_id,      mu,  rmss/%,  rmsd/Å, theta/°')
    
    if not args.noplot:
        pass
        if args.orientation != None:
            pass
    
    print(f"\nTotal time spent: {time()-t0:.2f} seconds")

if __name__ == "__main__":
    main()