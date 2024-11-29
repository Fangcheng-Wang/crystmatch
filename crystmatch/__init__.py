from .io import *
from .core import *
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
__epilog__ = 'The current version (' + __version__ + ') may contain bugs. To get the latest version, please see: \
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
    parser.add_argument("-I", "--initial", type=str, metavar='POSCAR_I',
                        help="POSCAR file of the initial crystal structure")
    parser.add_argument("-F", "--final", type=str, metavar='POSCAR_F',
                        help="POSCAR file of the final crystal structure")
    parser.add_argument("-E", "--enumerate", nargs=2, type=float, metavar=('MAX_MU', 'MAX_RMSS'),
                        help="enumeration mode, using MAX_MU and MAX_RMSS as upper bounds")
    parser.add_argument("-A", "--analyze", nargs='?', const=-1, type=str, metavar='CSM_LIST',
                        help="analysis mode, using the CSM given by POSCAR_I and POSCAR_F (or from CSM_LIST if provided)")
    parser.add_argument("-s", "--symprec", type=float, default=1e-5, metavar='TOL',
                        help="tolerance for determining crystal symmetry, default is 1e-5")
    parser.add_argument("-c", "--csv", type=str, metavar='CSV_FILE',
                        help="the CSV file to store the analysis results")
    parser.add_argument("-e", "--export", nargs='+', type=int, metavar=('index1', 'index2'),
                        help="export CSMs from CSM_LIST with the given indices")
    parser.add_argument("-o", "--orientation", nargs=12, type=float, metavar=('vix','viy','viz','vfx','vfy','vfz','wix','wiy','wiz','wfx','wfy','wfz'),
                        help="benchmark CSMs by the orientation relationship (OR) vi||vf, wi||wf")
    parser.add_argument("--fixusp", action='store_true', help="use USP-fixed manner [1] when determining OR")
    parser.add_argument("--nocsv", action='store_true', help="not creating CSV file")
    parser.add_argument("--noplot", action='store_true', help="not creating PDF plots")
    parser.add_argument("-v", "--version", action='version', version="%(prog)s " + __version__)
    parser.add_argument("-B", nargs='*', type=str, metavar='', help="deprecated")
    args = parser.parse_args()

    if args.B != None:
        print("\nError: '-A -B' is deprecated since version 1.1.0, use '-I -F' instead. See 'crystmatch -h' for usage.")
        return
        
    if args.enumerate != None and args.analyze != None:
        # mode: conflict
        print("\nError: Cannot choose both enumeration and analysis modes")
        return

    if args.enumerate == None and args.analyze == None:
        # mode: unspecified
        mode = input("\nPlease enter a mode ('enumeration' or 'analysis'): ")
        if mode == 'enumeration':
            args.enumerate = input("\nEnter the upper bounds of multiplicity and RMSS (separated by space, such as '2 0.4'): ")
        elif mode == 'analysis':
            args.analyze = input("\nEnter the path of the CSM_LIST file or press Enter to use POSCARs: ")
            if args.analyze == '': args.analyze = -1
        else:
            print("\nError: Invalid mode")
            return
    
    # main computation (enumeration or analysis)

    if args.enumerate != None:
        
        # mode: enumeration
        print("\nMode: Enumeration")
        mu_max, rmss_max = args.enumerate
        mu_max = np.rint(mu_max).astype(int)
        if args.initial == None: args.initial = input("Enter the path of the initial POSCAR file: ")
        if args.final == None: args.final = input("Enter the path of the final POSCAR file: ")
        crystA = load_poscar(args.initial, symprec=args.symprec)
        crystB = load_poscar(args.final, symprec=args.symprec)
        check_chem_comp(crystA[1], crystB[1])
        job = f"{args.initial}-{args.final}-m{mu_max:d}s{rmss_max:.2f}"
        
        # start enumeration
        slist = complete_slm_list(crystA, crystB, mu_max, rmss_max)
        mulist = multiplicity(crystA, crystB, slist)
        rmsslist = root_mean_square_strain(sing_val(crystA, crystB, slist)).round(decimals=4)
        t = time()
        rmsdlist = []
        makedirs(f".{sep}CSM_LIST({job})")  # in progress!
        print(f"\nOptimizing atomic correspondences for {slist.shape[0]} SLMs ...")
        for i in slist.shape[0]:
            d, crystA_sup, crystB_sup, _ = minimize_rmsd(crystA, crystB, slist[i])
            if args.export != None:
                save_poscar(f".{sep}CSM_{i:d}({job}){sep}POSCAR_I", crystA_sup, crystname=args.initial)
                save_poscar(f".{sep}CSM_{i:d}({job}){sep}POSCAR_F", crystB_sup, crystname=args.final)
            rmsdlist.append(d)
        rmsdlist = np.array(rmsdlist).round(decimals=4)
        print(f"{slist.shape[0]} CSMs (with lowest RMSDs under own SLMs) generated (in {time()-t:.2f} seconds)")
        ind = np.lexsort((rmsdlist, rmsslist, mulist))
        slist = slist[ind]
        mulist = mulist[ind]
        rmsslist = rmsslist[ind]
        rmsdlist = rmsdlist[ind]
        np.save(f'.' + sep + 'SLM_LIST({mode},{label}).npy', slist)
        ind = np.arange(slist.shape[0])

    elif args.analyze == -1:
        # mode: analysis (from POSCARs)
        mode = 'analysis'
    
    else:
        # mode: analysis (from SLM_LIST)
        mode = 'analysis'
        print("\nMode: Analysis (load SLMs from file)")
        slist = np.load(args.analyze)
        ind = np.arange(slist.shape[0])
        if args.export != None: ind = np.array(args.export, dtype=int)
        slist = slist[ind]
        print(f"\n{slist.shape[0]} SLMs loaded from '{args.analyze}'")
        if args.csv != None:
            mulist, rmsslist, rmsdlist = np.loadtxt(args.csv, delimiter=',', usecols=(0,1,2,3), unpack=True)
            print(f"\nScores loaded from '{args.csv}'")
        if args.csv == None or args.export != None:
            if args.export != None: makedirs(f'.' + sep + 'CSM({mode},{label})')
            mulist = multiplicity(crystA, crystB, slist)
            rmsslist = root_mean_square_strain(sing_val(crystA, crystB, slist))
            t = time()
            rmsdlist = []
            print(f"\nOptimizing atomic correspondences for {slist.shape[0]} SLMs ...")
            for i in range(slist.shape[0]):
                d, crystA_sup, crystB_sup, _ = minimize_rmsd(crystA, crystB, slist[i])
                if args.export != None:
                    save_poscar(f'.' + sep + 'CSM({mode},{label})' + sep + 'CSM_{ind[i]}_A', crystA_sup, crystname=args.initial)
                    save_poscar(f'.' + sep + 'CSM({mode},{label})' + sep + 'CSM_{ind[i]}_B', crystB_sup, crystname=args.final)
                rmsdlist.append(d)
            rmsdlist = np.array(rmsdlist)
            print(f"{slist.shape[0]} CSMs (with lowest RMSDs under own SLMs) generated (in {time()-t:.2f} seconds)")
    
    # save and plot results
    
    if not args.noplot:                         # plot RMSS-RMSD-multiplicity
        scatter_colored(rmsslist, rmsdlist, multiplicity(crystA, crystB, slist), cbarlabel='Multiplicity $\mu$', \
            save=f'.' + sep + 'SCORE({mode},{label}).pdf')
        print(f"\nMultiplicities, RMSSs, and RMSDs are saved as 'SCORE({mode},{label}).pdf'")
    if args.orientation != None:                         # benchmark CSMs by orientation relationship
        b = np.array(args.orientation)
        r = orientation_relation(b[:3], b[3:6], b[6:9], b[9:])
        if args.fixusp: or_type = 'USP-fixed'
        else: or_type = 'rotation-free'
        print(f"\nBenchmarking SLMs by the OR ({or_type}): {b[:3]}_A || {b[3:6]}_B,  {b[6:9]}_A || {b[9:]}_B ...")
        anglelist = compare_orientation(crystA, crystB, slist, r, fix_usp=args.fixusp)
        if not args.noplot:                     # plot RMSS-RMSD-deviation
            scatter_colored(rmsslist, rmsdlist, anglelist, cbarlabel=f'Deviation (radian, {or_type})', \
                save=f'.' + sep + 'DEVIATION({or_type},{mode},{label}).pdf')
            print(f"\nDeviations are saved as 'DEVIATION({or_type},{mode},{label}).pdf'")
    if not args.nocsv:                          # save CSV file
        if args.orientation == None:
            np.savetxt(f'.' + sep + 'SCORE({mode},{label}).csv', np.column_stack((ind, mulist, rmsslist, rmsdlist)), \
                fmt=('%d','%d','%.4f','%.4f'), delimiter=',', header="index,mu,kappa,rmsd")
            print(f"\nIndices, multiplicities, RMSSs, and RMSDs are saved as 'SCORE({mode},{slist.shape[0]}).csv'")
        else:
            np.savetxt(f'.' + sep + 'SCORE({mode},{label}).csv', np.column_stack((ind, mulist, rmsslist, rmsdlist, anglelist)), \
                fmt=('%d','%d','%.4f','%.4f','%.4f'), delimiter=',', header=f"index,mu,kappa,rmsd,deviation({or_type})")
            print(f"\nIndices, multiplicities, RMSSs, RMSDs, and deviations({or_type}) appended to 'SCORE({mode},{slist.shape[0]}).csv'")
    print(f"\nTotal time spent: {time()-t0:.2f} seconds")

if __name__ == "__main__":
    main()