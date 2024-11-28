"""
The command-line interface for crystmatch.
"""

import argparse
from .io import *
from .core import *
from .analysis import *
from os import makedirs

epilog = 'The current version of crystmatch may contain bugs; \
        please see https://github.com/Fangcheng-Wang/crystmatch for the latest version. \
        You are also welcome to report any issues or suggestions to wfc@pku.edu.cn. \
        If you use crystmatch in your research, please cite the following paper: \n\n\
        [1] F.-C. Wang et al., Physical Review Letters 132, 086101 (2024) (https://arxiv.org/abs/2305.05278).'

def main():
    t0 = time()
    parser = argparse.ArgumentParser(
        prog='crystmatch',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Enumerating and analyzing crystal-structure matches for solid-solid phase transitions.',
        epilog=epilog)
    parser.add_argument("-I", "--initial", type=str, metavar='POSCAR_I',
                        help="POSCAR file of the initial crystal structure")
    parser.add_argument("-F", "--final", type=str, metavar='POSCAR_F',
                        help="POSCAR file of the final crystal structure")
    parser.add_argument("-E", "--enumerate", nargs=2, metavar=('MAX_MU', 'MAX_RMSS'), type=float,
                        help="Enumeration mode with MAX_MU and MAX_RMSS")
    parser.add_argument("-A", "--analyze", nargs='?', metavar='SLM_LIST', const=-1, type=str,
                        help="Analysis mode, using the CSM given by the initial and final POSCARs, or SLMs loaded from SLM_LIST")
    parser.add_argument("-s", "--load_score", nargs='?', type=str, const='./SCORE.csv',
                        help="Load SCORE.csv")
    parser.add_argument("-i", "--index", nargs='+', type=int,
                        help="Specify CSM (SLM) indices")
    parser.add_argument("-p", "--poscar", action='store_true',
                        help="Create POSCAR files for each CSM")
    parser.add_argument("-x", "--xdatcar", nargs='?', type=int, const=50,
                        help="Create XDATCAR file from given CSM")
    parser.add_argument("-o", "--outdir", type=str, default='.',
                        help="Output directory")
    parser.add_argument("-b", "--bm", nargs=12, type=float,
                        help="Benchmark by an OR. See documentation")
    parser.add_argument("--fixusp", action='store_true',
                        help="Fix USP when determining OR")
    parser.add_argument("--nocsv", action='store_true',
                        help="Not creating SCORE.csv")
    parser.add_argument("--noplot", action='store_true',
                        help="Not plotting")
    parser.add_argument("-h", "--help", action='store_true',
                        help="Not plotting")
    parser.add_argument("-B", "--poscar_B", type=str)       # deprecated
    args = parser.parse_args()

    if args.poscar_B != None:
        print("Error: '-B' is deprecated, use '-F' instead. See 'crystmatch -h' for help.")
        return
    
    crystA = load_poscar(args.initial)
    crystB = load_poscar(args.final)
    check_chem_comp(crystA[1], crystB[1])
    makedirs(args.outdir, exist_ok=True)
    label = args.initial + '-to-' + args.final
    
    # Determine mode: Enumeration / Analysis (from POSCARs) / Analysis (from SLM_LIST)
    
    mode = None
    if args.enumerate != None:
        # mode: enumeration
        mode = 'enumerate'
        print("\nMode: Enumerate SLMs")
        mu_max, kappa_max = args.enumerate
        mu_max = np.rint(mu_max).astype(int)
        slist = complete_slm_list(crystA, crystB, mu_max, kappa_max)
        mulist = multiplicity(crystA, crystB, slist)
        kappalist = root_mean_square_strain(sing_val(crystA, crystB, slist)).round(decimals=4)
        t = time()
        rmsdlist = []
        if args.poscar: makedirs(f'{args.outdir}/CSM_POSCAR')
        print(f"\nOptimizing atomic correspondences for {slist.shape[0]} SLMs ...")
        for s in slist:
            d, crystA_sup, crystB_sup, _ = minimize_rmsd(crystA, crystB, s)
            if args.poscar:
                save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{len(rmsdlist)}_A', crystA_sup, crystname=args.initial)
                save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{len(rmsdlist)}_B', crystB_sup, crystname=args.final)
            rmsdlist.append(d)
        rmsdlist = np.array(rmsdlist).round(decimals=4)
        print(f"{slist.shape[0]} CSMs (with lowest RMSDs under own SLMs) generated (in {time()-t:.2f} seconds)")
        ind = np.lexsort((rmsdlist, kappalist, mulist))
        slist = slist[ind]
        mulist = mulist[ind]
        kappalist = kappalist[ind]
        rmsdlist = rmsdlist[ind]
        np.save(f'{args.outdir}/SLM_LIST({mode},{label}).npy', slist)
        ind = np.arange(slist.shape[0])
    elif args.analyze == -1:
        # mode: analysis (from POSCARs)
        mode = 'analyze'
    elif args.analyze != None:
        # mode: analysis (from SLM_LIST)
        mode = 'analyze'
        print("\nMode: Load SLMs from file")
        slist = np.load(args.analyze)
        ind = np.arange(slist.shape[0])
        if args.index != None: ind = np.array(args.index, dtype=int)
        slist = slist[ind]
        print(f"\n{slist.shape[0]} SLMs loaded from '{args.analyze}'")
        if args.load_score != None:
            mulist, kappalist, rmsdlist = np.loadtxt(args.load_score, delimiter=',', usecols=(0,1,2,3), unpack=True)
            print(f"\nScores loaded from '{args.load_score}'")
        if args.load_score == None or args.poscar or args.xdatcar:
            if args.poscar or args.xdatcar: makedirs(f'{args.outdir}/CSM({mode},{label})')
            mulist = multiplicity(crystA, crystB, slist)
            kappalist = root_mean_square_strain(sing_val(crystA, crystB, slist))
            t = time()
            rmsdlist = []
            print(f"\nOptimizing atomic correspondences for {slist.shape[0]} SLMs ...")
            for i in range(slist.shape[0]):
                d, crystA_sup, crystB_sup, _ = minimize_rmsd(crystA, crystB, slist[i])
                if args.poscar:
                    save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{ind[i]}_A', crystA_sup, crystname=args.initial)
                    save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{ind[i]}_B', crystB_sup, crystname=args.final)
                if args.xdatcar != None:
                    save_trajectory(f'{args.outdir}/CSM({mode},{label})/TRAJECTORY_{ind[i]}', crystA_sup, crystB_sup, num=args.xdatcar, \
                        crystname=f'{args.initial}-{args.final}')
                rmsdlist.append(d)
            rmsdlist = np.array(rmsdlist)
            print(f"{slist.shape[0]} CSMs (with lowest RMSDs under own SLMs) generated (in {time()-t:.2f} seconds)")
    else:
        print("Error: Must choose a mode using '-E' or '-L'")
        return
    
    # Output control
    
    if not args.noplot:                         # plot RMSS-RMSD-multiplicity
        scatter_colored(kappalist, rmsdlist, multiplicity(crystA, crystB, slist), cbarlabel='Multiplicity $\mu$', \
            save=f'{args.outdir}/SCORE({mode},{label}).pdf')
        print(f"\nMultiplicities, RMSSs, and RMSDs are saved as 'SCORE({mode},{label}).pdf'")
    if args.bm != None:                         # benchmark CSMs by orientation relationship
        b = np.array(args.bm)
        r = orientation_relation(b[:3], b[3:6], b[6:9], b[9:])
        if args.fixusp: or_type = 'USP-fixed'
        else: or_type = 'rotation-free'
        print(f"\nBenchmarking SLMs by the OR ({or_type}): {b[:3]}_A || {b[3:6]}_B,  {b[6:9]}_A || {b[9:]}_B ...")
        anglelist = compare_orientation(crystA, crystB, slist, r, fix_usp=args.fixusp)
        if not args.noplot:                     # plot RMSS-RMSD-deviation
            scatter_colored(kappalist, rmsdlist, anglelist, cbarlabel=f'Deviation (radian, {or_type})', \
                save=f'{args.outdir}/DEVIATION({or_type},{mode},{label}).pdf')
            print(f"\nDeviations are saved as 'DEVIATION({or_type},{mode},{label}).pdf'")
    if not args.nocsv:                          # save CSV file
        if args.bm == None:
            np.savetxt(f'{args.outdir}/SCORE({mode},{label}).csv', np.column_stack((ind, mulist, kappalist, rmsdlist)), \
                fmt=('%d','%d','%.4f','%.4f'), delimiter=',', header="index,mu,kappa,rmsd")
            print(f"\nIndices, multiplicities, RMSSs, and RMSDs are saved as 'SCORE({mode},{slist.shape[0]}).csv'")
        else:
            np.savetxt(f'{args.outdir}/SCORE({mode},{label}).csv', np.column_stack((ind, mulist, kappalist, rmsdlist, anglelist)), \
                fmt=('%d','%d','%.4f','%.4f','%.4f'), delimiter=',', header=f"index,mu,kappa,rmsd,deviation({or_type})")
            print(f"\nIndices, multiplicities, RMSSs, RMSDs, and deviations({or_type}) appended to 'SCORE({mode},{slist.shape[0]}).csv'")
    print(f"\nTotal time spent: {time()-t0:.2f} seconds")

if __name__ == "__main__":
    main()