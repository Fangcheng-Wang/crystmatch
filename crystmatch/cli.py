"""
The command-line interface for crystmatch.
"""

import argparse
from .io import *
from .core import *
from .analysis import *
from os import makedirs

def main():
    t0 = time()
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--cryst_a", type=str, default='./POSCAR_A', help="POSCAR file of crystal structure A")
    parser.add_argument("-B", "--cryst_b", type=str, default='./POSCAR_B', help="POSCAR file of crystal structure B")
    parser.add_argument("-E", "--enumerate", nargs=2, type=float, help="Enumeration: MU_MAX KAPPA_MAX")
    parser.add_argument("-L", "--load_slm", nargs='?', type=str, const='./SLM_LIST.npy', help="Load from file: SLM_LIST")
    parser.add_argument("-s", "--load_score", nargs='?', type=str, const='./SCORE.csv', help="Load SCORE.csv")
    parser.add_argument("-i", "--index", nargs='+', type=int, help="Specify CSM (SLM) indices")
    parser.add_argument("-p", "--poscar", action='store_true', help="Create POSCAR files for each CSM")
    parser.add_argument("-x", "--xdatcar", nargs='?', type=int, const=50, help="Create XDATCAR file from given CSM")
    parser.add_argument("-o", "--outdir", type=str, default='.', help="Output directory")
    parser.add_argument("-b", "--bm", nargs=12, type=float, help="Benchmark by an OR. See documentation")
    parser.add_argument("--fixusp", action='store_true', help="Fix USP when determining OR")
    parser.add_argument("--nocsv", action='store_true', help="Not creating SCORE.csv")
    parser.add_argument("--noplot", action='store_true', help="Not plotting")
    args = parser.parse_args()
    # enumerate or load SLMs, optimize atomic correspondences for CSMs
    crystA = load_poscar(args.cryst_a)
    crystB = load_poscar(args.cryst_b)
    check_chem_comp(crystA[1], crystB[1])
    makedirs(args.outdir, exist_ok=True)
    mode = 'none'
    label = args.cryst_a + ',' + args.cryst_b
    if args.enumerate != None:                  # mode: enumeration
        mode = 'enum'
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
                save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{len(rmsdlist)}_A', crystA_sup, crystname=args.cryst_a)
                save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{len(rmsdlist)}_B', crystB_sup, crystname=args.cryst_b)
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
    elif args.load_slm != None:                 # mode: load SLMs from file
        mode = 'load'
        print("\nMode: Load SLMs from file")
        slist = np.load(args.load_slm)
        ind = np.arange(slist.shape[0])
        if args.index != None: ind = np.array(args.index, dtype=int)
        slist = slist[ind]
        print(f"\n{slist.shape[0]} SLMs loaded from '{args.load_slm}'")
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
                    save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{ind[i]}_A', crystA_sup, crystname=args.cryst_a)
                    save_poscar(f'{args.outdir}/CSM({mode},{label})/CSM_{ind[i]}_B', crystB_sup, crystname=args.cryst_b)
                if args.xdatcar != None:
                    save_trajectory(f'{args.outdir}/CSM({mode},{label})/TRAJECTORY_{ind[i]}', crystA_sup, crystB_sup, num=args.xdatcar, \
                        crystname=f'{args.cryst_a}-{args.cryst_b}')
                rmsdlist.append(d)
            rmsdlist = np.array(rmsdlist)
            print(f"{slist.shape[0]} CSMs (with lowest RMSDs under own SLMs) generated (in {time()-t:.2f} seconds)")
    if mode == 'none':
        print("Error: Must choose a mode using '-E' or '-L'")
        return
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