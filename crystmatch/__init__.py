from .utilities import *
from .enumeration import *
from .analysis import *
from time import perf_counter
import sys
import argparse

__name__ = "crystmatch"
__version__ = "1.2.0"
__author__ = "Fang-Cheng Wang"
__email__ = "wfc@pku.edu.cn"
__description__ = 'Enumerating and analyzing crystal-structure matches for solid-solid phase transitions.'
__url__ = 'https://fangcheng-wang.github.io/crystmatch/'
__epilog__ = 'The current version is v' + __version__ + '. To get the latest version, please run:\
\n\n\t$ pip3 install --upgrade crystmatch\n\nWe also recommend you to see the documentation at:\
\n\n\t' + __url__ + '\n\nIf you use crystmatch in your research, please cite the following paper: \n\n\t\
[1] FC Wang, QJ Ye, YC Zhu, and XZ Li, Physical Review Letters 132, 086101 (2024) (https://arxiv.org/abs/2305.05278)\n\n\
You are also welcome to contact us at ' + __email__ + ' for any questions, feedbacks or comments.'

def main():
    time0 = perf_counter()
    sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)
    sys.stderr = open(sys.stdout.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)

    parser = argparse.ArgumentParser(
        prog = __name__,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        usage = "%(prog)s [-h] [-V] [-E POSCAR_I POSCAR_F MAX_MU MAX_STRAIN] [-R CSMLIST [IND1 IND2 ...]] [-D POSCAR_I POSCAR_F] [-l]"\
                + "[-e CSMCAR] [-a MAX_D] [-t TOL] [-i [SIZE]] [-o ASSUM] [-p [ASSUM]] [-x [ASSUM]] [-c] [-s]",
        description = __description__,
        epilog = __epilog__)
    
    parser.add_argument("-V", "-v", "--version", action='version', version="%(prog)s " + __version__)
    parser.add_argument("-E", "--enumerate", nargs=4, type=str, metavar=('POSCAR_I', 'POSCAR_F', 'MAX_MU', 'MAX_STRAIN'),
                        help="enumerate representative CSMs between POSCAR_I and POSCAR_F, with MAX_MU and MAX_STRAIN as upper bounds for the multiplicity and RMSS")
    parser.add_argument("-R", "--read", nargs='+', type=str, metavar=('CSMLIST', 'IND1 IND2'),
                        help="read CSMs from CSMLIST, with optional IND1 IND2 ... to specify the indices of the CSMs to be loaded")
    parser.add_argument("-D", "--direct", nargs=2, type=str, metavar=('POSCAR_I', 'POSCAR_F'),
                        help="directly read a single CSM from POSCAR_I and POSCAR_F")
    parser.add_argument("-l", "--literal", action='store_true',
                        help="when using --direct, do not add integers to fractional coordinates")
    parser.add_argument("-e", "--extra", type=str, metavar='CSMCAR',
                        help="load extra parameters from CSMCAR, including the elastic tensors, atomic weights, and orientation relationship used to benchmark CSMs")
    parser.add_argument("-a", "--all", type=float, metavar='MAX_D',
                        help="enumerate all CSMs for each SLM, with MAX_D as the upper bound for the shuffle distance")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-3, metavar='TOL',
                        help="tolerance in angstroms for detecting symmetry (default is 1e-3)")
    parser.add_argument("-i", "--interact", nargs='?', type=float, default=None, const=1.5, metavar='SIZE',
                        help="interactively visualize each CSM using a 3D plot, with SIZE controlling the radius of the cluster to display (default is 1.5)")
    parser.add_argument("-o", "--orientation", type=str, choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="calculate the deviation angle from the given orientation relationship for each CSM (must be used with --extra)")
    parser.add_argument("-p", "--poscar", nargs='?', type=str, default=None, const='norot', choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="export each enumerated/read CSM as a pair of POSCAR files (default is 'norot')")
    parser.add_argument("-x", "--xdatcar", nargs='?', type=str, default=None, const='norot', choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="export each enumerated/read CSM as an XDATCAR file with 100 frames (default is 'norot')")
    parser.add_argument("-c", "--csv", action='store_true',
                        help="create CSV file containing basic properties of all enumerated/read CSMs")
    parser.add_argument("-s", "--scatter", action='store_true',
                        help="create scatter plot in PDF format containing basic properties all enumerated/read CSMs")

    args = parser.parse_args()
    
    # check if exactly one mode is specified
    n_modes = sum(1 for x in [args.enumerate, args.read, args.direct] if x is not None)
    if n_modes >= 2: raise ValueError("Cannot choose more than one mode.")
    elif n_modes == 0: raise ValueError("No mode is specified.")
    mode = 'enumerate' if args.enumerate is not None else ('read' if args.read is not None else 'direct')
    
    # check if --literal is used with --direct
    if args.literal and mode != 'direct': raise ValueError("'--literal' can only be used with '--direct'.")
    
    # automatically enable --csv and --scatter
    if args.enumerate is not None or args.all is not None or args.orientation is not None:
        args.csv = True
        args.scatter = True
    
    # load global settings
    voigtA, voigtB, weight_func, ori_rel = None, None, None, None
    if args.extra is not None: voigtA, voigtB, weight_func, ori_rel = load_csmcar(args.extra)
    tol = args.tolerance
    
    # get CSMs
    
    print(f"\nUsing '{mode}' mode to get CSMs.")
    
    if mode == 'enumerate':

        fileA = args.enumerate[0]
        fileB = args.enumerate[1]
        max_mu = int(args.enumerate[2])
        max_strain = float(args.enumerate[3])
        jobname = f"m{args.enumerate[2].strip()}s{args.enumerate[3].strip()}"
        crystA = load_poscar(fileA, tol=tol)
        crystB = load_poscar(fileB, tol=tol)
        check_stoichiometry(crystA[1], crystB[1])
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))

        if voigtA is None and voigtB is None:
            strain = rmss
            print("Using the root-mean-squared strain (RMSS) to quantify deformation.")
        elif voigtA is not None and voigtB is not None:
            elasticA = voigt_to_tensor(voigtA, cryst=crystA, tol=tol)
            elasticB = voigt_to_tensor(voigtB, cryst=crystB, tol=tol)
            strain = strain_energy_func(elasticA, elasticB)
            print("Using the estimated strain energy (w) to quantify deformation.")
        else:
            raise ValueError("Both voigtA and voigtB must be provided if either is provided.")

        if max_mu < 1: raise ValueError("Multiplicity should be a positive integer.")
        elif max_strain <= 0: raise ValueError("Root-mean-squared strain should be a positive float.")
        if max_mu >= 8 or (strain==rmss and max_strain > 0.4):
            print(f"Warning: Current MAX_MU = {max_mu:d} and MAX_STRAIN = {max_strain:.2f} may result in a large number of SLMs, which may take "
                    + "a very long time to enumerate.")
        
        if args.all is None:
            slmlist, pct_arrs, mulist, strainlist, dlist = enumerate_rep_csm(crystA, crystB, max_mu, max_strain, strain=strain, weight_func=weight_func, tol=tol)
            slm_ind = np.arange(len(slmlist))
        else:
            max_d = args.all
            jobname += f"d{max_d:.2f}"
            if max_d > 2.0: print(f"Warning: Current MAX_D = {max_d:.2f} may result in a large number of CSMs, which may take a very long time to enumerate.")
            if zlcm * max_mu >= 12: print(f"Warning: Current MAX_MU = {max_mu:d} leads to a maxinum period of {zlcm * max_mu:d}, which may result in too many CSMs.")
            slmlist, slm_ind, pct_arrs, mulist, strainlist, dlist = enumerate_all_csm(crystA, crystB, max_mu, max_strain, max_d, strain=strain, weight_func=weight_func, tol=tol)
        
        if strain == rmss:
            header = "# csm_id,  slm_id,      mu,  period,    rmss,  rmsd/Å"
            data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, dlist]).T
        else:
            rmsslist = rmss(deformation_gradient(crystA, crystB, slmlist))[slm_ind]
            header = "# csm_id,  slm_id,      mu,  period,       w,    rmss,  rmsd/Å"
            data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, rmsslist, dlist]).T

    elif mode == 'read':
        
        jobname = 'read'
        filenpz = args.read[0]
        crystA, crystB, slmlist, slm_ind, pct_arrs, table = load_npz(filenpz, verbose=True)
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))
        max_mu = imt_multiplicity(crystA, crystB, slmlist[slm_ind]).max()
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        header = table.header
        data = table.data
        
        if len(args.read) > 1:
            indices = [int(i) for i in args.read[1:]]
            print(f"\tReading CSMs with indices:", end='')
            max_mu = imt_multiplicity(crystA, crystB, slmlist[slm_ind[indices]]).max()
            slmlist_temp = np.array([], dtype=ind).reshape(0,3,3,3)
            slm_ind_temp = []
            pct_arrs_temp = [NPZ_ARR_COMMENT] + [[] for _ in range(max_mu)]
            for i in indices:
                print(f" {i:4d}", end='')
                slm, p, ks = unzip_csm(i, crystA, crystB, slmlist, slm_ind, pct_arrs)
                where = np.nonzero((slm == slmlist).all(axis=(1,2,3)).any())[0]
                if where.shape[0] == 0:
                    slm_ind_temp.append(slmlist_temp.shape[0])
                    slmlist_temp = np.concatenate((slmlist_temp, slm.reshape(-1,3,3,3)), axis=0)
                else:
                    slm_ind_temp.append(where[0])
                mu = imt_multiplicity(crystA, crystB, slm)
                pct_arrs_temp[mu].append(zip_pct(p, ks))
            print('\n', end='')
            slmlist = slmlist_temp
            slm_ind = slm_ind_temp
            pct_arrs = pct_arrs_temp
            data = data[indices,:]
        
        if args.all is not None:
            
            if voigtA is None and voigtB is None:
                strain = rmss
            elif voigtA is not None and voigtB is not None:
                elasticA = voigt_to_tensor(voigtA, cryst=crystA, tol=tol)
                elasticB = voigt_to_tensor(voigtB, cryst=crystB, tol=tol)
                strain = strain_energy_func(elasticA, elasticB)
            else:
                raise ValueError("Both voigtA and voigtB must be provided if either is provided.")
                
            max_d = args.all
            jobname += f"-d{max_d:.2f}"
            if max_d > 2.0: print(f"Warning: Current MAX_D = {max_d:.2f} may result in a large number of CSMs, which may take a very long time to enumerate.")
            if zlcm * max_mu >= 12: print(f"Warning: Current MAX_MU = {max_mu:d} leads to a maxinum period of {zlcm * max_mu:d}, which may result in too many CSMs.")
            print(f"\nEnumerating all CSMs with MAX_D = {max_d:.2f} for {slmlist.shape[0]:d} SLMs ...")
            slm_ind = []
            pct_arrs = [NPZ_ARR_COMMENT] + [np.array([], dtype=int).reshape(0, mu * zlcm, 4) for mu in range(1, max_mu + 1)]
            dlist = np.array([], dtype=float)
            progress_bar = tqdm(total=slmlist.shape[0], position=0, desc=f"\r\tSLMs", ncols=60, mininterval=0.5,
                    bar_format=f'{{desc}}: {{n_fmt:>3s}}/{{total_fmt:>3s}} |{{bar}}| [elapsed {{elapsed:5s}}, remaining {{remaining:5s}}]')
            
            for i, slm in enumerate(slmlist):
                mu = imt_multiplicity(crystA, crystB, slm)
                pcts, ds = enumerate_pct(crystA, crystB, slm, max_d, weight_func=weight_func, verbose=0, warning_threshold=100000)
                slm_ind += [i] * pcts.shape[0]
                pct_arrs[mu] = np.concatenate([pct_arrs[mu], pcts], axis=0)
                dlist = np.concatenate([dlist, ds], axis=0)
                progress_bar.update(1)
            progress_bar.close()
            slm_ind = np.array(slm_ind, dtype=int)
            mulist = imt_multiplicity(crystA, crystB, slmlist)[slm_ind]
            strainlist = strain(deformation_gradient(crystA, crystB, slmlist))[slm_ind]

            if strain == rmss:
                header = "# csm_id,  slm_id,      mu,  period,    rmss,  rmsd/Å"
                data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, dlist]).T
            else:
                rmsslist = rmss(deformation_gradient(crystA, crystB, slmlist))[slm_ind]
                header = "# csm_id,  slm_id,      mu,  period,       w,    rmss,  rmsd/Å"
                data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, rmsslist, dlist]).T

    elif mode == 'direct':

        jobname = 'direct'
        fileA = args.direct[0]
        fileB = args.direct[1]
        crystA_sup = load_poscar(args.initial, tol=tol, to_primitive=False)
        crystB_sup = load_poscar(args.final, tol=tol, to_primitive=False)
        if not crystA_sup[1].shape[0] == crystB_sup[1].shape[0]: raise ValueError("The numbers of atoms in two POSCAR files are not equal!")
        check_stoichiometry(crystA_sup[1], crystB_sup[1])
        crystA, crystB, slm0, p0, ks0 = cryst_to_csm(crystA_sup, crystB_sup, tol=tol)
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))

        if not args.literal:
            d0 = csm_distance(crystA, crystB, slm0, p0, ks0, weight_func=weight_func)
            ks0 = optimize_ct(crystA, crystB, slm0, p0, weight_func=weight_func)[1]
        slm, p, ks = primitive_shuffle(crystA, crystB, slm0, p0, ks0)
        d = csm_distance(crystA, crystB, slm, p, ks, weight_func=weight_func)
        if not args.literal and d < d0 - tol:
            print("\nWarning: Integers are added to fractional coordinates to minimize the shuffle distance. Use '--literal' to disable this behavior.")
        
        if voigtA is None and voigtB is None:
            strain = rmss
        elif voigtA is not None and voigtB is not None:
            elasticA = voigt_to_tensor(voigtA, cryst=crystA, tol=tol)
            elasticB = voigt_to_tensor(voigtB, cryst=crystB, tol=tol)
            strain = strain_energy_func(elasticA, elasticB)
        else:
            raise ValueError("Both voigtA and voigtB must be provided if either is provided.")
        
        print("\nCSM properties:")
        mu = imt_multiplicity(crystA, crystB, slm)
        rms_strain = rmss(deformation_gradient(crystA, crystB, slm))
        print(f"\tmultiplicity: {mu:2d}")
        print(f"\tperiod: {zlcm * mu:2d}")
        print(f"\trmss: {100 * rms_strain:.1f} %")
        print(f"\trmsd: {d:.2f} Å")
        
        if strain == rmss:
            header = "# csm_id,  slm_id,      mu,  period,    rmss,  rmsd/Å"
            data = np.array([[0], [0], [mu], [zlcm * mu], [rms_strain], [d]]).T
        else:
            w = strain(deformation_gradient(crystA, crystB, slm))
            header = "# csm_id,  slm_id,      mu,  period,       w,    rmss,  rmsd/Å"
            data = np.array([[0], [0], [mu], [zlcm * mu], [w], [rms_strain], [d]]).T
        
        print("Using the primitive cells above, the SLM has the following IMT representation:")
        print("H_A =")
        print(slm[0])
        print("H_B =")
        print(slm[1])
        print("Q =")
        print(slm[2])
        mu0 = imt_multiplicity(crystA, crystB, slm0)
        if mu0 == mu: print(f"The input POSCAR files are already smallest supercells required to describe the CSM, with:")
        else: print(f"The input POSCAR files are {mu0/mu:d}-fold larger than the smallest supercells required to describe the CSM, with:")
        _, mA = decompose_cryst(crystA_sup, cryst_prim=crystA, tol=tol)
        _, mB = decompose_cryst(crystB_sup, cryst_prim=crystB, tol=tol)
        print("M_A =")
        print(mA)
        print("M_B =")
        print(mB)
        
        print(f"\nTo produce a list of CSMs that contains this CSM, run:\n"
                + f"\n\t$ crystmatch --enumerate '{fileA}' '{fileB}' {mu} {0.01 + rms_strain:.2f} --all {d:.2f}")
        print(f"\nIf you are not using a remote server and want to visualize this CSM, run:\n"
                + f"\n\t$ crystmatch --direct '{fileA}' '{fileB}'{' --literal' if args.literal else ''} --interact")
        
        if args.all is not None:
            max_d = args.all
            jobname += f"-d{max_d:.2f}"
            if max_d > 2.0: print(f"Warning: Current MAX_D = {max_d:.2f} may result in a large number of CSMs, which may take a very long time to enumerate.")
            if zlcm * mu0 >= 12: print(f"Warning: Current period {zlcm * mu0:d} may result in a large number of CSMs, which may take a very long time to enumerate.")
            print(f"\nEnumerating all CSMs with the same SLM as the input POSCAR files (mu = {mu0:d}) and MAX_D = {max_d:.2f} ...")
            pctlist, dlist = enumerate_pct(crystA, crystB, slm0, max_d, weight_func=weight_func, verbose=False, warning_threshold=100000)
            pct_arrs = [NPZ_ARR_COMMENT] + [np.array([], dtype=int).reshape(0, mu * zlcm, 4) for mu in range(1, mu0)]
            pct_arrs.append(pctlist)

    print('\n', end='')   # simply a line break

    # orientation analysis
    
    if args.orientation != None:
        vix, viy, viz, vfx, vfy, vfz, wix, wiy, wiz, wfx, wfy, wfz = args.orientation
        print(f"Comparing each CSM with the given OR:\n\t[{vix},{viy},{viz}] || [{vfx},{vfy},{vfz}]"
                + f"\n\t[{wix},{wiy},{wiz}] || [{wfx},{wfy},{wfz}]")
        r = orient_matrix([vix,viy,viz], [vfx,vfy,vfz], [wix,wiy,wiz], [wfx,wfy,wfz])
        anglelist = deviation_angle(crystA, crystB, slmlist, r, uspfix=args.uspfix)
        table = np.hstack((table, anglelist.reshape(-1,1)))
        if args.direct is not None: print(f"The deviation angle is {anglelist[0]:.6f} = {anglelist[0]*180/np.pi:.4f}°\n")

    # saving NPZ and POSCAR files
    
    if args.enumerate is not None:
        # saving enumeration results in NPZ format
        metadata = np.array([jobname, crystA, crystB, max_mu, max_strain, tol], dtype=object)
        np.savez(unique_filename("Saving CSMs to", f"CSM_LIST-{jobname}.npz"), *csm_bins, metadata=metadata, slmlist=slmlist, table=table)
    
    if args.direct is not None:
        # saving primitive cells
        dirname = unique_filename(f"Saving optimized final structure (rotation-free, translated to minimize RMSD) to", "CSM_single")
        makedirs(dirname)
        save_poscar(f"{dirname}{sep}{args.initial.split(sep)[-1]}", crystA_sup)
        save_poscar(f"{dirname}{sep}{splitext(args.final)[0].split(sep)[-1]}-optimized{splitext(args.final)[1]}", (cB_sup_final.T, crystB_sup[1], pB_sup_final.T))
    
    if args.poscar is not None:
        if mode == 'direct': pass
        
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
            crystA_sup, crystB_sup = csm_to_cryst(crystA, crystB, slm, p, ks)
            dirname = unique_filename(f"Saving the CSM with index {i:d} to", f"CSM_{i:d}")
            makedirs(dirname)
            save_poscar(f"{dirname}{sep}POSCAR_I", crystA_sup, crystname=f"CSM_{i:d} initial")
            save_poscar(f"{dirname}{sep}POSCAR_F", crystB_sup, crystname=f"CSM_{i:d} final")
            if args.interpolate is not None:
                print(f"Creating XDATCAR file with {args.interpolate} additional images ...")
                save_interpolation(f"{dirname}{sep}XDATCAR", crystA_sup, crystB_sup, args.interpolate, crystname=f"CSM_{i:d}")
    
    # saving CSV file and plotting

    if args.csv:
        if args.orientation is None:
            np.savetxt(unique_filename("Saving results to", f"TABLE-{jobname}.csv"), table * np.array([1,1,1,100,1]), fmt='%8d,%8d,%8d,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,    rmss,  rmsd/Å')
        else:
            np.savetxt(unique_filename("Saving results to", f"TABLE-OR-{jobname}.csv"), table * np.array([1,1,1,100,1,180/np.pi]), fmt='%8d,%8d,%8d,%8.4f,%8.4f,%8.4f',
                        header=' index,  slm_id,      mu,    rmss,  rmsd/Å, theta/°')
    
    if args.plot:
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        visualize_slmlist(unique_filename("Creating scatter plot in", f"PLOT-{jobname}.pdf"), table[:,3], table[:,4], table[:,2].round().astype(int),
                    cbarlabel=f"Multiplicity (× {zlcm:.0f} atom{'s' if zlcm>1 else ''})")
        if args.orientation is not None:
            visualize_slmlist(unique_filename("Creating scatter plot in", f"OR-{jobname}.pdf"),
                        table[:,3], table[:,4], table[:,5]*180/np.pi, cbarlabel="Deviation (°)", cmap=plt.cm.get_cmap('jet'))
    
    print(f"\nTotal time spent: {perf_counter()-time0:.2f} seconds")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n\nError: {e} Use '--help' for usage information, or see the documentation at \n\n\t" + __url__)