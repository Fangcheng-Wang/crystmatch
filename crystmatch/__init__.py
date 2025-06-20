# Copyright (C) 2024 Fangcheng Wang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Enumerating and analyzing crystal-structure matches."""

from .utilities import *
from .enumeration import *
from .analysis import *
from time import perf_counter
import sys
import argparse

__name__ = "crystmatch"
__version__ = "2.0.6"
__author__ = "Fang-Cheng Wang"
__email__ = "wfc@pku.edu.cn"
__description__ = "Enumerating and analyzing crystal-structure matches for solid-solid phase transitions."
__url__ = "https://fangcheng-wang.github.io/crystmatch/"
__epilog__ = "The current version is v" + __version__ + ". To get the latest version, please run:\
\n\n\t$ pip3 install --upgrade crystmatch\n\nWe also recommend you to see the documentation at:\
\n\n\t" + __url__ + "\n\nIf you use crystmatch in your research, please cite one of the following paper:\
\n\n\t[1] FC Wang, QJ Ye, YC Zhu, and XZ Li, Physical Review Letters 132, 086101 (2024) (https://arxiv.org/abs/2305.05278)\
\n\t[2] FC Wang, QJ Ye, YC Zhu, and XZ Li, arXiv:2506.05105 (2025) (https://arxiv.org/abs/2506.05105)\
\n\nYou are also welcome to contact us at " + __email__ + " for any questions, feedbacks or comments."

def main():

    time0 = perf_counter()
    sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)
    sys.stderr = open(sys.stdout.fileno(), mode='w', buffering=1, encoding='utf-8', closefd=False)

    parser = argparse.ArgumentParser(
        prog = __name__,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        usage = "%(prog)s [-h] [-V] [-E POSCAR_I POSCAR_F MAX_MU MAX_STRAIN] [-R CSMLIST [IND1 IND2 ...]] [-D POSCAR_I POSCAR_F] [-n IMAGES] [-r]"\
                + "[-e CSMCAR] [-a MAX_D] [-t TOL] [-i [SIZE]] [-o ASSUM] [-p [ASSUM]] [-x [ASSUM]] [-m] [-c] [-s]",
        description = __description__,
        epilog = __epilog__)
    
    parser.add_argument("-V", "-v", "--version", action='version', version="%(prog)s " + __version__)
    parser.add_argument("-E", "--enumerate", nargs=4, type=str, metavar=('POSCAR_I', 'POSCAR_F', 'MAX_MU', 'MAX_STRAIN'),
                        help="enumerate representative CSMs between POSCAR_I and POSCAR_F, with MAX_MU and MAX_STRAIN as upper bounds for the multiplicity and RMSS")
    parser.add_argument("-R", "--read", nargs='+', type=str, metavar=('CSMLIST', 'IND1 IND2'),
                        help="read CSMs from CSMLIST, with optional IND1 IND2 ... to specify the indices of the CSMs to be loaded")
    parser.add_argument("-D", "--direct", nargs=2, type=str, metavar=('POSCAR_I', 'POSCAR_F'),
                        help="directly read a single CSM from POSCAR_I and POSCAR_F")
    parser.add_argument("-e", "--extra", type=str, metavar='CSMCAR',
                        help="load extra parameters from CSMCAR, including the elastic tensors, atomic weights, and orientation relationship used to benchmark CSMs")
    parser.add_argument("-a", "--all", type=float, metavar='MAX_D',
                        help="enumerate all CSMs for each SLM, with MAX_D as the upper bound for the shuffle distance")
    parser.add_argument("-n", "--nebmake", type=int, metavar='IMAGES',
                        help="when using --direct, interpolate between POSCAR_I and POSCAR_F to generate images for NEB calculation")
    parser.add_argument("-r", "--restore", action='store_true',
                        help="when using --direct, add integers to fractional coordinates to minimize the shuffle distance")
    parser.add_argument("-t", "--tol", type=float, default=1e-3, metavar='TOL',
                        help="tolerance in angstroms for detecting symmetry (default is 1e-3)")
    parser.add_argument("-i", "--interact", nargs='?', type=float, default=None, const=1.2, metavar='SIZE',
                        help="interactively visualize each CSM using a 3D plot, with SIZE controlling the radius of the cluster to display (default is 1.2)")
    parser.add_argument("-o", "--orientation", type=str, choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="calculate the deviation angle from the given orientation relationship for each CSM (must be used with --extra)")
    parser.add_argument("-p", "--poscar", nargs='?', type=str, default=None, const='norot', choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="export each enumerated/read CSM as a pair of POSCAR files (default is 'norot')")
    parser.add_argument("-x", "--xdatcar", nargs='?', type=str, default=None, const='norot', choices=['norot', 'uspfixed'], metavar='ASSUM',
                        help="export each enumerated/read CSM as an XDATCAR file with 100 frames (default is 'norot')")
    parser.add_argument("-m", "--mediumcell", action='store_true',
                        help="when using --poscar or --xdatcar, output the average cell instead of the initial and final cells")
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
    
    # check if --read has duplicate indices
    if args.read is not None and len(set([int(i) for i in args.read[1:]])) != len(args.read[1:]): raise ValueError("Duplicate indices in --read.")
    
    # check if --nebmake and --restore is used with --direct without --all
    if args.nebmake is not None:
        if not (mode == 'direct' and args.all is None): raise ValueError("'--nebmake' can only be used with '--direct' without '--all'.")
    if args.restore:
        if not (mode == 'direct' and args.all is None): raise ValueError("'--restore' can only be used with '--direct' without '--all'.")
    
    # check if --mediumcell is used with (--poscar or --xdatcar) and ASSUM == 'norot'
    if args.mediumcell:
        if args.poscar is None and args.xdatcar is None: raise ValueError("'--mediumcell' can only be used with '--poscar' or '--xdatcar'.")
        if args.poscar == 'uspfixed' or args.xdatcar == 'uspfixed': raise ValueError("'--mediumcell' can only be used with ASSUM == 'norot'.")
    
    # automatically enable --csv and --scatter
    if args.enumerate is not None or args.all is not None or args.orientation is not None:
        args.csv = True
        args.scatter = True
    
    # load global settings
    voigtA, voigtB, weight_func, ori_rel = None, None, None, None
    if args.extra is not None: voigtA, voigtB, weight_func, ori_rel = load_csmcar(args.extra)
    tol = args.tol
    
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
        print("A CSM with multiplicity μ contains ({zlcm:d} * μ) atoms.")

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
            print(f"Warning: Current MAX_MU = {max_mu:d} and MAX_STRAIN = {max_strain:.2f} may result in a large number of SLMs, which may take a very long time to enumerate.")
        
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
            header = ['csm_id', 'slm_id', 'mu', 'period', 'rmss', 'd']
            data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, dlist]).T
        else:
            rmsslist = rmss(deformation_gradient(crystA, crystB, slmlist))[slm_ind]
            header = ['csm_id', 'slm_id', 'mu', 'period', 'w', 'rmss', 'd']
            data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, rmsslist, dlist]).T

    elif mode == 'read':
        
        jobname = 'read'
        filenpz = args.read[0]
        crystA, crystB, slmlist, slm_ind, pct_arrs, table = load_npz(filenpz, verbose=True)
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))
        print('\n', end='')
        max_mu = imt_multiplicity(crystA, crystB, slmlist[slm_ind]).max()
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        header = table.header
        data = table.data
        
        if len(args.read) > 1:
            indices = np.array([int(i) for i in args.read[1:]], dtype=int)
            valid = indices < len(slm_ind)
            if not valid.all():
                print(f"Warning: Indices out of range (> {len(slm_ind)-1:d}) will be ignored.")
                indices = indices[valid]
            print(f"Reading CSMs with indices:", end='')
            max_mu = imt_multiplicity(crystA, crystB, slmlist[slm_ind[indices]]).max()
            slmlist_temp = np.array([], dtype=int).reshape(0,3,3,3)
            slm_ind_temp = []
            pct_arrs_temp = [NPZ_ARR_COMMENT] + [np.array([], dtype=int).reshape(0, mu * zlcm, 4) for mu in range(1, max_mu + 1)]
            for i in indices:
                print(f" {i:4d}", end='')
                slm, p, ks = unzip_csm(i, crystA, crystB, slmlist, slm_ind, pct_arrs)
                where = np.nonzero((slm == slmlist_temp).all(axis=(1,2,3)))[0]
                if where.shape[0] == 0:
                    slm_ind_temp.append(slmlist_temp.shape[0])
                    slmlist_temp = np.concatenate([slmlist_temp, slm.reshape(1,3,3,3)], axis=0)
                else:
                    slm_ind_temp.append(where[0])
                mu = imt_multiplicity(crystA, crystB, slm)
                pct_arrs_temp[mu] = np.concatenate([pct_arrs_temp[mu], zip_pct(p, ks).reshape(1,-1,4)], axis=0)
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
                header = ['csm_id', 'slm_id', 'mu', 'period', 'rmss', 'd']
                data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, dlist]).T
            else:
                rmsslist = rmss(deformation_gradient(crystA, crystB, slmlist))[slm_ind]
                header = ['csm_id', 'slm_id', 'mu', 'period', 'w', 'rmss', 'd']
                data = np.array([np.arange(len(slm_ind)), slm_ind, mulist, zlcm * mulist, strainlist, rmsslist, dlist]).T
        
        else: strain = 'w' if 'w' in header else rmss

    elif mode == 'direct':

        jobname = 'direct'
        fileA = args.direct[0]
        fileB = args.direct[1]
        crystA_sup = load_poscar(fileA, tol=tol, to_primitive=False)
        crystB_sup = load_poscar(fileB, tol=tol, to_primitive=False)
        if not crystA_sup[1].shape[0] == crystB_sup[1].shape[0]: raise ValueError("The numbers of atoms in two POSCAR files are not equal!")
        check_stoichiometry(crystA_sup[1], crystB_sup[1])
        crystA, crystB, slm0, p0, ks0 = cryst_to_csm(crystA_sup, crystB_sup, tol=tol)
        zlcm = np.lcm(len(crystA[1]), len(crystB[1]))
        print(save_poscar(None, crystA, crystname="\nPrimitive cell of the initial structure (in POSCAR format):"))
        print(save_poscar(None, crystB, crystname="\nPrimitive cell of the final structure (in POSCAR format):"))

        if args.restore:
            d0 = csm_distance(crystA, crystB, slm0, p0, ks0, weight_func=weight_func)
            ks0 = optimize_ct(crystA, crystB, slm0, p0, weight_func=weight_func)[1]
        slm, p, ks = primitive_shuffle(crystA, crystB, slm0, p0, ks0)
        d = csm_distance(crystA, crystB, slm, p, ks, weight_func=weight_func)
        if args.restore and d < d0 - tol:
            print(f"\nBy adding integers to fractional coordinates (--restore), the shuffle distance is reduced by {d0 - d:.4f} Å.")
        
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
        ps = la.svd(deformation_gradient(crystA, crystB, slm), compute_uv=False) - 1
        rms_strain = rmss(deformation_gradient(crystA, crystB, slm))
        print(f"\tmultiplicity (μ): {mu:d}")
        print(f"\tperiod: {zlcm * mu:d}")
        print(f"\tprincipal strains: {100 * ps[0]:.2f} %, {100 * ps[1]:.2f} %, {100 * ps[2]:.2f} %")
        print(f"\troot-mean-squared strain (RMSS): {100 * rms_strain:.2f} %")
        print(f"\tshuffle distance (RMSD): {d:.4f} Å")
        if not strain == rmss:
            w = strain(deformation_gradient(crystA, crystB, slm))
            print(f"\testimated strain energy: {w:.3f} (same unit as in CSMCAR)")
        
        print("\nUsing the primitive cells above, the SLM has the following IMT representation:")
        print("H_A =")
        print(slm[0])
        print("H_B =")
        print(slm[1])
        print("Q =")
        print(slm[2])
        mu0 = imt_multiplicity(crystA, crystB, slm0)
        if mu0 == mu: print(f"\nThe input POSCAR files are already smallest supercells required to describe the CSM, with:")
        else: print(f"\nThe input POSCAR files are {mu0//mu:d}-fold larger than the smallest supercells required to describe the CSM, with:")
        print("M_A =")
        print(decompose_cryst(crystA_sup, cryst_prim=crystA, tol=tol)[1])
        print("M_B =")
        print(decompose_cryst(crystB_sup, cryst_prim=crystB, tol=tol)[1])
        
        print(f"\nTo produce a list of CSMs that contains this CSM, run:\n"
                + f"\n\t$ crystmatch --enumerate '{fileA}' '{fileB}' {mu} {0.01 + rms_strain:.2f} --all {0.01 + d:.2f}")
        print(f"\nIf you are not using a remote server and want to visualize this CSM, run:\n"
                + f"\n\t$ crystmatch --direct '{fileA}' '{fileB}' --interact")
        
        if args.all is None:
            dlist = np.array([d])
            mulist = np.array([mu], dtype=int)
            slmlist = np.array([slm], dtype=int)
            slm_ind = np.array([0], dtype=int)
            pct_arrs = [NPZ_ARR_COMMENT] + [np.array([], dtype=int).reshape(0, m * zlcm, 4) for m in range(1, mu)]
            pct_arrs.append(zip_pct(p, ks).reshape(1, -1, 4))
        else:
            max_d = args.all
            jobname += f"-d{max_d:.2f}"
            if max_d > 2.0: print(f"Warning: Current MAX_D = {max_d:.2f} may result in a large number of CSMs, which may take a very long time to enumerate.")
            if zlcm * mu0 >= 12: print(f"Warning: Current period {zlcm * mu0:d} may result in a large number of CSMs, which may take a very long time to enumerate.")
            print(f"\nEnumerating all CSMs with the same SLM as the input POSCAR files (μ = {mu0:d}) and MAX_D = {max_d:.2f} ...")
            pctlist, dlist = enumerate_pct(crystA, crystB, slm0, max_d, weight_func=weight_func, verbose=False, warning_threshold=100000)
            mulist = []
            slmlist = np.array([], dtype=int).reshape(0,3,3,3)
            slm_ind = []
            pct_arrs = [NPZ_ARR_COMMENT] + [np.array([], dtype=int).reshape(0, mu * zlcm, 4) for mu in range(1, mu0 + 1)]
            for pct in pctlist:
                p_tmp, ks_tmp = unzip_pct(pct)
                slm1, p1, ks1 = primitive_shuffle(crystA, crystB, slm0, p_tmp, ks_tmp)
                mu1 = imt_multiplicity(crystA, crystB, slm1)
                mulist.append(mu1)
                where = np.nonzero((slm1 == slmlist).all(axis=(1,2,3)))[0]
                if where.shape[0] == 0:
                    slm_ind.append(slmlist.shape[0])
                    slmlist = np.concatenate((slmlist, slm1.reshape(-1,3,3,3)), axis=0)
                else:
                    slm_ind.append(where[0])
                pct_arrs[mu1] = np.concatenate([pct_arrs[mu1], zip_pct(p1, ks1).reshape(1,-1,4)], axis=0)
            mulist = np.array(mulist, dtype=int)
            slm_ind = np.array(slm_ind, dtype=int)

        if strain == rmss:
            header = ['csm_id', 'slm_id', 'mu', 'period', 'rmss', 'd']
            data = np.array([np.arange(len(dlist)), slm_ind, mulist, zlcm * mulist, [rms_strain] * len(dlist), dlist]).T
        else:
            header = ['csm_id', 'slm_id', 'mu', 'period', 'w', 'rmss', 'd']
            data = np.array([np.arange(len(dlist)), slm_ind, mulist, zlcm * mulist, [w] * len(dlist), [rms_strain] * len(dlist), dlist]).T

    print('\n', end='')   # simply a line break

    # orientation analysis
    
    if args.orientation is not None:
        if ori_rel is None:
            print("Warning: Orientation relationship is not provided. Skipping orientation analysis.")
        else:
            print("Comparing each CSM with the given OR ...")
            hkl_i1, hkl_f1 = ori_rel[0]
            hkl_i2, hkl_f2 = ori_rel[1]
            vi = miller_to_vec(crystA, hkl_i1, tol=tol)
            vf = miller_to_vec(crystB, hkl_f1, tol=tol)
            wi = miller_to_vec(crystA, hkl_i2, tol=tol)
            wf = miller_to_vec(crystB, hkl_f2, tol=tol)
            r = orient_matrix(vi, vf, wi, wf)
            anglelist = deviation_angle(crystA, crystB, slmlist, r, orientation=args.orientation)[slm_ind]
            header = header + ['angle']
            data = np.hstack([data, anglelist.reshape(-1,1)])

    # save NPZ, CSV, POSCAR, and XDATCAR files
    
    table = Table(data, header)
    if mode == 'enumerate' or args.all is not None:
        save_npz(unique_filename("Saving SLMs and CSMs to", f"CSMLIST-{jobname}.npz"), crystA, crystB, slmlist, slm_ind, pct_arrs, table)
    if args.csv: save_csv(unique_filename("Saving summary table to", f"SUMMARY-{jobname}.csv"), table)
    
    if args.poscar is not None or args.xdatcar is not None:
        direxport = unique_filename(f"Saving POSCAR and/or XDATCAR files to", f"EXPORT-{jobname}", ext=False)
        makedirs(direxport)
        for i in range(len(slm_ind)):
            slm, p, ks = unzip_csm(i, crystA, crystB, slmlist, slm_ind, pct_arrs)
            if mode == 'read' and len(args.read) > 1: ind = indices[i]
            else: ind = i
            makedirs(f"{direxport}{sep}CSM_{ind:d}")
            if args.poscar == 'norot' or args.xdatcar == 'norot':
                crystA_csm, crystB_csm_norot = csm_to_cryst(crystA, crystB, slm, p, ks, tol=tol, orientation='norot', use_medium_cell=args.mediumcell,
                                                            min_t0=True, weight_func=weight_func)
                if args.poscar == 'norot':
                    save_poscar(f"{direxport}{sep}CSM_{ind:d}{sep}POSCAR_I", crystA_csm, crystname=f"CSM_{ind:d} initial")
                    save_poscar(f"{direxport}{sep}CSM_{ind:d}{sep}POSCAR_F", crystB_csm_norot, crystname=f"CSM_{ind:d} final (no rotation)")
                if args.xdatcar == 'norot':
                    save_xdatcar(f"{direxport}{sep}CSM_{ind:d}{sep}XDATCAR", crystA_csm, crystB_csm_norot, n_im=50, crystname=f"CSM_{ind:d} (no rotation)")
            if args.poscar == 'uspfixed' or args.xdatcar == 'uspfixed':
                crystA_csm, crystB_csm_usp1, crystB_csm_usp2 = csm_to_cryst(crystA, crystB, slm, p, ks, tol=tol, orientation='uspfixed',
                                                                            min_t0=True, weight_func=weight_func)
                if args.poscar == 'uspfixed':
                    save_poscar(f"{direxport}{sep}CSM_{ind:d}{sep}POSCAR_I", crystA_csm, crystname=f"CSM_{ind:d} initial")
                    save_poscar(f"{direxport}{sep}CSM_{ind:d}{sep}POSCAR_F1", crystB_csm_usp1, crystname=f"CSM_{ind:d} final (USP1 fixed)")
                    save_poscar(f"{direxport}{sep}CSM_{ind:d}{sep}POSCAR_F2", crystB_csm_usp2, crystname=f"CSM_{ind:d} final (USP2 fixed)")
                if args.xdatcar == 'uspfixed':
                    save_xdatcar(f"{direxport}{sep}CSM_{ind:d}{sep}XDATCAR1", crystA_csm, crystB_csm_usp1, n_im=50, crystname=f"CSM_{ind:d} (USP1 fixed)")
                    save_xdatcar(f"{direxport}{sep}CSM_{ind:d}{sep}XDATCAR2", crystA_csm, crystB_csm_usp2, n_im=50, crystname=f"CSM_{ind:d} (USP2 fixed)")
    
    if args.nebmake is not None:
        n_im = args.nebmake
        if n_im <= 0 or n_im >= 99:
            print("Warning: Invalid number of images for '--nebmake'. Skipping interpolation.")
        else:
            print(f"Interpolating between the initial and final structures with {n_im} images ...")
            for i in range(n_im+2):
                if exists(f"{i:02d}"):
                    raise ValueError(f"Please remove all folders named {', '.join([f'{j:02d}' for j in range(n_im+2)])} before using '--nebmake'.")
            for i, cryst_im in enumerate(nebmake(crystA_sup, crystB_sup, n_im)):
                makedirs(f"{i:02d}")
                save_poscar(f"{i:02d}{sep}POSCAR", cryst_im, crystname=f"crystmatch-nebmake (image {i:02d})")
    
    # create scatter plot

    if args.scatter:
        visualize_slmlist(unique_filename("Creating scatter plot in", f"SCATTER-{jobname}.pdf"), data[:,header.index('rmss' if strain == rmss else 'w')],
                            data[:,header.index('d')], data[:,header.index('mu')].round().astype(int),
                            cbarlabel=f"Multiplicity (× {zlcm:.0f} atoms)")
        if 'angle' in header:
            visualize_slmlist(unique_filename("Creating scatter plot in", f"ANGLE-{jobname}.pdf"), data[:,header.index('rmss' if strain == rmss else 'w')],
                            data[:,header.index('d')], data[:,header.index('angle')],
                            cbarlabel="Deviation angle (rad)", cmap=plt.cm.get_cmap('magma'))
        print(f"\tThe horizontal axis 'Strain' represents {'RMSS' if strain == rmss else 'estimated strain energy (with the same unit as in CSMCAR)'}.")
    
    # interactive visualization
    
    time1 = perf_counter()
    print(f"\nTotal time spent{' (excluding interactive visualization)' if args.interact is not None else ''}: {time1-time0:.2f} s")
    
    if args.interact is not None:
        for i in range(len(slm_ind)):
            slm, p, ks = unzip_csm(i, crystA, crystB, slmlist, slm_ind, pct_arrs)
            print(f"Displaying the CSM with csm_id = {indices[i] if (mode == 'read' and len(args.read) > 1) else i} ...")
            visualize_csm(crystA, crystB, slm, p, ks, weight_func=weight_func, tol=tol, cluster_size=args.interact, label = 
                            f"csm_id = {indices[i] if (mode == 'read' and len(args.read) > 1) else i}, "
                            + f"μ = {imt_multiplicity(crystA, crystB, slm):d}, "
                            + f"RMSS = {100 * rmss(deformation_gradient(crystA, crystB, slm)):.2f}%, "
                            + f"RMSD = {csm_distance(crystA, crystB, slm, p, ks, weight_func=weight_func):.3f}Å")

    return

if __name__ == "__main__":
    main()