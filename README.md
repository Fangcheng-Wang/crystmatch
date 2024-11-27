If you use this code in your research, please cite our paper:

- [F.-C. Wang *et al.*, Physical Review Letters **132**, 086101 (2024)](https://doi.org/10.1103/PhysRevLett.132.086101)

# crystmatch

Enumerating and benchmarking crystal-structure matches for solid-solid phase transitions.

## Introduction

A solid-solid phase transition establishes an *atom-to-atom correspondence* between crystal structures $\mathcal A$ and $\mathcal B$. Such correspondence is called a crystal-structure match (CSM) [[1]](https://arxiv.org/abs/2305.05278). A CSM can be described by a pair of POSCAR files, which specifies how the lattice deforms from $\mathcal A$ to $\mathcal B$ and the correspondence between the atoms in supercells of $\mathcal A$ and $\mathcal B$.

The main goals of `crystmatch` are:

- Provide a comprehensive list of candidate CSMs between two given POSCAR files of $\mathcal A$ and $\mathcal B$.
- Score CSMs by root-mean-square strain (RMSS) and root-mean-square displacement (RMSD).
- Benchmark CSMs by the deviation from an orientation relationship (OR).
- Save CSMs as pairs of POSCAR files.

## Installation

In the directory where `setup.py` is located, run:

```bash
$ pip3 install .
```

Check whether `crystmatch` is successfully installed:

```bash
$ crystmatch --version
```

## Usage (as command-line tool)

### Commonly Used

- To generate a comprehensive list of CSMs (scored and plotted):

  ```bash
  $ crystmatch -A POSCAR_A -B POSCAR_B -E MU_MAX KAPPA_MAX
  ```

  Append `-p` to create POSCAR files for each CSM.

- To benchmark a list of SLMs by an OR, for example, $(111)_A\parallel(110)_B$ and $[1\overline{1}0]_A\parallel[001]_B$:

  ```bash
  $ crystmatch -A POSCAR_A -B POSCAR_B -L SLM_LIST -s SCORE_CSV -b 1 1 1 1 1 0 1 -1 0 0 0 1
  ```

  Please use *Cartesian coordinates* to specify the OR. Append `--fixusp` to determine OR with fixed uniformly scaled plane.

- To create POSCAR files for CSMs generated from `SLM_LIST`:

  ```bash
  $ crystmatch -A POSCAR_A -B POSCAR_B -L SLM_LIST --poscar
  ```

  Append `-i INDEX1 INDEX2 ...` to specify certain CSMs in the list.

- To get help:

  ```bash
  $ crystmatch --help
  ```

### Mandatory Arguments

Please specify one of the following modes:

1. `-E MU_MAX KAPPA_MAX`: Enumerate SLMs. Results are saved as `SLM_LIST.npy`
2. `-L SLM_LIST`: Load SLMs from file

Either mode needs to specify crystal structures by `-A POSCAR_A` and `-B POSCAR_B`

### Optional Arguments

- By default, `SCORE.csv` and `score.pdf` are created, which can be disabled by `--nocsv` and `--noplot`
- Change output directory: `-o OUTDIR`
- Load scores from file (only in `-L` mode): `-s SCORE`
- Load CSMs with specified indices (only in `-L` mode): `-i INDEX1 INDEX2 ...`
- Benchmark SLMs by an OR: `-b v1x v1y v1z w1x w1y w1z v2x v2y v2z w2x w2y w2z`
- Determine OR with fixed uniformly scaled plane (only use together with `-b`): `--fixusp`

## Usage (as Python package)

`crystmatch` is more flexible in Python scripts. Please see the [documentation (in progress)]().
