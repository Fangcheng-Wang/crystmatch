# Tutorial

!!! abstract

    A solid-solid phase transition establishes an *atom-to-atom correspondence* between crystal structures $\mathcal A$ and $\mathcal B$. Such correspondence is called a *crystal-structure match* (CSM)[^1]. A CSM can be described by a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, which specifies how the lattice deforms from $\mathcal A$ to $\mathcal B$ and the correspondence between atoms in a supercell of $\mathcal A$ and those in $\mathcal B$.

    The main functions of `crystmatch` are as follows:

    - **Enumeration**:
        - Provide a complete list of *representative*[^1] CSMs between two given crystal structures, with user-specified upper bounds on the *multiplicity*[^1] and root-mean-square strain (RMSS).
        - (In progress) Provide a complete list of CSMs with user-specified upper bounds on the multiplicity, RMSS, and root-mean-square displacement (RMSD).

    - **Analysis**:
        - Read a CSM from a pair of POSCAR files, and save CSMs in the same format.
        - Score CSMs by RMSS and RMSD.
        - Benchmark a CSM by its deviation angle from an orientation relationship (OR).
    
    Congruent CSMs (those differ only by a space-group transformation) are identified and excluded from the enumeration using the `spglib` library by [Atsushi Togo *et al.*](https://www.tandfonline.com/doi/full/10.1080/27660400.2024.2384822)

## Installation

Make sure you have **Python 3.9 or later** installed. You can check it by running:

```
$ python3 --version
```

Clone [**this repository**](https://github.com/fangcheng-wang/crystmatch) and navigate to the directory where `setup.py` is located, run:

```
$ pip3 install .
```

Check whether `crystmatch` is successfully installed:

```
$ crystmatch --version
```

## Usage

To run `crystmatch`, one of the following modes must be selected:

1. **Enumeration mode**: Use `-E` or `--enumeration` to generate a list of CSMs, save them to a `CSM_LIST.npz` file, and perform preliminary analysis. The initial and final crystal structures must be specified in the form of POSCAR files.
2. **Read mode**: Use `-R` or `--read` to read CSMs from a `CSM_LIST.npz` file. You can export specific CSMs to POSCARs using `--export index1 [index2 ...]`, perform OR analysis using `--orientation` (see example below), generate CSV tables using `--csv`, and/or visualize the RMSD-RMSS-multiplicity distribution using `--plot`.
3. **Single-CSM mode**: Use `-S` or `--single` to determine a single CSM directly by two POSCAR files (must have the same number and species of atoms) and perform detailed analysis.

If you are confused about the usage, simply run:

```
$ crystmatch
```

and it will ask for input. To see all available options, run:

```
$ crystmatch --help
```

## Examples

### Enumerating CSMs

To generate a list of representative CSMs between two crystal structures stored in `./fcc` and `./bcc`, with a maximum multiplicity of `4` and a maximum RMSS of `0.2`:

```
$ crystmatch --initial fcc --final bcc --enumeration 4 0.2
```

!!! note
    We strongly recommend you to try small multiplicity (`2` or `4`) and RMSS between `0.2` and `0.5` first, and then gradually adjust these upper bounds to obtain desired results. Otherwise, the enumeration may take a very long time, or find no CSMs at all.

The following files will be created in the current directory:

```
./
├── CSM_LIST-fcc-bcc-m4s0.20.npz       # stores the enumerated CSMs and metadata.
├── PLOT-fcc-bcc-m4s0.20.pdf           # shows the RMSD-RMSS distribution of the CSMs.
└── TABLE-fcc-bcc-m4s0.20.csv          # shows the multiplicity, RMSS, and RMSD of each CSM.
```

### Exporting CSMs from an NPZ file

After enumeration, you can see the properties of CSMs in the CSV file, which also contains their indices in the NPZ file. If you want to export the CSMs with indices `7` and `10` in `CSM_LIST-foo.npz`, run:

```
$ crystmatch --read CSM_LIST-foo.npz --export 7 10
```

Two folders will be created in the current directory, each containing a pair of POSCAR files representing the CSM. The current directory will look like this:

```
./
├── CSM_7/
│   ├── POSCAR_I
│   └── POSCAR_F
└── CSM_10/
    ├── POSCAR_I
    └── POSCAR_F
```

### Orientation relationship analysis

To benchmark CSMs in `CSM_LIST-foo.npz` by their deviation angles from the OR $(111)_A\parallel(110)_B,[1\bar{1}0]_A\parallel[001]_B$, run:

```
$ crystmatch --read CSM_LIST-foo.npz --orientation 1 1 1 1 1 0 1 -1 0 0 0 1
```

Note that the arguments after `--orientation` must be **cartesian coordinates**.

The ORs are determined via the *rotation-free*[^1] manner by default, and you can also use `--uspfix` to determine ORs via the *USF-fixed*[^1] manner.

### Single-CSM analysis

To analyze a single CSM defined by two POSCAR files, run:

```
$ crystmatch --initial POSCAR1 --final POSCAR2 --single
``

[^1]: [FC Wang, QJ Ye, YC Zhu, and XZ Li, *Physical Review Letters* **132**, 086101 (2024)](https://arxiv.org/abs/2305.05278)