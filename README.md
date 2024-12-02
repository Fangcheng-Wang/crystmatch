# Crystmatch

If you use this code in your research, please cite the following paper:

\[1\] [F.-C. Wang et al., *Physical Review Letters* **132**, 086101 (2024)](https://arxiv.org/abs/2305.05278)

You are also welcome to contact me at `wfc@pku.edu.cn` for any questions or comments.

## Introduction

A solid-solid phase transition establishes an *atom-to-atom correspondence* between crystal structures $\mathcal A$ and $\mathcal B$. Such correspondence is called a *crystal-structure match* (CSM) [[1]](https://arxiv.org/abs/2305.05278). A CSM can be described by a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, which specifies how the lattice deforms from $\mathcal A$ to $\mathcal B$ and the correspondence between atoms in a supercell of $\mathcal A$ and those in $\mathcal B$.

The main functions of `crystmatch` are as follows:

- **Enumeration**:
    - Provide a complete list of *representative* [[1]](https://arxiv.org/abs/2305.05278) CSMs between two given crystal structures, with user-specified upper bounds on the multiplicity [[1]](https://arxiv.org/abs/2305.05278) and root-mean-square strain (RMSS).
    - (In progress) Provide a complete list of CSMs with user-specified upper bounds on the multiplicity, RMSS, and root-mean-square displacement (RMSD).

- **Analysis**:
    - Read a CSM from a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, and save CSMs in the same format.
    - Score CSMs by RMSS and RMSD.
    - Benchmark a CSM by its deviation angle from an orientation relationship (OR).

## Installation

Make sure you have **Python 3.9 or later** installed. You can check it by running:

```
$ python3 --version
```

Clone this repository and navigate to the directory where `setup.py` is located, run:

```
$ pip3 install .
```

Check whether `crystmatch` is successfully installed:

```
$ crystmatch --version
```

## Usage

To run `crystmatch`, one of the following modes must be selected:

1. **Enumeration mode**: Use `-E` or `--enumerate` to generate a list of CSMs, then perform preliminary analysis.
2. **Analysis mode**: Use `-A` or `--analyze` to read CSM(s) from [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files (or `CSM_LIST.npz` if provided), then perform detailed analysis. If `CSM_LIST.npz` is provided, you can export specific CSMs using `--export index1 [index2 ...]`.

The initial and final crystal structures should be specified in the form of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, unless in the analysis mode and `CSM_LIST.npz` is provided.

If you don't specify any mode or crystal structures, `crystmatch` will ask for input.

To see all available options, run:

```
$ crystmatch --help
```

## Examples

### Enumerating CSMs

To generate a list of representative [[1]](https://arxiv.org/abs/2305.05278) CSMs between two crystal structures stored in `./fcc` and `./bcc`, with a maximum multiplicity of `4` and a maximum RMSS of `0.2`:

```
$ crystmatch --initial fcc --final bcc --enumeration 4 -0.2
```

> We **strongly recommend** you to try small multiplicity (`2` or `4`) and RMSS between `0.2` and `0.5` first, and then gradually adjust these upper bounds to obtain desired results. Otherwise, the enumeration may take a very long time, or find no CSMs at all.

The following files will be created in the `fcc-bcc-m4s0.20` directory:

```
fcc-bcc-m4s0.20/
├── CSM_LIST(fcc-bcc-m4s0.20).npz       # stores the enumerated CSMs and metadata.
├── rmsd-rmss-mu(fcc-bcc-m4s0.20).pdf   # shows the RMSD-RMSS distribution of the CSMs.
└── TABLE(fcc-bcc-m4s0.20).csv          # organizes the multiplicity, RMSS, and RMSD of each CSM.
```

### Exporting CSMs from an NPZ file

To export the CSMs with indices `7` and `10` in `TABLE(foo).csv` from `CSM_LIST(foo).npz`:

```
$ crystmatch --analysis CSM_LIST(foo).npz --export 7 10
```

Two folders will be created in the `foo` directory, each containing a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files representing the CSM as:

```
foo/
├── CSM_7(m3s0.15d0.83)/
│   ├── POSCAR_I
│   └── POSCAR_F
└── CSM_10(m4s0.11d0.90)/
    ├── POSCAR_I
    └── POSCAR_F
```

The values in the brackets represent the basic attributes of CSMs. For example, `m3s0.15d0.83` indicates a multiplicity of `3`, RMSS of `0.15`, and RMSD of `0.83`.

### Orientation relationship analysis

To benchmark CSMs in `CSM_LIST(foo).npz` by their deviation angles from the OR $(111)_A\parallel(110)_B,[1\bar{1}0]_A\parallel[001]_B$:

```
$ crystmatch --analysis CSM_LIST(foo).npz --orientation 1 1 1 1 1 0 1 -1 0 0 0 1
```

The arguments after `--orientation` must be **cartesian coordinates**.

The ORs are determined via the rotation-free manner by default, and you can also use `--fixusp` to determine ORs via the USF-fixed manner; see Ref. [[1]](https://arxiv.org/abs/2305.05278) for their definitions.

## Python API

You can also use `crystmatch` as a Python module in your own Python scripts. See the [documentation](https://fangcheng-wang.github.io/crystmatch/) for details.