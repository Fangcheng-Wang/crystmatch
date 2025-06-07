# Tutorial

!!! abstract

    A solid-solid phase transition establishes an **atom-to-atom correspondence** between the initial and final crystal structures $\mathcal A$ and $\mathcal B$. Such correspondence is called a **crystal-structure match (CSM)**. A CSM can be described by a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, which specifies how the lattice deforms from $\mathcal A$ to $\mathcal B$ and the correspondence between atoms in a supercell of $\mathcal A$ and those in $\mathcal B$.

    The main functions of `crystmatch` are as follows:

    - **Enumeration**:
        - Provide a complete list of **representative** CSMs between $\mathcal A$ and $\mathcal B$, with user-specified upper bounds on the **multiplicity** and **strain**.
        - Provide a complete list of CSMs with user-specified upper bounds on the **multiplicity**, **strain**, and **shuffle distance**.

    - **Analysis**:
        - Read a CSM from a pair of POSCAR files, and save CSMs in the same format.
        - Calculate the root-mean-squared strain (RMSS), estimated strain energy density, and shuffle distance (RMSD) for each CSM.
        - Benchmark each CSM by its deviation angle from a given orientation relationship.
        - Visualize the distribution of strain, shuffle distance, and multiplicity of CSMs in a 2D scatter plot.
        - Visualize a CSM in a 3D interactive plot.
    
    Congruent CSMs (those differ only by a space-group transformation) are identified and excluded from the enumeration using the **[Spglib](https://spglib.readthedocs.io/en/stable/python-interface.html)** library by [Atsushi Togo *et al.*](https://www.tandfonline.com/doi/full/10.1080/27660400.2024.2384822)

## Installation

Make sure you have **Python 3.9 or later** installed. You can check it by running:

```
python3 --version
```

To install the latest version of `crystmatch`, run:

```
pip3 install --upgrade crystmatch
```


Check whether `crystmatch` is successfully installed:

```
crystmatch --version
```

!!! tip
    If you prefer using `conda`, you can install `crystmatch` by running:
    ```
    conda install -c conda-forge crystmatch
    ```
    and update it by running:
    ```
    conda update -c conda-forge crystmatch
    ```

## How to cite

If you use `crystmatch` in your research, please cite one of the following paper:

- **[Crystal-Structure Matches in Solid-Solid Phase Transitions](https://arxiv.org/abs/2305.05278)**

    *Physical Review Letters* **132**, 086101 (2024)

    ```
    @article{wang2024crystal,
        title={Crystal-Structure Matches in Solid-Solid Phase Transitions},
        author={Wang, Fang-Cheng and Ye, Qi-Jun and Zhu, Yu-Cheng and Li, Xin-Zheng},
        journal={Phys. Rev. Lett.},
        volume={132},
        number={8},
        pages={086101},
        year={2024},
        publisher={APS},
        doi={10.1103/PhysRevLett.132.086101}
    }
    ```

- **[Classification and Enumeration of Solid-Solid Phase Transition Mechanisms](https://arxiv.org/abs/2506.05105)**

    *Under review* (2025)

    ```
    @unpublished{wang2025classification,
        title={Classification and Enumeration of Solid-Solid Phase Transition Mechanisms},
        author={Wang, Fang-Cheng and Ye, Qi-Jun and Zhu, Yu-Cheng and Li, Xin-Zheng},
        year={2025}
    }
    ```

## Usage

To run `crystmatch`, one of the following modes must be selected:

1. **`--enumerate`**
    
    Enumerate a list of representative CSMs, save them to a `CSMLIST.npz` file, and perform preliminary analysis. The initial and final crystal structures must be specified in two separate [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files.

2. **`--read`**
    
    Read CSMs from a `CSMLIST.npz` file. You can export specific CSMs to POSCARs or XDATCARs, perform orientation-relationship analysis, and visualize them interactively.

3. **`--direct`**
    
    Directly determine a single CSM by two POSCAR files (must have the same number of atoms) and perform detailed analysis.

**We strongly recommend starting with the [examples](#examples) provided below.** To see [all available options](https://fangcheng-wang.github.io/crystmatch/cli/), run:

```
crystmatch --help
```

## Examples

### Enumerating CSMs

To generate a list of representative CSMs between two crystal structures given by `graphite.poscar` and `diamond.poscar`, with upper bounds `MAX_MU = 2` for multiplicity and `MAX_STRAIN = 0.4` for RMSS, run:

```
crystmatch --enumerate graphite.poscar diamond.poscar 2 0.4
```

The following files will be created in the current directory:

```powershell
./
├── CSMLIST-m2s0.4.npz       # stores the enumerated CSMs and metadata.
├── SUMMARY-m2s0.4.csv       # lists the multiplicity, RMSS, and RMSD.
└── SCATTER-m2s0.4.pdf       # shows the RMSD-RMSS distribution of the CSMs.
```

!!! warning "Caveat"
    We recommend you to **try `MAX_MU = 2` and `MAX_W = 0.4` first**, and then gradually adjust these upper bounds, usually by **increasing `MAX_MU` and decreasing `MAX_W`**, to obtain desired results. Otherwise, the enumeration may take a very long time, or find no CSMs at all.

### Reading and exporting CSMs from an NPZ file

After enumeration, you can see the properties of CSMs in the CSV file, which also contains their indices in the NPZ file. If you want to export the CSMs with indices `7` and `10` in `CSMLIST-foo.npz`, run:

```
crystmatch --read CSMLIST-foo.npz --export 7 10
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

To benchmark CSMs in `CSMLIST-foo.npz` by their deviation angles from the OR $(111)_A\parallel(110)_B,[1\bar{1}0]_A\parallel[001]_B$, run:

```
crystmatch --read CSMLIST-foo.npz --orientation 1 1 1 1 1 0 1 -1 0 0 0 1
```

Note that the arguments after `--orientation` must be **Cartesian coordinates**.

The ORs are determined via the *rotation-free*[^1] manner by default, and you can also use `--uspfix` to determine ORs via the *USF-fixed*[^1] manner.

### Single-CSM analysis

To analyze a single CSM defined by two POSCAR files, run:

```
crystmatch --initial POSCAR1 --final POSCAR2 --single
```

`crystmatch` will also save the rigid-transformation optimized (with rotation-free orientation and RMSD-minimized overall position) CSM in the current directory like this:

```
./
└── CSM_single/
    ├── POSCAR1
    └── POSCAR2-optimized
```
