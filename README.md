# Crystmatch

## Introduction

A solid-solid phase transition establishes an **atom-to-atom correspondence** between the initial and final crystal structures $\mathcal A$ and $\mathcal B$. Such correspondence is called a **crystal-structure match (CSM)**. A CSM can be described by a pair of [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) files, which specifies how the lattice deforms from $\mathcal A$ to $\mathcal B$ and the correspondence between atoms in a supercell of $\mathcal A$ and those in $\mathcal B$.

The main functions of `crystmatch` are as follows:

- **Enumeration**:
    - Provide a complete list of **representative** CSMs between $\mathcal A$ and $\mathcal B$, with user-specified upper bounds on the **multiplicity** and **strain**.
    - Provide a complete list of CSMs with user-specified upper bounds on the **multiplicity**, **strain**, and **shuffle distance**.

- **Analysis**:
    - Read a CSM from a pair of POSCAR files, and save CSMs in the same format.
    - Calculate the root-mean-square strain (RMSS), estimated strain energy density, and shuffle distance (RMSD) for each CSM.
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
pip3 install --upgrade numpy scipy spglib numba matplotlib
pip3 install --upgrade crystmatch
```

Check whether `crystmatch` is successfully installed:

```
crystmatch --version
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
    @article{wang2025classification,
        title={Classification and Enumeration of Solid-Solid Phase Transition Mechanisms},
        author={Wang, Fang-Cheng and Ye, Qi-Jun and Zhu, Yu-Cheng and Li, Xin-Zheng},
        journal={arXiv:2506.05105},
        year={2025},
        doi={10.48550/arXiv.2506.05105}
    }
    ```

## Usage

`crystmatch` can be used both as a command-line tool and as a Python module; see the [documentation](https://fangcheng-wang.github.io/crystmatch/) for details.