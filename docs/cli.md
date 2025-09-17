# Command-line interface

## Flowchart

If you are not using `--interpolate` or `--polish`, the flowchart is as follows:

``` mermaid
flowchart TB
  config["if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--extra</code>,<br>read **CSMCAR** file"] --> mode;
  mode{how to<br>**get CSMs**} --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--enumerate</code>"| modeE[enumerate a list of<br>**representative** CSMs];
  modeE --> all;
  mode --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--read</code>"| modeR["read a list of CSMs<br>from an **NPZ** file<br>"];
  modeR --> all;
  mode --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--direct</code>"| modeD["read a single CSM from<br>a pair of **POSCAR** files"];
  modeD --> all["if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--all</code>, enumerate<br>**all** CSMs for each SLM"];
  all --> export["**export** individual CSMs as<br><code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--poscar --xdatcar</code>"];
  export --> visualize["**visualize** CSMs using<br><code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; font-family: monospace; word-break: break-word;'>--scatter --interact</code>"]
```

## CSMCAR

CSMCAR is a plain text file that provides extra information for generating and analyzing CSMs. `crystmatch` will **not** read a CSMCAR file unless it is specified in the `--extra` option, as:

```
crystmatch --enumerate graphite.poscar diamond.poscar 2 0.3 --extra CSMCAR.txt
```

CSMCAR may contain at most three fields:

- The elastic tensors
- The atomic weights used to defined the shuffle distance
- The orientation relationship used to benchmark CSMs

### Elastic tensors

If you want to use the strain energy density $w$, instead of the root-mean-square strain $\bar\varepsilon$, as the pruning criterion, the CSMCAR file should contain elastic tensors of the initial and final structures. They must be specified as $6\times 6$ symmetric matrices ([Voigt notation](https://en.wikipedia.org/wiki/Voigt_notation)) as:

```python
# ./CSMCAR.txt

initial elastic tensor:
    XX    268.6236     46.1585     47.5328      0.6739      0.9388     -1.8586
    YY     46.1585    268.5367     48.5034      0.1023      1.8528     -3.2724
    ZZ     47.5328     48.5034    268.7187     -0.0438      0.8634     -1.6993
    XY      0.6739      0.1023     -0.0438     44.0620     -1.7433      1.4076
    YZ      0.9388      1.8528      0.8634     -1.7433     42.5860     -0.5327
    ZX     -1.8586     -3.2724     -1.6993      1.4076     -0.5327     42.2492

final elastic tensor:
    XX    262.9004    128.7075    102.3784    -22.5882      5.6997      7.9351
    YY    128.7075    262.8757    102.3694    -21.7608      8.0470      5.6492
    ZZ    102.3784    102.3694    295.9447     32.4762    -19.2804    -19.2551
    XY    -22.5882    -21.7608     32.4762    138.7953      5.1622      5.4153
    YZ      5.6997      8.0470    -19.2804      5.1622    104.6801     34.9596
    ZX      7.9351      5.6492    -19.2551      5.4153     34.9596    104.6258
```

The axis name must be given in each row. The above CSMCAR is in the order of `XX`, `YY`, `ZZ`, `XY`, `YZ`, `ZX`, which is the same as [VASP](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences). However, the elastic tensors provided by [Materials Project](https://docs.materialsproject.org/methodology/materials-methodology/elasticity) are in the order of `XX`, `YY`, `ZZ`, `YZ`, `ZX`, `XY`, in which case you should write the axis name in the same order.

!!! note
    Once you use `--extra` or `-e` and provide the elastic tensors, the `--enumerate` mode will interpret its last parameter `MAX_STRAIN` as the maximum strain energy density, not the maximum root-mean-square strain (RMSS). The unit of `MAX_STRAIN` is the same as the unit of the elastic tensors in the CSMCAR file.

    If you want to use the RMSS as the pruning criterion, and then calculate the strain energy density $w$ for each CSM, please use `crystmatch` as a [Python module](https://fangcheng-wang.github.io/crystmatch/api/).

### Atomic weights

If the initial and final structures involve more than one atomic species, you can specifie their atomic weights in the CSMCAR file as:

```python
# ./CSMCAR.txt

atomic weights:
    Cs      Cl
    132.9   35.45
```

`crystmatch` will normalize these weights and use them to define the shuffle distance

$$
d(\mathcal{J})=\min_{\mathbf{t}_0 \in \mathbb{R}^{3}}\sqrt{\sum_{i=1}^{\tilde Z}\theta_i\left|\tilde{\mathcal{J}}(\tilde{\mathbf{a}}_i)+\mathbf{t}_0-\tilde{\mathbf{a}}_i\right|^2},
$$

where $\mathcal{J}\colon\mathcal{A}\to\mathcal{B}$ is a CSM, $\tilde{\mathcal{J}}\colon\tilde{\mathcal{A}}\to\tilde{\mathcal{B}}$ its standard shuffle, $\tilde Z$ its period, $\tilde{\mathbf{a}}_i\in\tilde{\mathcal{A}}$ the position of the $i$-th atom, and $\theta_i$ the normalized weight of the $i$-th atom. See [this paper](https://arxiv.org/abs/2506.05105) for more details.

### Orientation relationship

To benchmark CSMs with a given orientation relationship, two parameters must be specified:

1. The orientation relationship, provided in the CSMCAR file as:

    ```python
    # ./CSMCAR.txt

    orientation relationship:
        1 0 0 || -1 0 0
        0 1 0 || 0 -1 0
    ```

    Here, exactly two rows after `orientation relationship` will be read. Each row specifies a parallelism between two **crystal directions** of the initial and final structures. These crystal directions are denoted as [Miller indices](https://en.wikipedia.org/wiki/Miller_index)

    $$
    [hk\ell]=h\mathbf{a}+k\mathbf{b}+\ell\mathbf{c},
    $$
    
    where $\mathbf{a}, \mathbf{b}, \mathbf{c}$ are the lattice vectors of the **conventional cell**. They are determined by **[Spglib](https://spglib.readthedocs.io/en/stable/definition.html#def-idealize-cell)** and depend only on the crystal systems of the initial and final structures.

    !!! warning "Caveat"
        **The conventional cell is not the primitive cell if the [Bravais lattice](https://en.wikipedia.org/wiki/Bravais_lattice) is base-centered, face-centered, or body-centered.** We recommend using
        
        ```
        crystmatch --read CSMLIST.npz 0 --interact
        ```

        to check if the conventional $\mathbf{a}, \mathbf{b}, \mathbf{c}$ is the same as you expect.

2. The assumption used to determine the orientation of the final structure, must be one of the following:

    - `'norot'`: The final lattice is deformed from the initial one **without rotation**.

    - `'uspfixed'`: One of the two **uniformly-scaled planes (USPs)** is fixed.

    For more details, see the supplemental material of [this paper](https://arxiv.org/abs/2305.05278).

If you have a `CSMLIST.npz` file and a `CSMCAR.txt` file containing an orientation relationship, you can run:

```
crystmatch --read CSMLIST.npz --extra CSMCAR.txt --orientation 'norot'
```
or

```
crystmatch --read CSMLIST.npz --extra CSMCAR.txt --orientation 'uspfixed'
```

## All options

```bash
crystmatch [--help] [--version]
           [--extra CSMCAR] [--tol TOL]
           [--enumerate POSCAR_I POSCAR_F MAX_MU MAX_STRAIN]
           [--read CSMLIST [IND1 IND2 ...]]
           [--direct POSCAR_I POSCAR_F] [--fix-integer]
           [--interpolate POSCAR_I POSCAR_F N_IMAGES]
           [--all MAX_D]
           [--l ELL]
           [--interact [SIZE]]
           [--poscar [ASSUM]] [--xdatcar [ASSUM]] [--mediumcell]
           [--orientation ASSUM]
           [--csv] [--scatter]
```

If you are not familiar with `crystmatch`, we recommend you first read the [tutorial](https://fangcheng-wang.github.io/crystmatch/) and then come back to this page for a detailed explanation of the options.

- `-h`, `--help`

    Show help message and exit.

- `-V`, `-v`, `--version`

    Show program's version number and exit.

- `-E POSCAR_I POSCAR_F MAX_MU MAX_STRAIN`, `--enumerate POSCAR_I POSCAR_F MAX_MU MAX_STRAIN`

    Enumerate representative CSMs (those with the minimal shuffle distance among all CSMs that share the same deformation gradient) between the initial and final structures defined in `POSCAR_I` and `POSCAR_F`, with `MAX_MU` and `MAX_STRAIN` as upper bounds for the multiplicity and RMSS. If elastic tensors are provided, `MAX_STRAIN` is instead interpreted as the maximum strain energy density (with the same unit as in the [CSMCAR](#elastic-tensors) file).

- `-R CSMLIST [IND1 IND2 ...]`, `--read CSMLIST [IND1 IND2 ...]`

    Read CSMs from an existing NPZ file `CSMLIST`. If the indices `IND1 IND2 ...` are provided, only load the CSMs with the given indices.

- `-D POSCAR_I POSCAR_F`, `--direct POSCAR_I POSCAR_F`

    Directly read a single CSM from a pair of POSCAR files `POSCAR_I` and `POSCAR_F`. The CSM is determined such that each cell vector or atomic position in `POSCAR_I` is mapped to the corresponding one (with the same row index) in `POSCAR_F`. **Make sure that the cell vectors and atomic positions in `POSCAR_I` and `POSCAR_F` are ordered as you expect.**

    !!! warning "Caveat"
        When writing the POSCAR, some softwares (e.g. [VASP](https://www.vasp.at/)) may add or subtract integers (cell vectors) to the fractional (cartesian) coordinates of the atoms to make them within the cell, thus altering the CSM. To avoid this issue, `crystmatch` will add and subtract integers to the fractional coordinates to restore the original CSM. If you do not want to do this, use the `--fix-integer` option described below.

- `-I POSCAR_I POSCAR_F N_IMAGES`, `--interpolate POSCAR_I POSCAR_F N_IMAGES`
    
    Create `IMAGES` interpolated structures between `POSCAR_I` and `POSCAR_F` when using `--direct`, similar to `nebmake.pl` provided by [VTST scripts](https://theory.cm.utexas.edu/vtsttools/scripts.html). The interpolated structures can be used for subsequent NEB calculations.

    !!! danger "Important"
        The original `nebmake.pl` does not preserve the CSM between `POSCAR_I` and `POSCAR_F`. Specifically, it will first reduce each fractional coordinate to `[0,1)` and then interpolate the reduced coordinates. This may lead to unexpected CSMs and thus incorrect transition paths as well as energy barriers. Therefore, we always recommend generating interpolated structures using
        
        ```
        crystmatch -I POSCAR_I POSCAR_F N_IMAGES
        ```

- `-P N_IMAGES`, `--polish N_IMAGES`

    Create and evenly distribute `N_IMAGES` images along the transition path stored in `00`, `01`, `02`, ..., similar to `--interpolate` but not limited to the two endpoints. The new images will be stored in `IMAGES/00`, `IMAGES/01`, `IMAGES/02`, ..., and can be used for subsequent NEB calculations.

- `-f`, `--fix-integer`

    Do not add or subtract integers to the fractional coordinates when using `--direct` or `--interpolate`.

    !!! warning "Caveat"
        If your input POSCAR files are generated by `crystmatch`, you should always use this option, or the CSM may be altered.

- `-l ELL`

    The norm used to quantify the shuffle distance; see [this paper](https://arxiv.org/abs/2506.05105) for more details. Default is `2.0` (RMSD).

- `-e CSMCAR`, `--extra CSMCAR`

    Load extra parameters from a [CSMCAR](#csmcar) file, including the [elastic tensors](#elastic-tensors), the [atomic weights](#atomic-weights), and the [orientation relationship](#orientation-relationship) used to benchmark CSMs.

- `-a MAX_D`, `--all MAX_D`

    Additionally enumerate all CSMs (instead of only the representative or read ones) for each SLM, with `MAX_D` as the upper bound for the shuffle distance.

- `-t TOL`, `--tol TOL`

    Tolerance in angstroms for detecting symmetry. Default is `1e-3`.

- `-i [SIZE]`, `--interact [SIZE]`

    Interactively visualize each CSM using a 3D plot. `SIZE` specifies the radius of the cluster to display, which we recommend to be smaller than `2.0`. Use `SIZE = 0` to display the cell contents only. Default is `1.2`.

    !!! warning "Caveat"
        The interactive plot may not work properly in some environments, such as Windows Subsystem for Linux (WSL) or remote servers.

- `-o ASSUM`, `--orientation ASSUM`

    Calculate the deviation angle from the given [orientation relationship](#orientation-relationship) for each CSM. Only works when the orientation relationship is provided by `--extra`. The [orientation of the final structure](#orientation-relationship) is determined by `ASSUM`, which must be specified as `'norot'` or `'uspfixed'`.

- `-p [ASSUM]`, `--poscar [ASSUM]`

    Export each enumerated/read CSM as a pair of POSCAR files. `ASSUM` must be `'norot'` or `'uspfixed'`. Default is `'norot'`.

    !!! note
        If `ASSUM = 'uspfixed'`, two final structures (instead of one) are created since there are two USPs.

- `-x [ASSUM]`, `--xdatcar [ASSUM]`

    Export each CSM as an XDATCAR file with 50 frames, similar to `--poscar`.

- `-m`, `--mediumcell`

    Output the average supercell for `POSCAR_I` and `POSCAR_F` when using `--poscar 'norot'` and `--xdatcar 'norot'`. This option is useful when subsequent calculations require the initial and final structures to be described by the same unit cell, but should only be used when the deformation is negligible.
    
**These options are automatically enabled when using `--enumerate`, `--all` or `--orientation`:**

- `-c`, `--csv`

    Create CSV file containing the multiplicities, RMSS, shuffle distances, estimated strain energy densities (if elastic tensors are provided) and deviation angles (if `--orientation` is used) of the enumerated/read CSMs. 

- `-s`, `--scatter`

    Create scatter plot in PDF format containing all enumerated/read CSMs.
