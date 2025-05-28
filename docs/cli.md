# Command-line interface

## Flowchart

``` mermaid
flowchart TB
  mode{how to **get CSMs**?} --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--enumeration</code>"| modeE[enumerate a list of<br>**representative** CSMs];
  modeE --> full;
  mode --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--read</code>"| modeR["read a list of CSMs<br>from an **NPZ** file<br>"];
  modeR --> full;
  mode --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--direct</code>"| modeD["read a single CSM from<br>a pair of **POSCAR** files"];
  modeD --> full["if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--full</code>, enumerate<br>**all** CSMs for each SLM"];
  full --> result["save NPZ file and visualize<br>by <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--plot --interactive</code>"];
  result --> export["save CSMs in other formats<br>by <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--poscar --xdatcar</code>"]
```

## CSMCAR

CSMCAR is the configuration file of `crystmatch`. It is a plain text file specified in the `--config` or `-c` option:

```
crystmatch --enumeration Rocksalt.poscar CsCl.poscar 2 0.3 --config CSMCAR
```

CSMCAR may contain the elastic tensors, the atomic weights used to defined the shuffle distance, and/or the orientation relationship used to benchmark CSMs.

### Elastic tensors

If you want to use the strain energy density $w$, instead of the root-mean-square strain $\bar\varepsilon$, as the pruning criterion, your CSMCAR should contain elastic tensors of the initial and final structures. They should be given in as $6\times 6$ symmetric matrices according to the Voigt notation, as:

```
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

The axis name must be given in each row. The above CSMCAR is in the order of `XX`, `YY`, `ZZ`, `XY`, `YZ`, `ZX`, which is the same as the [OUTCAR of VASP](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences). However, the elastic tensors on [Materials Project](https://docs.materialsproject.org/methodology/materials-methodology/elasticity) are in the order of `XX`, `YY`, `ZZ`, `YZ`, `ZX`, `XY`, in which case you should write the axis name in the same order.

!!! warning "Caveat"
    Once you use `--config` or `-c` and provide the elastic tensors, the `--enumeration` mode will interpret its last parameter `MAX_W` as the maximum strain energy density, not the maximum root-mean-square strain (RMSS). The unit of `MAX_W` is the same as the unit of the elastic tensors in the CSMCAR file.

    If you want to use the RMSS as the pruning criterion, and then calculate the strain energy density $w$ for each CSM, you can use `--enumeration` without providing the elastic tensors, and then use `--read` and `--config` to calculate $w$ for each CSM in the NPZ file.

### Atomic weights



### Orientation relationship

## All options

The `crystmatch` command-line interface is already introduced in [Tutorial](https://fangcheng-wang.github.io/crystmatch/). Here, we provide a table of all available options:

| Option | Description | Mode |
| --- | --- | --- |
| `-h`, `--help` | Show help message and exit. | *none* |
| `-V`, `-v`, `--version` | Show program's version number and exit. | *none* |
| `-E MAX_MU MAX_RMSS`<br>`--enumeration MAX_MU MAX_RMSS` | Use 'enumeration' mode, with `MAX_MU` and `MAX_RMSS` as the multiplicity and RMSS upper bounds. | 'enumeration' |
| `-R CSM_LIST`<br>`--read CSM_LIST` | Use 'read' mode, with CSMs loaded from an existing NPZ file `CSM_LIST`. | 'read' |
| `-S`<br>`--single` | Use 'single-CSM' mode, with the CSM uniquely determined by `-I` and `-F`. | 'single-CSM' |
| `-I POSCAR_I`<br>`--initial POSCAR_I` | POSCAR file of the initial crystal structure. | 'enumeration', 'single-CSM' |
| `-F POSCAR_F`<br>`--final POSCAR_F` | POSCAR file of the final crystal structure. | 'enumeration', 'single-CSM' |
| `-e INDEX1 [INDEX2 ...]`<br>`--export INDEX1 [INDEX2 ...]` | Export CSMs from NPZ file with the given indices. | 'read' |
| `-i [IMAGES]`<br>`--interpolate [IMAGES]` | Create XDATCAR files when `-e` or `--export` is used. `IMAGES` is the number of images to be added; default is 10. | 'read' |
| `-t TOL`<br>`--tolerance TOL` | Tolerance for determining crystal symmetry; default is 1e-3. | 'enumeration', 'single-CSM' |
| `-o vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz`<br>`--orientation vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz` | Benchmark CSMs by their deviation angles from the orientation relationship: "\(\mathbf{v}_\text{i}\parallel\mathbf{v}_\text{f}\) and \(\mathbf{w}_\text{i}\parallel\mathbf{w}_\text{f}\)".<br>Note: arguments must be given in **Cartesian** coordinates. | *all* |
| `-c`, `--csv` | Create CSV file. | *all* |
| `-p`, `--plot` | Create scatter plot. | *all* |
| `-u`, `--uspfix` | Use USP-fixed manner in `--orientation`. | *all* |