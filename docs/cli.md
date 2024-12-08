# Command-line interface


The `crystmatch` command-line interface is already introduced in [Tutorial](https://fangcheng-wang.github.io/crystmatch/). Here, we provide a table of all available options, and a flowchart to show how they interact with each other.

## Table of options

| Option | Description | Mode |
| --- | --- | --- |
| `-h`, `--help` | Show help message and exit. | *none* |
| `-v`, `--version` | Show program's version number and exit. | *none* |
| `-E MAX_MU MAX_RMSS`<br>`--enumeration MAX_MU MAX_RMSS` | Use 'enumeration' mode, with `MAX_MU` and `MAX_RMSS` as the multiplicity and RMSS upper bounds. | 'enumeration' |
| `-R CSM_LIST`<br>`--read CSM_LIST` | Use 'read' mode and read CSMs from NPZ file. | 'read' |
| `-S`<br>`--single` | Use 'single-CSM' mode and analyze a single CSM defined by the POSCARs specified by `-I` and `-F`. | 'single-CSM' |
| `-I POSCAR_I`<br>`--initial POSCAR_I` | POSCAR file of the initial crystal structure. | 'enumeration', 'single-CSM' |
| `-F POSCAR_F`<br>`--final POSCAR_F` | POSCAR file of the final crystal structure. | 'enumeration', 'single-CSM' |
| `-t TOL`<br>`--tolerance TOL` | Tolerance for determining crystal symmetry; default is 1e-5. | 'enumeration', 'single-CSM' |
| `-e index1 [index2 ...]`<br>`--export index1 [index2 ...]` | Export CSMs from NPZ file with the given indices. | 'read' |
| `-l`, `--list` | Generate a list of CSMs in 'single-CSM' mode. | 'single-CSM' |
| `-o vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz`<br>`--orientation vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz` | Benchmark CSMs by their deviation angles from the orientation relationship: "\(\mathbf{v}_\text{i}\parallel\mathbf{v}_\text{f}\) and \(\mathbf{w}_\text{i}\parallel\mathbf{w}_\text{f}\)".<br>Note: arguments must be given in Cartesian coordinates. | *all* |
| `-c`, `--csv` | Create CSV file. | *all* |
| `-p`, `--plot` | Create scatter plot. | *all* |
| `-u`, `--uspfix` | Use USP-fixed manner in `--orientation`. | *all* |

## Flowchart

``` mermaid
flowchart TD
  B{mode?} --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--enumeration</code>"| C([enumerate CSMs]);
  B --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--single</code>"| E(["read a single CSM<br>from <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>-I</code> and <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>-F</code>"]);
  B --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--read</code>"| D([read CSMs from NPZ]);
  C --> F([save CSMs to NPZ]);
  E --> |"if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--list</code>"| C;
  E --> |"if not <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--list</code>"| H(["analyze orientation<br>if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--orientation</code>"]);
  D --> I(["export CSM to POSCAR<br>if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--export</code>"]);
  I --> H;
  F --> H;
  H --> J(["generate reports according to <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--csv</code> and <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--plot</code>"]);
```

