# Command-line interface

``` mermaid
flowchart TB
  A{mode?} --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--enumeration</code>"| B[enumerate a list of<br>representative CSMs];
  B --> E[save CSMs to NPZ file];
  E --> F;
  A --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--read</code>"| C["read CSMs from NPZ file<br>(after <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--enumeration</code>)"];
  C --> G["export CSM to POSCAR<br>files if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--export</code>"];
  G --> F;
  A --> |"<code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--single</code>"| D["read a single CSM from<br><code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>-I</code> and <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>-F</code> POSCAR files"];
  D --> F["analyze orientation<br>if <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--orientation</code>"];
  F --> I["generate reports according<br>to <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--csv</code> and <code style='color: #dd0099 !important; background-color: var(--md-code-bg-color); border-radius: .1rem; font-size: .85em; padding: 0 .2941176471em; word-break: break-word;'>--plot</code>"];
```

The `crystmatch` command-line interface is already introduced in [Tutorial](https://fangcheng-wang.github.io/crystmatch/). Here, we provide a table of all available options:

| Option | Description | Mode |
| --- | --- | --- |
| `-h`, `--help` | Show help message and exit. | *none* |
| `-v`, `--version` | Show program's version number and exit. | *none* |
| `-E MAX_MU MAX_RMSS`<br>`--enumeration MAX_MU MAX_RMSS` | Use 'enumeration' mode, with `MAX_MU` and `MAX_RMSS` as the multiplicity and RMSS upper bounds. | 'enumeration' |
| `-R CSM_LIST`<br>`--read CSM_LIST` | Use 'read' mode, with CSMs loaded from an existing NPZ file `CSM_LIST`. | 'read' |
| `-S`<br>`--single` | Use 'single-CSM' mode, with the CSM uniquely determined by `-I` and `-F`. | 'single-CSM' |
| `-I POSCAR_I`<br>`--initial POSCAR_I` | POSCAR file of the initial crystal structure. | 'enumeration', 'single-CSM' |
| `-F POSCAR_F`<br>`--final POSCAR_F` | POSCAR file of the final crystal structure. | 'enumeration', 'single-CSM' |
| `-a`<br>`--accurate` | Use more accurate algorithm for RMSD minimization, taking about 4x longer time. | 'enumeration' |
| `-e index1 [index2 ...]`<br>`--export index1 [index2 ...]` | Export CSMs from NPZ file with the given indices. | 'read' |
| `-t TOL`<br>`--tolerance TOL` | Tolerance for determining crystal symmetry; default is 1e-3. | 'enumeration', 'single-CSM' |
| `-o vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz`<br>`--orientation vix viy viz vfx vfy vfz wix wiy wiz wfx wfy wfz` | Benchmark CSMs by their deviation angles from the orientation relationship: "\(\mathbf{v}_\text{i}\parallel\mathbf{v}_\text{f}\) and \(\mathbf{w}_\text{i}\parallel\mathbf{w}_\text{f}\)".<br>Note: arguments must be given in **Cartesian** coordinates. | *all* |
| `-c`, `--csv` | Create CSV file. | *all* |
| `-p`, `--plot` | Create scatter plot. | *all* |
| `-u`, `--uspfix` | Use USP-fixed manner in `--orientation`. | *all* |