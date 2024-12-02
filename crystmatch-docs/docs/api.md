
# Python API

## Data Structures

| Type | Structure | Description |
| --- | --- | --- |
| `cryst` | `(lattice, species, positions)` | A crystal structure described by a primitive cell `lattice`, atomic species `species`, and fractional coordinates `positions`. Cell vectors are **rows** of `lattice`. |
| `slm` | `(hA, hB, q)` | A sublattice match described by an integer-matrix triplet; see our [paper](https://arxiv.org/abs/2305.05278). |

## Available Functions

::: crystmatch.enumeration
    handler: python
    options:
        docstring_style: numpy
        show_root_heading: true
        show_source: false

::: crystmatch.analysis
    handler: python
    options:
        docstring_style: numpy
        show_root_heading: true
        show_source: false

::: crystmatch.utilities
    handler: python
    options:
        docstring_style: numpy
        show_root_heading: true
        show_source: false