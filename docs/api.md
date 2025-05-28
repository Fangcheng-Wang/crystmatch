# Python API

!!! warning "Under development"
    **The API reference for `crystmatch` v1.2.0 is still under development. Please do not rely on it yet.**


## Data structures

| Type | Structure | Description |
| --- | --- | --- |
| `cryst` | `(lattice, species, positions)` | A crystal structure described by a primitive cell `lattice`, atomic species `species`, and fractional coordinates `positions`. Cell vectors are **rows** of `lattice`. |
| `slm` | `(hA, hB, q)` | A [sublattice match](https://arxiv.org/abs/2305.05278) (SLM) described by an integer-matrix triplet. |

## Functions

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