
# Python API

## Data Structures

| Type | Structure | Description |
| --- | --- | --- |
| `Cryst` | 3-tuple, as `(lattice, species, positions)` | A crystal structure described by a cell matrix `lattice`, a list of atomic species `species`, and a list of fractional coordinates `positions`. |
| `SLM` | 3-tuple of (3, 3) arrays, as `(hA, hB, q)` | A sublattice match described by an IMT; see our [paper](https://arxiv.org/abs/2305.05278). |
| `CSM` | 3-tuple, as `(slm, permutation, translations)` | A crystal-structure match described by a sublattice match `x`, an integer vector `permutation`, and an integer matrix `translations`. |

## Available Functions

::: crystmatch
    handler: python
    options:
        docstring_style: numpy
        show_root_heading: false
        show_source: false