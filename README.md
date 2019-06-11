# doebase

Basic Optimal DoE (Design of Experiments) routines for metabolic pathway design.

Modules:

* `OptDes`: main module of optimal DoE routines.
* `doebase`: routines to interface the current SynBioChem DoE input sheet template with the `OptDes` module.

The DoE specifications can be provided in tabular form (`csv` or `xlsx`) with the following columns:

 * `DoE position`: Position in the genetic construct.
 * `DoE permutation class`: When labeled as `1`, it means that the genetic part can be rearranged in any other labeled position. This option is typically used for gene coding regions. It is assumed that each position forms a unit, i.e., rearrangements always take place by selecting one coding region from each labeled position.
 * `DoE designation`: type of genetic part: origin | resistance | promoter | gene.
 * `Part number`: identifier linking to the parts repository (`SynBioHub`, 'ICE`, etc.).




Basic call from within the library (see `evaldes` for details):

```python
factors, fnames, diagnostics = makeDoeOptDes(fact, size, seed=None, starts=1040, makeFullFactorial=False, RMSE=1, alpha=0.05, verbose=False, random=False)
```

## Input parameters:

* `fact`: a dictionary of sortable keys (the positions in the construct) containing `doebase.spec` objects with the following attributes:

  *  `positional`: `float`

      `1.0` or `None` depending if the genetic part can be rearranged

  *  `component`: `str`

      origin | resistance | promoter | gene

  *  `levels`: `list`

      levels of the genetic part (see Note)

### Note

  *  Origin levels (plasmid copy numbers)                   

        [`pl1`, `pl2`, ... ]

  *  Resistance levels

        [`res1`, ... ]

  *  Promoter levels (and blanks)

        [`prom1`, `prom2`, ..., '-', '-', ... ]

  *  Gene levels (gene variants)

        [`g1_1`, `g1_2`, ... ]

* `size`: size of the library.
* `seed`: random seed.
* `starts`: number of the starts for the DoE algorithm.
* `makeFullFactorial`: make a full factorial rather than a library of the given `size`.
* `RMSE`, `alpha`: parameters of the DoE algorithm.
* `verbose`
* `random`: make a random rather than an optimal design.
