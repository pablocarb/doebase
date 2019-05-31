# doebase

Basic Optimal DoE (Design of Experiments) routines for metabolic pathway design.

Modules:

* `OptDes`: main module of optimal DoE routines.
* `doebase`: routines to interface the current SynBioChem DoE input sheet template with the `OptDes` module.

Basic call (see `evaldes` for details):

```python
factors, fnames, diagnostics = makeDoeOptDes(fact, size=libsize, seed=seed, starts=starts, RMSE= RMSE, alpha=alpha, random=random )
```

*To do:*

* Specifications for the DoE sheet.
* Connexion to synbiohub (use pysbol library).
