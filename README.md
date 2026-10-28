# NAFDEMG

This repository contains all the necessary files to reproduce the results presented in the paper available on [arXiv](https://arxiv.org/abs/0000.00000).

For a detailed description of this repository, follow this [link](https://leninrafaelrierasegura.github.io/NAFDEMG/index.html).

## Direct solver with compressed reference errors

On `reference-error-speedup-diffusion`, the three convergence notebooks keep
main's FEM evolution, package coefficient tables, physical grids, reference
meshes, forcing loads, and endpoint error quadrature. Each solver step is
still advanced through shifted sparse solves. Fixed shifted systems are
factorized once and their factors reused.

Set these controls in the parameter chunk of each notebook:

- `diffusion = 1` selects the coefficient nu in `u_t + nu L^(alpha/2) u = f`.
  It multiplies the operator outside the fractional power, while the physical
  time step and forcing stay unchanged. Use any finite positive value, fixed
  along a convergence curve. Exact reference rates use the same coefficient.
- `error_method = "compressed"` groups reference samples by solver interval
  and evaluates the original reference-matrix norm through centered coefficient
  moments. `"reference"` retains the original block accumulator for comparisons.
  The derivation is in `fun/fractional_evolution.Rmd`; it requires no discrete
  eigenfunction assumptions or Gaussian quadrature.
- `workers = 2L` runs independent alpha/order/time-step cases in fork workers
  on Unix. Set `1L` for serial execution or on other platforms. Only the
  coordinator writes results. `error_n_blocks` controls coefficient-block
  memory and does not alter reference samples or weights.
- `run_id = NULL` creates a fresh private run. Set a saved ID to reload it
  with exactly matching scientific parameters. Computation chunks refuse to
  overwrite existing result files.

Research runs and figures go to
`~/Desktop/DOCS/NAFDEMG/diffusion_experiments/<study>/<run_id>/`, following
this workspace's private-output convention. `NAFDEMG_PRIVATE_DIR` can override
that private project root. No Slack setup/uploads or presentation copies are
performed by these notebooks. Existing publication outputs are untouched.

Run the numerical regression suite from the repository root:

```sh
Rscript tests/check_reference_speedup.R
Rscript tests/check_reference_notebooks.R
Rscript tests/benchmark_reference_speedup.R
```

Use the project's renv environment for the actual notebook smoke test, whose
Graph API depends on the locked development version of MetricGraph. The
regression suite independently checks a nonuniform FEM problem with a
nonidentity mass matrix and nonnested reference mesh. It covers reference
weights, intervals crossing block boundaries, diffusion scaling, cached
solves, tiny errors, independent-worker ordering, and run metadata. The
benchmark compares both accumulators on identical inputs and includes time
stepping in its timings; shared spatial compression setup is excluded.

The editable helpers are `.Rmd` sources. Their generated `.R` and `.hash`
files are kept synchronized for `source_if_updated()`.


## Isolated rational-order experiment

`conv_in_m_isolated.Rmd` measures the difference between exact-discrete and
rational FEM/backward Euler solutions, using the same mesh, loads, diffusion,
and time grid. The primary rational solution still uses the shifted sparse
solver. The new functions are documented in
`fun/isolated_rational_convergence.Rmd`.

The defaults retain T=0.001 and the original physical data and package tables.
They use a moderate base h=0.01 and target step 0.00001, then halve h, tau,
and both. A complete numerical generalized eigensystem is computed once per
mesh; exact reference histories are shared across rational orders. Dense
allocation and history-size guards protect against accidentally using the
fine meshes of the original notebook. This experiment needs no reference
space-time quadrature and no analytical discrete-mode identities.

The notebook independently checks generalized eigen-residuals and mass
orthogonality, compares sparse evolution with the scalar rational product in
the complete numerical FEM eigenbasis, runs an alpha=2 zero-error control,
and reports mesh/time sensitivity, resolution indicators, package beta
interpolation brackets, curve monotonicity, and finite-order fitted slopes.
A plateau or nonmonotone computed curve is retained and diagnosed. A fitted
slope is not required to equal the theoretical bound exponent.

Run `Rscript tests/check_isolated_m.R` in the project environment for the
regression checks. The new notebook's RStudio Knit button uses a custom
renderer that keeps the report private. From R, run:

```r
source("sourceifupdated.R")
source_if_updated("fun/reference_experiments")
source_if_updated("fun/isolated_rational_convergence")
settings <- isolated_m_parameters()
# settings$diffusion <- 2
render_isolated_m_notebook(parameters = settings, workers = 4L)
```

Private metadata, results, CSV diagnostics, plots, and the rendered notebook
are saved under `diffusion_experiments/m/<run_id>/`. Reload with
`render_isolated_m_notebook(parameters = settings, run_id = "<saved ID>",
execute = FALSE)`. Scientific metadata and package versions must match.
Reloading regenerates figures and the report without rerunning trajectories.
Complete result files are never overwritten. Existing `conv_in_m.Rmd` runs
are separate and are not overwritten by this experiment.


### Computed coefficients in the isolated study

Both sources are available in `conv_in_m_isolated.Rmd`. Change its YAML
`coefficient_source: package` to `coefficient_source: computed` and click Knit.
The package choice retains orders 1--4; the computed choice selects orders
1--8. The YAML `coefficient_file` points to
this branch's `matlab/computed_fractional_roots.csv` by default. An external
CSV, including the original NAP file, can still be passed explicitly.
From R, an explicit computed run is:

```r
settings <- isolated_m_parameters("computed")
render_isolated_m_notebook(parameters = settings, workers = 4L)
```

Computed coefficients are matched at exact beta values (alpha = 0.5, 0.9,
1.3 for the supplied artifact); there is no beta interpolation, extrapolation
in m, or package fallback. CSV schema, complete blocks, factor positivity,
root signs and interlacing are checked before expensive setup. The source,
file path, SHA-256 and generator provenance distinguish saved runs. The
verified artifact was generated by `tests/audit_rational_generation.m` in
NAFDEMG's `diffusion-parameter` branch; a verbatim copy of that generator
and the CSV/provenance JSON is included here. The audit explored orders through ten;
the production CSV retains orders 1--8. MATLAB is not called by the notebook.

The sparse supplied-root update uses direct polynomial multiplication and
stable residue evaluation; its sampled reconstruction and the independent
complete-eigenbasis calculation check the implementation. Higher-order errors
can fall below numerical resolution. All points remain plotted, while the
diagnosis reports resolution flags and both all-order and resolved-only fits.
Existing package runs and the exact discrete reference are unchanged.

The copied audit writes to a new `audit_output_dir` and does not overwrite
production coefficient files. It requires MATLAB and Chebfun only when
regenerating/auditing coefficients. The notebook uses the supplied CSV and
does not require MATLAB. Selecting validated orders 1--8 and preserving the
17-digit root/factor export are separate artifact-preparation steps; running
the audit alone does not replace the production CSV.

Computed-source regression checks: `Rscript tests/check_computed_isolated_m.R`.
The plot uses open circles for errors below the diagnostic resolution
indicator; these points are retained rather than silently omitted.
