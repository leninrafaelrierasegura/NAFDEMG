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
