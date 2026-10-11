# spin-system
Generalized XXZ \
code1: https://www.yu.edu/faculty-bios/santos/computer-codes


## Exercise map

| Path | Content |
| --- | --- |
| `code1.txt` | MATLAB spin basis, fixed-magnetization Ising/XXZ Hamiltonians and parity/rotation exercises |
| `Generalized XXZ model` | Extensionless MATLAB long-range XXZ/entanglement draft, with an unfinished local OTOC function |
| `AAStarkPlot/PlotAAIpr.jl` | Average AA inverse participation ratio versus potential strength |
| `AAStarkPlot/PlotMFScaling.jl` | AA IPR finite-size fit |
| `AAStarkPlot/PlotStark.jl` | Stark-chain eigenstate probability |
| `AAStarkPlot/Ploteigenvalues.jl` | AA-chain eigenstate probability |

The Julia scripts use `LinearAlgebra` from the standard library and the external
`Plots` package; `PlotMFScaling.jl` also uses `LsqFit`.
No original Julia/package versions are recorded. Create an environment from the
repository root with:

```sh
julia --project=. -e 'using Pkg; Pkg.add(["Plots", "LsqFit"])'
```

Then run a chosen script from `AAStarkPlot/` with that environment, for example
`julia --project=.. PlotStark.jl`. PNG outputs are written to the current
directory and can overwrite the existing example images. The IPR scaling
script uses dense diagonalization up to L=1597; inspect the sizes before a sweep.

For the text-format MATLAB exercises, open `code1.txt` in the MATLAB editor
and run individual sections, or copy it to a local `.m` file.
The separate extensionless XXZ draft is not a polished standalone OTOC program:
its local `otoc` function refers to variables that are not passed into it, and
its two identical J2 terms need model clarification before changing them.

## Corrected selections and small check

Julia's `eigen(H)` stores the eigenvectors in columns. The two eigenstate plot
scripts now choose `Ev[:, 50]` and `Ev[:, 100]` respectively.
See the [Julia linear-algebra documentation](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.eigen).
Existing PNGs have not been regenerated; their plotted states may differ from
a new run after this correction.

The periodic Ising exercise now adds its closing bond rather than overwriting
the open-chain bonds. In MATLAB, from this repository directory:

```matlab
run('tests/check_ising_exercise.m')
```

The check compares all basis-state diagonal energies in that exercise with an
independent periodic-spin sum. It does not validate the unfinished XXZ/OTOC
draft or reproduce the Julia plots.

