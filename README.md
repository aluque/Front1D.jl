# One-dimensional model of a streamer front

This code implements the model described in the publication "Collective dynamics of a dense streamer front", by M. Gomes, F.J. Gordillo-Vázquez and A. Luque.

## Install
The code is written in julia and distributed as a julia package.  To install it:

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/aluque/Front1D.jl")
```

## Run
To run the code and reproduce the figure presented in the paper:

```julia
julia> using Front1D
julia> Front1D.main()
```

You can inspect the code to make changes.  The default simulation parameters are in `src/model.jl`.

## Differences with the paper
Note that the variables $\bar{q}$ and $\bar{n}_e$ in the code are defined differently than in the paper. In the code they are *not* normalized by the number of streamers $N$. This gives slightly different electrostatic field expressions.






