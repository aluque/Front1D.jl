module Front1D

using Constants: co
using RecursiveArrayTools
using SpecialFunctions
using LoopVectorization
using Polyester
using LinearAlgebra
using SparseArrays
using OrdinaryDiffEq
using QuadGK
using PyPlot
using LaTeXStrings

include("model.jl")
include("naidis.jl")
include("plot.jl")
include("prettyprint.jl")
include("main.jl")

end # module Front1D
