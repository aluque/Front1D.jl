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

using Logging: global_logger, with_logger
using TerminalLoggers: TerminalLogger


include("model.jl")
include("naidis.jl")
include("plot.jl")
include("prettyprint.jl")
include("main.jl")

end # module Front1D
