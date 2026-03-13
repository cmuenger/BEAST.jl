module BEASTCUDAExt

using CUDA
using CUDA.Adapt

using BEAST
import BEAST: assemble!, Threading, Operator, Space, IntegralOperator
using BEAST.CompScienceMeshes
using BEAST.SauterSchwabQuadrature
using BEAST.StaticArrays
using BEAST.LinearAlgebra


Adapt.@adapt_structure CommonVertex
Adapt.@adapt_structure CommonEdge
Adapt.@adapt_structure CommonFace


include("gpu_utils.jl")
include("gpu_assemble_integralop.jl")

end