module BEASTCUDAExt

using CUDA
using CUDA.Adapt
using CUDA.CUSPARSE

using BEAST
import BEAST: assemble!, Threading, Operator, Space, IntegralOperator
import BEAST: _integrands, _integrands_gen, Integrand, pulledback_integrand
import BEAST: LagrangeRefSpace, RTRefSpace, GWPDivRefSpace, GWPCurlRefSpace
using BEAST.CompScienceMeshes
using BEAST.SauterSchwabQuadrature
using BEAST.StaticArrays
using BEAST.SparseArrays
using BEAST.LinearAlgebra
using BEAST.ProgressMeter

Adapt.@adapt_structure CommonVertex
Adapt.@adapt_structure CommonEdge
Adapt.@adapt_structure CommonFace

# Adapt.@adapt_structure GWPDivRefSpace

function Adapt.adapt_structure(to, obj::GWPDivRefSpace{T,Degree}) where {T,Degree}
    GWPDivRefSpace{T,Degree}() #Adapt.adapt_structure(to, obj.storage))
end

function Adapt.adapt_structure(to, obj::GWPCurlRefSpace{T,Degree}) where {T,Degree}
    GWPCurlRefSpace{T,Degree}() #Adapt.adapt_structure(to, obj.storage))
end

# include("BEASTCUDAExt/gpu_refspace.jl")

include("BEASTCUDAExt/tiling.jl")

include("BEASTCUDAExt/cpu_assemble.jl")

include("BEASTCUDAExt/gpu_utils.jl")
include("BEASTCUDAExt/gpu_basis.jl")
include("BEASTCUDAExt/gpu_integrals.jl")
include("BEASTCUDAExt/gpu_assemble_integralop.jl")

end