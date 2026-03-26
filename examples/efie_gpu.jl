using CompScienceMeshes
using BEAST

using LinearAlgebra
# using SparseArrays
# using Profile

using CUDA


Γ = meshcuboid(1.0,1.0,1.0,1.0)


#Γ = meshsphere(1.0,0.1;generator=:gmsh)


X = raviartthomas(Γ)
Y = buffachristiansen(Γ)
Z = lagrangec0d1(Γ)

@show numcells(Γ)
@show numcells(geometry(Y))

@show numfunctions(X)
@show numfunctions(Y)
@show numfunctions(Z)

κ, η = 1.0, 1.0
T = Maxwell3D.singlelayer(wavenumber=κ)

V = Helmholtz3D.singlelayer(wavenumber=κ)


qstrat = BEAST.DoubleNumSauterQstrat(4, 4, 6, 6, 6, 6)


CUDAExt = Base.get_extension(BEAST, :BEASTCUDAExt)

gpu_tstrat = CUDAExt.TilingStrategy(CUDAExt.EqualTiling(1), CUDAExt.EqualTiling(1))

CUDA.@time Th_gpu = assemble(T,Y,Y;threading=:gpu,tilingstrat=gpu_tstrat,quadstrat=qstrat)

cpu_tstrat = CUDAExt.TilingStrategy(CUDAExt.WorksizeTiling(128), CUDAExt.WorksizeTiling(128))

@time Th_cpu2 = assemble(T,Y,Y;threading=:cellsplitting,tilingstrat=cpu_tstrat,quadstrat=qstrat)

@time Th_cpu3 = assemble(T,Y,Y;threading=:dofsplitting,quadstrat=qstrat)

@time Th_cpu = assemble(T,Y,Y;threading=:cellcoloring,quadstrat=qstrat)

@show numfunctions(X)
@show Threads.nthreads()
@show eps(real(eltype(Th_cpu))) maximum(norm.(Th_cpu2-Th_cpu)) 



