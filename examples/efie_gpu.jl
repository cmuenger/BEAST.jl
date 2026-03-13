using CompScienceMeshes
using BEAST

using CUDA

# Γ = meshcuboid(1.0,1.0,1.0,1.0)
Γ = meshsphere(1.0,0.1;generator=:gmsh)

X = raviartthomas(Γ)
#Y = buffachristiansen(Γ)

@show numcells(Γ)
@show numfunctions(X)

κ, η = 1.0, 1.0
T = Maxwell3D.singlelayer(wavenumber=κ)

qstrat = BEAST.DoubleNumSauterQstrat(4, 4, 6, 6, 6, 6)


@time Th_gpu = assemble(T,X,X;threading=:gpu,gpu_tiling=(1,1),quadstrat=qstrat)

@time Th_cpu = assemble(T,X,X;threading=:cellcoloring,quadstrat=qstrat)
@show numfunctions(X)
@show Threads.nthreads()
@show sum(Th_gpu-Th_cpu)