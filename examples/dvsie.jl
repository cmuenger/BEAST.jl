using CompScienceMeshes, BEAST
using LinearAlgebra
using Profile
using LiftedMaps
using BlockArrays
using StaticArrays



ntrc = X->ntrace(X,Γ)

T = tetmeshsphere(1.0,0.25)
X = nedelecd3d(T)
Γ = boundary(T)
Y = raviartthomas(Γ)

@show numfunctions(X)
@show numfunctions(Y)

Z = BEAST.buffachristiansen2(Γ)



module Material

    const ϵ_0 = 1.0
    const ϵ_b = 1.5
    const ϵ_r = 2.0
    
    const μ_0 = 1.0
    const μ_b = 3.0
    const μ_r = 2.0

    wavenumber_freespace(ω) = ω*√(ϵ_0*μ_0)
    impedance_freespace() = √(μ_0/ϵ_0)
    
    wavenumber_background(ω) = ω*√(ϵ_0*ϵ_b*μ_0*μ_b)
    impedance_background() = √(μ_0*μ_b/(ϵ_0*ϵ_b))

    wavenumber_inside(ω) = ω*√(ϵ_0*ϵ_b*ϵ_r*μ_0*μ_b*μ_r)
    impedance_inside() = √(μ_0*μ_b*μ_r/(ϵ_0*ϵ_b*ϵ_r))
end

ω = 0.3
#wavenumbers and impedance
κ,  η  = Material.wavenumber_freespace(ω),Material.impedance_freespace()
κ′, η′ = Material.wavenumber_background(ω),Material.impedance_background()
κ_in, η_in = Material.wavenumber_inside(ω),Material.impedance_inside()

#contrast currents
χ_E = x->(1.0-1.0/Material.ϵ_r)
χ_H = x->(1.0-1.0/Material.μ_r)

#=
function tau(x::SVector{U,T}) where {U,T}
    1.0-1.0/Material.ϵ_r
end
χ = tau
=#


N = NCross()
#Volume-Volume
I = Identity()
Le,Be,Ke = VIE.singlelayer(wavenumber=κ′, tau=χ_E),  VIE.boundary(wavenumber=κ′, tau=χ_E), VIE.doublelayer(wavenumber=κ′, tau=χ_E)
Lh,Bh,Kh = VIE.singlelayer(wavenumber=κ′, tau=χ_H),  VIE.boundary(wavenumber=κ′, tau=χ_H), VIE.doublelayer(wavenumber=κ′, tau=χ_H)

#Volume-Surface 
Lt,Bt,Kt = transpose(VSIE.singlelayer(wavenumber=κ′)), transpose(VSIE.boundary(wavenumber=κ′)), transpose(VSIE.doublelayer(wavenumber=κ′))

Lse,Bse,Kse = VSIE.singlelayer(wavenumber=κ′, tau=χ_E), VSIE.boundary(wavenumber=κ′, tau=χ_E), VSIE.doublelayer(wavenumber=κ′, tau=χ_E)
Lsh,Bsh,Ksh = VSIE.singlelayer(wavenumber=κ′, tau=χ_H), VSIE.boundary(wavenumber=κ′, tau=χ_H), VSIE.doublelayer(wavenumber=κ′, tau=χ_H)
#Surface-Surface
T  = Maxwell3D.singlelayer(wavenumber=κ)  #Outside
T′ = Maxwell3D.singlelayer(wavenumber=κ′) #Inside
K  = Maxwell3D.doublelayer(wavenumber=κ)  #Outside
K′ = Maxwell3D.doublelayer(wavenumber=κ′) #Inside

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*κ*η)*curl(E)

e, h = (n × E) × n, (n × H) × n

@hilbertspace D B j m
@hilbertspace k l o p
#=
β = 1/(ϵ_r*ϵ_b)
ν = 1/ϵ_b
θ = 1/(μ_b*μ_r)

ρ = 1/μ_b
α, α′ = 1/η, 1/η′
γ′ = im*η′/κ′
λ′ = im/(η′*κ′)
ζ′ = im*κ′/(ϵ_b*η′)
δ′ = im*κ′/ϵ_b

#= D-VSIE (D-VIE combined with PMCHWT) =#

eq = @varform (β*I[k,D]-ν*Le[k,D]-ν*Be[ntrc(k),D] + η′*Lt[k,j]-γ′*Bt[ntrc(k),j] -          Kt[k,m] +
                    -δ′*Lse[l,D]-ν*Bse[l,ntrc(D)]            + (η*T+η′*T′)[l,j] -      (K+K′)[l,m] +
                                    -ζ′*Kse[o,D] +                 (K+K′)[o,j] + (α*T+α′*T′)[o,m] == 
                                    -e[l] - h[o])
=#

β = 1/(Material.ϵ_r*Material.ϵ_b)
δ = 1/Material.ϵ_b
σ = 1/(Material.μ_b*Material.μ_r)
ν = 1/Material.μ_b

γ′ = im*κ′

α, α′ = 1/η, 1/η′


#= DB-VSIE (DB-VIE combined with PMCHWT) =#  
eq = @varform (β*I[k,D]-δ*Le[k,D]-δ*Be[ntrc(k),D]  +     γ′*η′*ν*Kh[k,B]                  +       η′*Lt[k,j]+ η′/γ′*Bt[ntrc(k),j] +         -Kt[k,m]                     +
                 -γ′*δ*α′*Ke[l,D]                  +  σ*I[l,B]-ν*Lh[l,B]-ν*Bh[ntrc(l),B]  +          Kt[l,j]                      +       α′*Lt[l,m]+α′/γ′*Bt[ntrc(l),m] +
                   -γ′*δ*Lse[o,D]-δ*Bse[o,ntrc(D)] +    γ′*ν*η′*Ksh[o,B]                  + (η*T+η′*T′)[o,j]                      +     -(K+K′)[o,m]                     +
                -γ′*δ*α′*Kse[p,D]                  +      -γ′*ν*Lsh[p,B]-ν*Bsh[p,ntrc(B)] +      (K+K′)[p,j]                      + (α*T+α′*T′)[p,m]                     == 
                                      -e[o] -h[p])


dbvsie = @discretise eq  D∈X B∈X j∈Y m∈Y k∈X l∈X o∈Y p∈Y

u_n  = solve(dbvsie)


#=
#preconditioner
Mxx = assemble(I,X,X)
iMxx = inv(Matrix(Mxx))

Tzz = assemble(T,Z,Z); println("dual discretisation assembled.")
Nyz = Matrix(assemble(N,Y,Z)); println("duality form assembled.")

iNyz = inv(Nyz); println("duality form inverted.")
NTN = iNyz' * Tzz * iNyz 

M = zeros(Int, 3)
N = zeros(Int, 3)

M[1] = numfunctions(X)
N[1] = numfunctions(X)
M[2] = M[3] = numfunctions(Y)
N[2] = N[3] = numfunctions(Y)

U = BlockArrays.blockedrange(M)
V = BlockArrays.blockedrange(N)

precond = BEAST.ZeroMap{Float32}(U, V)

z1 = LiftedMap(iMxx,Block(1),Block(1),U,V)
z2 = LiftedMap(NTN,Block(2),Block(2),U,V)
z3 = LiftedMap(NTN,Block(3),Block(3),U,V)
precond = precond +ϵ_r*z1 + z2 + z3

A_dvsie_precond = precond*A_dvsie

#GMREs
import IterativeSolvers
cT = promote_type(eltype(A_dvsie), eltype(rhs))
x = PseudoBlockVector{cT}(undef, M)
fill!(x, 0)
x, ch = IterativeSolvers.gmres!(x, A_dvsie, Rhs, log=true,  reltol=1e-4)
fill!(x, 0)
x, ch = IterativeSolvers.gmres!(x, precond*A_dvsie, precond*Rhs, log=true,  reltol=1e-8)
=#

#Post processing
#=
Θ, Φ = range(0.0,stop=2π,length=100), 0.0
ffpoints = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for θ in Θ for ϕ in Φ]

# Don't forget the far field comprises two contributions
ffm = potential(MWFarField3D(κ*im), ffpoints, u_n[m], Y)
ffj = potential(MWFarField3D(κ*im), ffpoints, u_n[j], Y)
ff = -η*im*κ*ffj + im*κ*cross.(ffpoints, ffm)


ffd = potential(VIE.farfield(wavenumber=κ, tau=χ), ffpoints, u_n[D], X)
ff2 = im*κ*ffd

using Plots
plot(xlabel="theta")
plot!(Θ,norm.(ff),label="far field",title="DVSIE")
plot!(Θ,√6.0*norm.(ff),label="far field",title="DVSIE 2")
=#

using Plots
import Plotly
using LinearAlgebra
fcrj, _ = facecurrents(u_n[j],Y)
fcrm, _ = facecurrents(u_n[m],Y)
vsie_j = Plotly.plot([patch(Γ, norm.(fcrj))], Plotly.Layout(title="j DB-VSIE"))
vsie_m = Plotly.plot(patch(Γ, norm.(fcrm)), Plotly.Layout(title="m DB-VSIE"))

#NearField
function nearfield(um,uj,Xm,Xj,κ,η,points,
    Einc=(x->point(0,0,0)),
    Hinc=(x->point(0,0,0)))

    K = BEAST.MWDoubleLayerField3D(wavenumber=κ)
    T = BEAST.MWSingleLayerField3D(wavenumber=κ)

    Em = potential(K, points, um, Xm)
    Ej = potential(T, points, uj, Xj)
    E = -Em + η*Ej + Einc.(points)

    Hm = potential(T, points, um, Xm)
    Hj = potential(K, points, uj, Xj)
    H = 1/η*Hm + Hj + Hinc.(points)

    return E, H
end

Zz = range(-4.0,4.0,length=100)
Yy = range(-4.0,4.0,length=100)
nfpoints = [point(0.0,y,z) for  z in Zz, y in Yy]


import Base.Threads: @spawn
task1 = @spawn nearfield(u_n[m],u_n[j],Y,Y,κ,η,nfpoints,E,H)
task2 = @spawn nearfield(-u_n[m],-u_n[j],Y,Y,κ_in,η_in,nfpoints)

E_ex, H_ex = fetch(task1)
E_in, H_in = fetch(task2)



Enear = BEAST.grideval(nfpoints,β.* u_n[D],X)
Enear = reshape(Enear,100,100)

Hnear = BEAST.grideval(nfpoints,σ.* u_n[B],X)
Hnear = reshape(Hnear,100,100)

contour(real.(getindex.(Enear,1)))
heatmap(Zz, Yy,  real.(getindex.(Enear,1)))

contour(real.(getindex.(Hnear,2)))
heatmap(Zz, Yy,  real.(getindex.(Hnear,2)))

heatmap(Zz, Yy, real.(getindex.(E_in,1)))
heatmap(Zz, Yy, real.(getindex.(E_ex,1)))

contour(real.(getindex.(E_ex,1)))

contour(real.(getindex.(E_in,1)))

contour(real.(getindex.(H_ex,2)))

contour(real.(getindex.(H_in,2)))


Dd = range(1,100,step=1)
plot!(Yy,real.(getindex.(E_in[Dd,50],1)))
plot!(collect(Yy)[2:99],real.(getindex.(E_in[2:99,50],1)),label="DB-VSIE")
plot!(collect(Yy)[2:99],real.(getindex.(H_in[2:99,50],2)),label="DB-VSIE")