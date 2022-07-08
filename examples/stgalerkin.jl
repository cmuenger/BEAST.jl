using CompScienceMeshes, BEAST, LinearAlgebra
o, x, y, z = euclidianbasis(3)

# D, Δx = 1.0, 0.35
# Γ = meshsphere(D, Δx)

Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))
Γ = meshsphere(1.0, 0.25)
X = raviartthomas(Γ)
#Δt, Nt = 0.08, 400
Δt, Nt = 0.3, 100
T = timebasisc0d1(Δt, Nt)
U = timebasiscxd0(Δt, Nt)

#T = timebasisshiftedlagrange(Δt, Nt, 1)
#U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

width, delay, scaling = 8.0, 12.0, 1.0
gaussian = creategaussian(width, delay, scaling)

direction, polarisation = ẑ, x̂
E = BEAST.planewave(polarisation, direction, derive(gaussian), 1.0)


BEAST.defaultquadstrat(::MWSingleLayerTDIO, tfs, bfs) = nothing

#BEAST.defaultquadstrat(::MWSingleLayerTDIO, tfs, bfs) = BEAST.SauterQStrat(3,3,3,3,3,3,3,3)


@hilbertspace j; @hilbertspace j′
T = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
tdefie = @discretise T[j′,j] == -1E[j′]   j∈V  j′∈W
xefie = solve(tdefie)


using PlotlyJS

n_slices = 100

fcr, geo = facecurrents(xefie[:,1], X)
initial_slice = patch(geo, norm.(fcr))


frames  = Vector{PlotlyFrame}(undef, n_slices)
for k in 1:n_slices
    fcr, geo = facecurrents(xefie[:,k], X)

    frames[k] = PlotlyJS.frame(data=[attr(intensity=norm.(fcr),cmin=0.0,cmax=0.45)],
                                 name="fr$k",
                                 traces=[0])
end    

sliders = [attr(steps = [attr(method= "animate",
                              args= [["fr$k"],                           
                              attr(mode= "immediate",
                                   frame= attr(duration=20, redraw= true),
                                   transition=attr(duration= 0))
                                 ],
                              label="$k"
                             ) for k in 1:n_slices], 
                active=17,
                transition= attr(duration= 0 ),
                x=0, # slider starting position  
                y=0, 
                currentvalue=attr(font=attr(size=12), 
                                  prefix="Timestep: ", 
                                  visible=true, 
                                  xanchor= "center"
                                 ),  
               len=1.0) #slider length
           ];
layout = Layout(title_text="Pulse SauterSchwab", title_x=0.5,
                width=600,
                height=600,
                sliders=sliders
            )
pl= Plot(initial_slice, layout, frames)



#=
Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1]/ fouriertransform(gaussian)(ω1)
=#


s1 = PlotlyJS.scatter(y=abs.(xefie[123,:]),mode="lines",name="Wilton")
s2 = PlotlyJS.scatter(y=abs.(xefie2[123,:]),mode="lines",name="SauterSchwab")
plot([s1,s2], Layout(yaxis_type="log"))

