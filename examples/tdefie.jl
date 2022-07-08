using CompScienceMeshes, BEAST, LinearAlgebra
Γ = readmesh(joinpath(@__DIR__,"sphere2.in"))

X = raviartthomas(Γ)

Δt, Nt = 0.3, 200
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U

duration = 20 * Δt * 2
delay = 1.5 * duration
amplitude = 1.0
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)

@hilbertspace j
@hilbertspace j′

SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
tdefie = @discretise SL[j′,j] == -1.0E[j′]   j∈V  j′∈W
xefie = solve(tdefie)

import Plots
Plots.plot(xefie[1,:])

import Plotly
fcr, geo = facecurrents(xefie[:,61], X)
Plotly.plot(patch(geo, norm.(fcr)))


using PlotlyJS

n_slices = 150

fcr, geo = facecurrents(xefie[:,1], X)
initial_slice = patch(geo, norm.(fcr))


frames  = Vector{PlotlyFrame}(undef, n_slices)
for k in 1:n_slices
    fcr, geo = facecurrents(xefie[:,k], X)

    frames[k] = PlotlyJS.frame(data=[attr(intensity=norm.(fcr),cmin=0.0,cmax=0.35)],
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
layout = Layout(title_text="Pulse", title_x=0.5,
                width=600,
                height=600,
                sliders=sliders
            )
pl= Plot(initial_slice, layout, frames)

Xefie, Δω, ω0 = fouriertransform(xefie, Δt, 0.0, 2)
ω = collect(ω0 .+ (0:Nt-1)*Δω)
_, i1 = findmin(abs.(ω.-1.0))

ω1 = ω[i1]
ue = Xefie[:,i1] / fouriertransform(gaussian)(ω1)
