mutable struct VIESingleLayerTDIO{T,C} <: RetardedPotential{T}
    "speed of light in medium"
    speed_of_light::T
    "weight of the weakly singular term"
    ws_weight::T
    "weight of the hypersingular term"
    hs_weight::T
    "number of temporal differentiations in the weakly singular term"
    ws_diffs::Int
    "number of temporal differentiations in the hyper singular term"
    hs_diffs::Int
    "contrast current parameter"
    c_current::C
end

function integrand(viop::VIESingleLayerTDIO,kerneldata, tvals, tgeo, bvals, bgeo, )

    gx = @SVector[tvals[i].value for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green
    gradG = kerneldata.gradgreen

    Ty = kerneldata.tau

    α = viop.α
    β = viop.β

    @SMatrix[α * dot(gx[i],Ty*fy[j]) * G - β * dot(dgx[i] * (Ty * fy[j]), gradG) for i in 1:4, j in 1:4]
end



module TDVIE
import ...BEAST

function singlelayer(;speedoflight, numdiffs=0)
	@assert numdiffs == 0
	return BEAST.integrate(BEAST.VIESingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0))
	return BEAST.VIESingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0)
end

end # module TDVIE

export TDVIE

defaultquadstrat(::VIESingleLayerTDIO, tfs, bfs) = nothing


function quaddata(op::VIESingleLayerTDIO, testrefs, trialrefs, timerefs,
    testels, trialels, timeels, ::Nothing)

end


quadrule(op::VIESingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::Nothing) = nothing #WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

function momintegrals!(z, op::VIESingleLayerTDIO,
    g, f, T, #test refspace, trial refspace, time refspace 
    τ, σ, ι, #test element, trial element, time interval
    qr::SauterSchwab3DTimeDomainStrat)

    
end