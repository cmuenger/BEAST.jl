

mutable struct MWSingleLayerTDIO{T} <: RetardedPotential{T}
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
end

function Base.:*(a::Number, op::MWSingleLayerTDIO)
	@info "scalar product a * op (SL)"
	MWSingleLayerTDIO(
		op.speed_of_light,
		a * op.ws_weight,
		a * op.hs_weight,
		op.ws_diffs,
		op.hs_diffs)
end

mutable struct MWDoubleLayerTDIO{T} <: RetardedPotential{T}
    speed_of_light::T
    weight::T
    num_diffs::Int
end

function Base.:*(a::Number, op::MWDoubleLayerTDIO)
	@info "scalar product a * op (DL)"
	MWDoubleLayerTDIO(
		op.speed_of_light,
		a * op.weight,
		op.num_diffs)
end

mutable struct MWDoubleLayerTransposedTDIO{T} <: RetardedPotential{T}
	speed_of_light::T
    weight::T
    num_diffs::Int
end

function Base.:*(a::Number, op::MWDoubleLayerTransposedTDIO)
	@info "scalar product a * op (DL)"
	MWDoubleLayerTransposedTDIO(
		op.speed_of_light,
		a * op.weight,
		op.num_diffs)
end

MWSingleLayerTDIO(;speedoflight) = MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0)
MWDoubleLayerTDIO(;speedoflight) = MWDoubleLayerTDIO(speedoflight, one(speedoflight), 0)


module TDMaxwell3D
import ...BEAST

function singlelayer(;speedoflight, numdiffs=0)
	@assert numdiffs >= 0
	numdiffs == 0 && return BEAST.integrate(BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2,0))
	return BEAST.MWSingleLayerTDIO(speedoflight,-1/speedoflight,-speedoflight,2+numdiffs-1,numdiffs-1)
end

function doublelayer(;speedoflight, numdiffs=0)
	@assert numdiffs >= -1
	numdiffs == -1 && BEAST.integrate(BEAST.MWDoubleLayerTDIO(speedoflight,1.0,0))
	return BEAST.MWDoubleLayerTDIO(speedoflight,1.0,numdiffs)
end

end # module TDMaxwell3D

export TDMaxwell3D

defaultquadstrat(::MWSingleLayerTDIO, tfs, bfs) = nothing

defaultquadstrat(::MWSingleLayerTDIO, tfs, bfs) = SauterQStrat(2,3,6,7,5,5,4,3)

function quaddata(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels, ::Nothing)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    # quadpoints(testrefs, testels, (3,)), bn, ws
      
    println("using Wilton Integration")
    quadpoints(testrefs, testels, (3,)), bn, ws
end


quadrule(op::MWSingleLayerTDIO, testrefs, trialrefs, timerefs,
        p, testel, q, trialel, r, timeel, qd, ::Nothing) = WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])


function quaddata(op::MWSingleLayerTDIO,
                  test_local_space::RefSpace, trial_local_space::RefSpace, time_refs,
                  test_charts, trial_charts, time_els, qs::SauterQStrat)
        
    tqd = quadpoints(test_local_space,  test_charts,  (qs.outer_rule_far,qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far,qs.inner_rule_near))
    leg = (
           _legendre(qs.sauter_schwab_common_vert,0,1),
           _legendre(qs.sauter_schwab_common_edge,0,1),
           _legendre(qs.sauter_schwab_common_face,0,1),)
        
        
    # High accuracy rules (use them e.g. in LF MFIE scenarios)
    # tqd = quadpoints(test_local_space, test_charts, (8,8))
    # bqd = quadpoints(trial_local_space, trial_charts, (8,9))
    # leg = (_legendre(8,a,b), _legendre(10,a,b), _legendre(5,a,b),)
    dmax = numfunctions(time_refs)-1
    bn = binomial.((0:dmax),(0:dmax)')   
    println("using Sauter Schwab Technique")
    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg),bn
end

function quadrule(op::MWSingleLayerTDIO, g::RTRefSpace, f::RTRefSpace, timerefs,  i, τ, j, σ, k, ι, qd,
    qs::SauterQStrat)

  hits = 0
  dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))
  dmin2 = floatmax(eltype(eltype(τ.vertices)))
  for t in τ.vertices
      for s in σ.vertices
          d2 = LinearAlgebra.norm_sqr(t-s)
          dmin2 = min(dmin2, d2)
          hits += (d2 < dtol)
      end
  end

  hits == 3 && return SauterQR(SauterSchwabQuadrature.CommonFace(qd[1].gausslegendre[3]),qd[2])
  hits == 2 && return SauterQR(SauterSchwabQuadrature.CommonEdge(qd[1].gausslegendre[2]),qd[2])
  hits == 1 && return SauterQR(SauterSchwabQuadrature.CommonVertex(qd[1].gausslegendre[1]),qd[2])
  hits == 0 && return SauterQR(SauterSchwabQuadrature.PositiveDistance(qd[1].gausslegendre[1]),qd[2])

  error("Should be never be reached")

  h2 = volume(σ)
  xtol2 = 0.2 * 0.2
  k2 = abs2(op.gamma)
  max(dmin2*k2, dmin2/16h2) < xtol2 && return WiltonSERule(
      qd.tpoints[2,i],
      DoubleQuadRule(
          qd.tpoints[2,i],
          qd.bpoints[2,j],),)
  return DoubleQuadRule(
      qd.tpoints[1,i],
      qd.bpoints[1,j],)
end

struct TransposedStorage{F}
	store::F
end

@inline (f::TransposedStorage)(v,m,n,k) = f.store(v,n,m,k)


# function allocatestorage(op::MWDoubleLayerTDIO, testST, basisST,
# 	::Type{Val{:bandedstorage}},
# 	::Type{LongDelays{:ignore}})

#     # tfs = spatialbasis(testST)
#     # bfs = spatialbasis(basisST)
# 	X, T = spatialbasis(testST), temporalbasis(testST)
# 	Y, U = spatialbasis(basisST), temporalbasis(basisST)

# 	if CompScienceMeshes.refines(geometry(Y), geometry(X))
# 		testST = Y⊗T
# 		basisST = X⊗U
# 	end

#     M = numfunctions(X)
#     N = numfunctions(Y)

#     K0 = fill(typemax(Int), M, N)
#     K1 = zeros(Int, M, N)

#     function store(v,m,n,k)
#         K0[m,n] = min(K0[m,n],k)
#         K1[m,n] = max(K1[m,n],k)
#     end

#     aux = EmptyRP(op.speed_of_light)
#     print("Allocating memory for convolution operator: ")
#     assemble!(aux, testST, basisST, store)
#     println("\nAllocated memory for convolution operator.")

# 	maxk1 = maximum(K1)
# 	bandwidth = maximum(K1 .- K0 .+ 1)
# 	data = zeros(eltype(op), bandwidth, M, N)
# 	Z = SparseND.Banded3D(K0, data, maxk1)
#     store1(v,m,n,k) = (Z[m,n,k] += v)
#     return ()->Z, store1
# end

function assemble!(dl::MWDoubleLayerTDIO, W::SpaceTimeBasis, V::SpaceTimeBasis, store,
    threading=Threading{:multi}; quadstrat=defaultquadstrat(dl,W,V))

	X, T = spatialbasis(W), temporalbasis(W)
	Y, U = spatialbasis(V), temporalbasis(V)
	if CompScienceMeshes.refines(geometry(Y), geometry(X))
		@assert !CompScienceMeshes.refines(geometry(X), geometry(Y))
		W = Y⊗T
		V = X⊗U
		op = MWDoubleLayerTransposedTDIO(dl.speed_of_light, dl.weight, dl.num_diffs)
		assemble!(op, W, V, store)
		return
	end

	P = Threads.nthreads()
	Y, S = spatialbasis(W), temporalbasis(W)
	splits = [round(Int,s) for s in range(0, stop=numfunctions(Y), length=P+1)]

	@info "Starting assembly with $P threads:"
	Threads.@threads for i in 1:P
		lo, hi = splits[i]+1, splits[i+1]
		lo <= hi || continue
		Y_p = subset(Y, lo:hi)
		store2 = (v,m,n,k) -> store(v,lo+m-1,n,k)
		assemble_chunk!(dl, Y_p ⊗ S, V, store2; quadstrat)
	end

	# return assemble_chunk!(dl, W, V, store1)
end

defaultquadstrat(::MWDoubleLayerTDIO, tfs, bfs) = nothing

function quaddata(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::Nothing)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    # quadpoints(testrefs, testels, (3,)), bn, ws
    quadpoints(testrefs, testels, (3,)), bn, ws
end

quadrule(op::MWDoubleLayerTDIO, testrefs, trialrefs, timerefs,
    p, testel, q, trialel, r, timeel, qd, quadstrat::Nothing) =
        WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])


defaultquadstrat(::MWDoubleLayerTransposedTDIO, tfs, bfs) = nothing

function quaddata(op::MWDoubleLayerTransposedTDIO,
		testrefs, trialrefs, timerefs,
        testels, trialels, timeels, quadstrat::Nothing)

    dmax = numfunctions(timerefs)-1
    bn = binomial.((0:dmax),(0:dmax)')

    V = eltype(testels[1].vertices)
    ws = WiltonInts84.workspace(V)
    quadpoints(testrefs, testels, (3,)), bn, ws
end

quadrule(op::MWDoubleLayerTransposedTDIO, testrefs, trialrefs, timerefs,
    p, testel, q, trialel, r, timeel, qd, quadstrat::Nothing) =
        WiltonInts84Strat(qd[1][1,p],qd[2],qd[3])

function momintegrals!(z, op::MWDoubleLayerTransposedTDIO,
	g, f, T, τ, σ, ι, qr::WiltonInts84Strat)

	op1 = MWDoubleLayerTDIO(op.speed_of_light, op.weight, op.num_diffs)
	momintegrals!(z, op1, g, f, T, τ, σ, ι, qr::WiltonInts84Strat)
	w = similar(z)
	permutedims!(w, z, [2,1,3])

end

@inline function tmRoR(d, iG, bns)
    sgn = isodd(d) ? -1 : 1
    r = sgn * iG[d+2]
end

# build
# ``\int (D-R)^d/R (y-b) dy`` from
# ``(ξ-b) \int R^k dy`` and
# ``\int R^k (y-ξ) dy``
@inline function tmRoRf(d, iG, iGξy, bξ, bns)
    sgn = isodd(d) ? -1 : 1
    iGf = iGξy[d+2] + bξ * iG[d+2]
    r = sgn * iGf
end

"""
    Q = qd(T,dh,::Val{N})

Q[k] is the factor in front resulting from differentiating t^(k-1) dh times.
"""
@generated function qh(::Type{T},dh,n::Type{Val{N}}) where {N,T}
    xp = quote end
    for k in 1:N
        qk = Symbol(:q,k)
        d = k-1
        xp1 = quote
            $(qk) = one($T)
            for v in 0 : dh-1
                $(qk) *= ($d-v)
            end
        end
        append!(xp.args, xp1.args)
    end
    xp1 = :(())
    for k in 1:N
        qk = Symbol(:q,k)
        push!(xp1.args, :($qk))
    end
    push!(xp.args, xp1)
    return xp
end


function innerintegrals!(zl, op::MWSingleLayerTDIO,
        p, # test_point, test_time
        U, V, W, # local_test_space, local_trial_space, local_temporal_space
        τ, σ, ι, # test_element, trial_element, spherial_shell
        qr, w)   # inner_quadrature_rule, outer_quadrature_weight

	T = typeof(w)

    sol = op.speed_of_light
    #Rmax = sol * tmax

    dx = w
    x = cartesian(p)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n #projection of test point on the plane of trial element

    r = ι[1]
    R = ι[2]

    @assert r < R
    @assert degree(W) <= 3

    ∫G, ∫Gξy, = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{2},qr.workspace)

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4π
	α = αg * αf * αG * op.ws_weight * dx
	β = 4 * αg * αf * αG * op.hs_weight * dx

	ds = op.ws_diffs
	dh = op.hs_diffs

    qhs = qh(T,dh,Val{4})
    qss = qh(T,ds,Val{4})

    bn = qr.binomials

    #solpowers = collect(sol^p for p ∈ 0:numfunctions(W)-1)
    sol2 = sol*sol
    sol3 = sol2*sol
    sol4 = sol3*sol
    sol5 = sol4*sol
    solpowers = (one(sol), sol, sol2, sol3, sol4, sol5)

    for i in 1 : numfunctions(U)
        a = τ[i]
        g = (x-a)
        for j in 1 : numfunctions(V)
            b = σ[j]; bξ = ξ-b
            for k in 1 : numfunctions(W)
				d = k-1 # ranges from 0 to numfunctions(W)-1
				sgn = isodd(d) ? -1 : 1
				# hyper singular contribution
				if d >= dh
                    @assert dh == 0
                    q = qhs[k]
                    Ih = tmRoR(d-dh, ∫G, bn) # \int (-R)^(d-dh)/R dy
                    #zl[i,j,k] += β * q * Ih / sol^(d-dh)
                    zl[i,j,k] += β * q * Ih / solpowers[d-dh+1]
				end
				# weakly singular contribution
				if d >= ds
                    q = qss[k]
                    Is = tmRoRf(d-ds, ∫G, ∫Gξy, bξ, bn) # \int (cTmax-R)^(d-ds)/R (y-b) dy
                    #zl[i,j,k] += α * q * (g ⋅ Is) / sol^(d-ds)
                    zl[i,j,k] += α * q * (g ⋅ Is) / solpowers[d-ds+1]
				end
            end
        end
    end

end

struct MWTDSL3DIntegrand{C,D,O,L,M}
    test_triangular_element::C
    trial_triangular_element::C
    time_element::D
    op::O
    test_local_space::L
    trial_local_space::L
    time_local_space::M
end

function (igd::MWTDSL3DIntegrand)(u,v)

        x = neighborhood(igd.test_triangular_element,u)
        y = neighborhood(igd.trial_triangular_element,v)

        r = cartesian(x) - cartesian(y)
        R = norm(r)

        T = typeof(R)
    
        rmin = igd.time_element[1]
        rmax = igd.time_element[2]

        if R <= rmin || R >= rmax
            return @SArray[ 0.0 for i in 1:3, j in 1:3, k in 1:3]
        end

        iR = 1 / R

        Rc = R*igd.op.speed_of_light

        ds = igd.op.ws_diffs
        dh = igd.op.hs_diffs
    
        qhs = qh(T,dh,Val{4})
        qws = qh(T,ds,Val{4})
        
        fx = igd.test_local_space(x)
        gy = igd.trial_local_space(y)

        j = jacobian(x) * jacobian(y)

        αG = 1 / 4π
        α = αG*iR * igd.op.ws_weight*j 
        β = αG*iR * igd.op.hs_weight*j

        Gws =   @SVector[d>= ds ?  qws[d+1]*(-Rc)^(d-ds) : 0.0 for d in 0:2]
        Ghs =   @SVector[d>= dh ?  qhs[d+1]*(-Rc)^(d-dh) : 0.0 for d in 0:2]
       
        @SArray[ α*fx[i].value'*Gws[k]*gy[j].value + β*fx[i].divergence*Ghs[k]*gy[j].divergence for i in 1:3, j in 1:3, k in 1:3]
end

function momintegrals!(z, op::MWSingleLayerTDIO, U, V, W, τ, σ, ι, qr::SauterQR)

    I, J, K, L = SauterSchwabQuadrature.reorder(τ.vertices,σ.vertices, qr.strat)

    tgeo  = simplex(
        τ.vertices[I[1]],
        τ.vertices[I[2]],
        τ.vertices[I[3]])

    bgeo = simplex(
        σ.vertices[J[1]],
        σ.vertices[J[2]],
        σ.vertices[J[3]])


    r = ι[1]
    R = ι[2]

    @assert r < R

    @assert degree(W) == 2


    igd = MWTDSL3DIntegrand(tgeo, bgeo, ι, op, U, V, W)

    G =  SauterSchwabQuadrature.sauterschwab_parameterized(igd, qr.strat)

    for j ∈ 1:3, i ∈ 1:3
        z[i,j,:] += G[K[i],L[j],:]
    end

end


function innerintegrals!(z, op::MWDoubleLayerTDIO,
    p,
    U, V, W,
    τ, σ, ι,
    qr, w)

	T = typeof(w)

    sol = op.speed_of_light

    dx = w
    x = cartesian(p)
    n = cross(σ[1]-σ[3],σ[2]-σ[3])
    n /= norm(n)
    ξ = x - ((x-σ[1]) ⋅ n) * n

    r = ι[1]
    R = ι[2]

    @assert r < R
    @assert degree(W) <= 3

    #N = max(degree(W), 0)
    ∫G, ∫Gξy, ∫∇G = WiltonInts84.wiltonints(σ[1],σ[2],σ[3],x,r,R,Val{2},qr.workspace)
	@assert isapprox(∫∇G[2] , point(0,0,0), atol=1e-8)

    αg = 1 / volume(τ) / 2
	αf = 1 / volume(σ) / 2
	αG = 1 / 4 / π
	α = αG * op.weight * dx # * αg * αf

	ds = op.num_diffs

    @inline function tmGRoR(d, iGG)
        sgn = isodd(d) ? -1 : 1
        r = sgn * iGG[d+1]
    end

    Ux = U(p)
    Vx = αf * @SVector[(x-σ[1]), (x-σ[2]), (x-σ[3])]

    for i in 1 : numfunctions(U)
        # a = τ[i]
        # g = αg * (x-τ[i])
        g = Ux[i].value
        for j in 1 : numfunctions(V)
            # b = σ[j]
            # f = αf * (x-σ[j])
            # f = Vx[j].value
            f = Vx[j]
            fxg = f × g
            for k in 1 : numfunctions(W)
				d = k-1
				sgn = isodd(d) ? -1 : 1
				if d >= ds
					q = one(T)
					for p in 0 : ds-1
						q *= (d-p)
					end
                    # @assert q == 1
                    z[i,j,k] += -α * q * ( fxg ⋅ tmGRoR(d-ds, ∫∇G) ) / sol^(d-ds)
				end
            end
        end
    end

end


