using LinearAlgebra


function addf!(fn::Vector{<:Shape}, x::Vector, space::Space, idcs::Vector{Int})
    for (m,bf) in enumerate(space.fns)
        for sh in bf
            cellid = idcs[sh.cellid]
            BEAST.add!(fn, cellid, sh.refid, sh.coeff * x[m])
        end
    end
end

function Base.in(mesh::CompScienceMeshes.AbstractMesh)
    cells_mesh = sort.(mesh)
    function f(cell)
        sort(cell) in cells_mesh
    end
end


function builddual2form(support, port, dirichlet, prt_fluxes)

    faces = CompScienceMeshes.skeleton_fast(support,2)
    edges = CompScienceMeshes.skeleton_fast(support,1)
    verts = CompScienceMeshes.skeleton_fast(support,0)
    @assert length(verts) - length(edges) + length(faces) - length(support) == 1
    bndry = boundary(support)

    cells_bndry = sort.(bndry)
    dirbnd = submesh(dirichlet) do face
        sort(face) in cells_bndry ? true : false
    end

    bnd_dirbnd = boundary(dirbnd)
    edges_dirbnd = CompScienceMeshes.skeleton_fast(dirbnd,1)
    cells_bnd_dirbnd = sort.(bnd_dirbnd)
    int_edges_dirbnd = submesh(edges_dirbnd) do edge
        sort(edge) in cells(bnd_dirbnd) ? false : true
    end

    num_faces_on_port = 0
    num_faces_on_dirc = 0

    cells_port = sort.(port)
    cells_dirbnd = sort.(dirbnd)
    int_faces = submesh(faces) do face
        sort(face) in cells_port && return false
        sort(face) in cells_dirbnd && return true
        sort(face) in cells_bndry && return false
        return true
    end

    bnd_edges = CompScienceMeshes.skeleton_fast(bndry,1)
    prt_edges = CompScienceMeshes.skeleton_fast(port,1)

    cells_int_edges_dirbnd = sort.(int_edges_dirbnd)
    cells_bnd_edges = sort.(bnd_edges)
    cells_prt_edges = sort.(prt_edges)

    int_edges = submesh(edges) do edge
        sort(edge) in cells_int_edges_dirbnd && return true
        sort(edge) in cells_bnd_edges && return false
        sort(edge) in cells_prt_edges && return false
        return true
    end

    RT_int = nedelecd3d(support, int_faces)
    RT_prt = nedelecd3d(support, port)
    L0_int = nedelecc3d(support, int_edges)

    Id = BEAST.Identity()
    div_RT_int = divergence(RT_int)
    div_RT_prt = divergence(RT_prt)
    D = assemble(Id, div_RT_int, div_RT_int)
    Q = assemble(Id, div_RT_int, div_RT_prt)
    d = -Q * prt_fluxes

    curl_L0_int = curl(L0_int)
    div_curl_L0_int = divergence(curl_L0_int)
    ZZ = real(assemble(Id, div_curl_L0_int, div_curl_L0_int))
    @assert isapprox(norm(ZZ), 0.0, atol=1e-8)
    C = assemble(Id, curl_L0_int, RT_int)
    c = real(assemble(Id, curl_L0_int, RT_prt)) * prt_fluxes

    T = eltype(D)
    nz = length(c)
    QQ = [D C'; C zeros(T,nz,nz)]
    qq = [d;c]
    x = (QQ \ qq)[1:end-nz]

    if !isapprox(C*x, c, atol=1e-8) || !isapprox(D*x, d, atol=1e-6)
        @show norm(D*x-d)
        @show norm(C*x-c)
        @show rank(C)
        @show size(C,1)
        error("error")
    end

    return RT_int, RT_prt, x, prt_fluxes

end


function dual2forms_init(Tetrs)

    tetrs = barycentric_refinement(Tetrs)
    v2t, v2n = CompScienceMeshes.vertextocellmap(tetrs)
    bnd = boundary(tetrs)

    return tetrs, bnd, v2t, v2n
end


function dual2forms(Tetrs, Edges, Dir)
    tetrs, bnd, v2t, v2n = dual2forms_init(Tetrs)
    dual2forms_body(Tetrs, Edges, Dir, tetrs, bnd, v2t, v2n)
end

function dual2forms_body(Tetrs, Edges, Dir, tetrs, bnd, v2t, v2n)

    T = coordtype(Tetrs)
    bfs = Vector{Vector{Shape{T}}}(undef, numcells(Edges))
    pos = Vector{vertextype(Edges)}(undef, numcells(Edges))

    gpred = CompScienceMeshes.overlap_gpredicate(Dir)
    dirichlet = submesh(bnd) do face
        gpred(chart(bnd,face))
    end

    for (F,Edge) in enumerate(cells(Edges))

        println("Constructing dual 2-forms: $F out of $(length(Edges)).")
        bfs[F] = Vector{Shape{T}}()

        pos[F] = cartesian(center(chart(Edges,Edge)))
        port_vertex_idx = argmin(norm.(vertices(tetrs) .- Ref(pos[F])))

        ptch_idcs1 = v2t[Edge[1],1:v2n[Edge[1]]]
        ptch_idcs2 = v2t[Edge[2],1:v2n[Edge[2]]]

        patch1 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs1])
        patch2 = Mesh(vertices(tetrs), cells(tetrs)[ptch_idcs2])
        patch = CompScienceMeshes.union(patch1, patch2)

        patch_bnd = boundary(patch1)

        bnd_patch1 = boundary(patch1)
        bnd_patch2 = boundary(patch2)
        port = submesh(in(bnd_patch2), bnd_patch1)

        prt_fluxes = ones(T, numcells(port)) / numcells(port)
        tgt = vertices(Edges)[Edge[1]] - vertices(Edges)[Edge[2]]
        for (i,face) in enumerate(cells(port))
            chrt = chart(port, face)
            prt_fluxes[i] *= sign(dot(normal(chrt), tgt))
        end

        # RT_int, RT_prt, x_int, x_prt = builddual2form(patch, port, dirichlet, prt_fluxes)
        # ptch_idcs = [ptch_idcs1; ptch_idcs2]
        # addf!(bfs[F], x_int, RT_int, ptch_idcs)
        # addf!(bfs[F], x_prt, RT_prt, ptch_idcs)

        RT1_int, RT1_prt, x1_int, x_prt = builddual2form(
            patch1, port, dirichlet, +prt_fluxes)
        RT2_int, RT2_prt, x2_int, x_prt = builddual2form(
            patch2, port, dirichlet, prt_fluxes)

        addf!(bfs[F], x1_int, RT1_int, ptch_idcs1)
        addf!(bfs[F], +x_prt, RT1_prt, ptch_idcs1)

        addf!(bfs[F], x2_int, RT2_int, ptch_idcs2)
        addf!(bfs[F], x_prt, RT2_prt, ptch_idcs2)

    end

    NDLCDBasis(tetrs, bfs, pos)
end




function builddual1form(supp, port, dir, x0)

    Id = BEAST.Identity()

    rim = boundary(dir)
    bnd = boundary(supp)
    supp_edges = CompScienceMeshes.skeleton_fast(supp,1)
    supp_nodes = CompScienceMeshes.skeleton_fast(supp,0)
    dir_edges = CompScienceMeshes.skeleton_fast(dir,1)
    bnd_edges = CompScienceMeshes.skeleton_fast(bnd,1)
    bnd_nodes = CompScienceMeshes.skeleton_fast(bnd,0)
    prt_nodes = CompScienceMeshes.skeleton_fast(port,0)

    srt_rim = sort.(rim)
    srt_port = sort.(port)
    srt_dir_edges = sort.(dir_edges)
    srt_bnd_edges = sort.(bnd_edges)
    srt_bnd_verts = sort.(bnd_nodes)
    srt_prt_verts = sort.(prt_nodes)

    int_edges = submesh(supp_edges) do edge
        sort(edge) in srt_port && return false
        sort(edge) in srt_rim && return false
        sort(edge) in srt_dir_edges && return true
        sort(edge) in srt_bnd_edges && return false
        return true
    end
    # @assert length(supp_nodes) - length(supp_edges) + length(skeleton(supp,2)) - length(supp) == 1

    Nd_prt = BEAST.nedelecc3d(supp, port)
    Nd_int = BEAST.nedelecc3d(supp, int_edges)

    curl_Nd_int = curl(Nd_int)
    A = assemble(Id, curl_Nd_int, curl_Nd_int)
    N = nullspace(A)
    a = -assemble(Id, curl(Nd_int), curl(Nd_prt)) * x0
    x1 = pinv(A) * a
    # x1 = A \ a

    int_verts = submesh(supp_nodes) do vert
        sort(vert) in srt_bnd_verts && return false
        sort(vert) in srt_prt_verts && return false
        return true
    end

    L0_int = lagrangec0d1(supp, int_verts)
    grad_L0_int = BEAST.gradient(L0_int)
    @assert numfunctions(grad_L0_int) == numfunctions(L0_int)

    B = assemble(Id, grad_L0_int, Nd_int)
    b = -assemble(Id, grad_L0_int, Nd_prt) * x0
    p = (B*N) \ (b-B*x1)
    x1 = x1 + N*p

    @assert isapprox(A*x1, a, atol=1e-8)
    @assert isapprox(B*x1, b, atol=1e-8)

    return Nd_int, Nd_prt, x1, x0
end

function dual1forms(Tetrs, Faces)

    T = coordtype(Tetrs)
    P = vertextype(Tetrs)
    S = Shape{T}

    tetrs = CompScienceMeshes.barycentric_refinement(Tetrs)
    bnd_tetrs = boundary(tetrs)
    srt_bnd_tetrs = sort.(bnd_tetrs)

    fns = Vector{Vector{S}}()
    pos = Vector{P}()
    for (F,Face) in enumerate(Faces)
        @show F

        po = cartesian(center(chart(Faces, Face)))

        idcs1 = [i for (i,tet) in enumerate(tetrs) if Face[1] in tet]
        idcs2 = [i for (i,tet) in enumerate(tetrs) if Face[2] in tet]
        idcs3 = [i for (i,tet) in enumerate(tetrs) if Face[3] in tet]
        idcs = vcat(idcs1, idcs2, idcs3)
        # @show idcs

        supp1 = Mesh(vertices(tetrs), cells(tetrs)[idcs1])
        supp2 = Mesh(vertices(tetrs), cells(tetrs)[idcs1])
        supp3 = Mesh(vertices(tetrs), cells(tetrs)[idcs1])

        supp1 = submesh(tet -> Face[1] in tet, tetrs.mesh)
        supp2 = submesh(tet -> Face[2] in tet, tetrs.mesh)
        supp3 = submesh(tet -> Face[3] in tet, tetrs.mesh)

        supp = CompScienceMeshes.union(supp1, supp2)
        supp = CompScienceMeshes.union(supp, supp3)

        dir = submesh(face -> sort(face) in srt_bnd_tetrs, boundary(supp))

        port = skeleton(supp1, 1)
        port = submesh(edge -> sort(edge) in sort.(skeleton(supp2,1)), port)
        port = submesh(edge -> sort(edge) in sort.(skeleton(supp3,1)), port)
        @assert 1 ≤ length(port) ≤ 2

        x0 = ones(length(port)) / length(port)
        for (i,edge) in enumerate(port)
            tgt = tangents(center(chart(port,edge)),1)
            dot(normal(chart(Faces,Face)), tgt) < 0 && (x0[i] *= -1)
        end

        Nd_int, Nd_prt, x_int, x_prt = builddual1form(supp, port, dir, x0)
        @show norm(x_int)
        @show norm(x_prt)

        fn = Vector{S}()
        addf!(fn, x_prt, Nd_prt, idcs)
        addf!(fn, x_int, Nd_int, idcs)

        push!(fns, fn)
        push!(pos, po)
    end

    NDLCCBasis(tetrs, fns, pos)
end
