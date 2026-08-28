#!/usr/bin/env julia
using PhyloNetworks
using JSON

# ============================================================
# Helpers: endpoints, direction, traversal
# ============================================================

edge_endpoints(e) = (e.node[1], e.node[2])

function other_endpoint(e, u)
    n1, n2 = edge_endpoints(e)
    return (u === n1) ? n2 : n1
end

function is_directed(u, v, e)
    n1, n2 = edge_endpoints(e)
    if e.ischild1
        return (u === n2) && (v === n1)
    else
        return (u === n1) && (v === n2)
    end
end

function directed_children_all(u)
    kids = PhyloNetworks.Node[]
    for e in u.edge
        v = other_endpoint(e, u)
        if is_directed(u, v, e)
            push!(kids, v)
        end
    end
    return kids
end

function descendant_leaves(start::PhyloNetworks.Node; edges_allowed::Union{Nothing,Set}=nothing)
    stack = PhyloNetworks.Node[start]
    seen = Set{PhyloNetworks.Node}()
    leaves = String[]

    while !isempty(stack)
        u = pop!(stack)
        if u in seen
            continue
        end
        push!(seen, u)

        if u.leaf
            push!(leaves, u.name)
            continue
        end

        for e in u.edge
            if edges_allowed !== nothing && !(e in edges_allowed)
                continue
            end
            v = other_endpoint(e, u)
            if is_directed(u, v, e)
                push!(stack, v)
            end
        end
    end

    sort!(unique!(leaves))
    return leaves
end

function component_leaves_undirected(start::PhyloNetworks.Node, blocked_edge)
    stack = PhyloNetworks.Node[start]
    seen = Set{PhyloNetworks.Node}()
    leaves = String[]

    while !isempty(stack)
        u = pop!(stack)
        if u in seen
            continue
        end
        push!(seen, u)

        if u.leaf
            push!(leaves, u.name)
        end

        for e in u.edge
            if e === blocked_edge
                continue
            end
            v = other_endpoint(e, u)
            push!(stack, v)
        end
    end

    sort!(unique!(leaves))
    return leaves
end

# ============================================================
# Blob extraction + clades
# ============================================================

function nodes_in_blob(blob_edges::Set)
    ns = Set{PhyloNetworks.Node}()
    for e in blob_edges
        n1, n2 = edge_endpoints(e)
        push!(ns, n1)
        push!(ns, n2)
    end
    return ns
end

function blobs_as_edgesets(net::PhyloNetworks.HybridNetwork; min_vertices::Int=5)
    bcc = PhyloNetworks.biconnectedcomponents(net)
    blobs = Vector{Set}()

    for comp in bcc
        bedges = Set(comp)
        ns = nodes_in_blob(bedges)
        if length(ns) >= min_vertices
            push!(blobs, bedges)
        end
    end

    return blobs
end

function hybrid_nodes_in_blob(blob_edges::Set)
    ns = nodes_in_blob(blob_edges)
    return [u for u in ns if u.hybrid]
end

function hybrid_clade(net, blob_edges::Set, u::PhyloNetworks.Node)
    kids = directed_children_all(u)
    @assert length(kids) == 1 "Hybrid node did not have a unique directed child in the full network"
    return descendant_leaves(kids[1])
end

function boundary_clades(net, blob_edges::Set)
    ns = nodes_in_blob(blob_edges)
    out = Vector{Tuple{Int,Vector{String}}}()

    for u in ns
        for e in u.edge
            if e in blob_edges
                continue
            end
            v = other_endpoint(e, u)
            ls = component_leaves_undirected(v, e)
            push!(out, (u.number, ls))
        end
    end

    uniq = Dict{String,Tuple{Int,Vector{String}}}()
    for (anch, ls) in out
        key = join(ls, ",")
        uniq[key] = (anch, ls)
    end

    return collect(values(uniq))
end

# ============================================================
# Level-1: cycle extraction + circular order of boundary clades
# ============================================================

function blob_adj(blob_edges::Set)
    adj = Dict{PhyloNetworks.Node, Vector{PhyloNetworks.Node}}()

    for e in blob_edges
        u, v = edge_endpoints(e)
        get!(adj, u, PhyloNetworks.Node[])
        get!(adj, v, PhyloNetworks.Node[])
        push!(adj[u], v)
        push!(adj[v], u)
    end

    return adj
end

function extract_unique_cycle_nodes(adj::Dict{PhyloNetworks.Node, Vector{PhyloNetworks.Node}})
    start = first(keys(adj))
    parent = Dict{PhyloNetworks.Node, PhyloNetworks.Node}()
    seen = Set{PhyloNetworks.Node}()

    NodeT = PhyloNetworks.Node
    StackItem = Tuple{NodeT, Union{Nothing,NodeT}}
    stack = StackItem[(start, nothing)]

    while !isempty(stack)
        u, p = pop!(stack)

        if u in seen
            continue
        end

        push!(seen, u)

        if p !== nothing
            parent[u] = p
        end

        for v in adj[u]
            if p !== nothing && v === p
                continue
            end

            if v in seen
                au = Set{NodeT}()
                x = u
                push!(au, x)

                while haskey(parent, x)
                    x = parent[x]
                    push!(au, x)
                end

                path_v = NodeT[v]
                y = v

                while !(y in au) && haskey(parent, y)
                    y = parent[y]
                    push!(path_v, y)
                end

                lca = y

                path_u = NodeT[u]
                x = u

                while x !== lca
                    x = parent[x]
                    push!(path_u, x)
                end

                cyc = vcat(path_u, reverse(path_v[1:end-1]))
                return cyc
            else
                push!(stack, (v, u))
            end
        end
    end

    error("No cycle found; check blob extraction / level-1 assumption.")
end

function boundary_order_level1(net, blob_edges::Set)
    adj = blob_adj(blob_edges)
    cycle_nodes = extract_unique_cycle_nodes(adj)

    boundary = boundary_clades(net, blob_edges)
    by_anchor = Dict{Int, Vector{String}}()

    for (anch, ls) in boundary
        lab = join(ls, ",")
        push!(get!(by_anchor, anch, String[]), lab)
    end

    for v in values(by_anchor)
        sort!(v)
    end

    order = String[]

    for u in cycle_nodes
        append!(order, get(by_anchor, u.number, String[]))
    end

    extra = setdiff(collect(keys(by_anchor)), [u.number for u in cycle_nodes])
    sort!(extra)

    for a in extra
        append!(order, by_anchor[a])
    end

    return order
end

# ============================================================
# Circular Kendall-tau distance
# ============================================================

function kendall_tau(p::Vector{String}, q::Vector{String})
    n = length(p)
    @assert length(q) == n

    pos = Dict{String,Int}()
    for (i, lab) in enumerate(q)
        pos[lab] = i
    end

    a = [pos[lab] for lab in p]
    bit = zeros(Int, n)

    function bit_add(i, v)
        while i <= n
            bit[i] += v
            i += i & -i
        end
    end

    function bit_sum(i)
        s = 0
        while i > 0
            s += bit[i]
            i -= i & -i
        end
        return s
    end

    inv = 0
    for i in 1:n
        inv += (i - 1) - bit_sum(a[i])
        bit_add(a[i], 1)
    end

    return inv
end

function circular_kendall_tau(p::Vector{String}, q::Vector{String})
    n = length(p)
    best = typemax(Int)

    for rev in (false, true)
        qq = rev ? reverse(q) : q
        for shift in 0:n-1
            rot = vcat(qq[shift+1:end], qq[1:shift])
            best = min(best, kendall_tau(p, rot))
        end
    end

    return best
end

function circular_kendall_tau_common(p::Vector{String}, q::Vector{String})
    common = intersect(Set(p), Set(q))
    p2 = [x for x in p if x in common]
    q2 = [x for x in q if x in common]

    if length(p2) <= 1
        return 0, length(p2)
    end

    return circular_kendall_tau(p2, q2), length(p2)
end

function normalize_kendall(d::Int, n::Int)
    if n < 2
        return 0.0
    end
    return d / (n * (n - 1) / 2)
end

# ============================================================
# Main analysis / comparison
# ============================================================

function analyze_network(path::String)
    net = readnewick(read(path, String))
    blobs = blobs_as_edgesets(net; min_vertices=5)

    results = []

    for (bi, bedges) in enumerate(blobs)
        hybrids = hybrid_nodes_in_blob(bedges)
        hclades = [hybrid_clade(net, bedges, h) for h in hybrids]

        bcl = boundary_clades(net, bedges)
        border_labels = boundary_order_level1(net, bedges)

        push!(results, (
            blob_index = bi,
            num_blob_edges = length(bedges),
            num_blob_vertices = length(nodes_in_blob(bedges)),
            boundary_count = length(border_labels),
            hybrid_nodes = [h.number for h in hybrids],
            hybrid_clades = hclades,
            boundary_clades = [x[2] for x in bcl],
            boundary_order_labels = border_labels
        ))
    end

    return net, results
end

function ensure_json_path(p::String)
    endswith(lowercase(p), ".json") ? p : (p * ".json")
end

function write_summary_json(outpath::String, summary::Dict{String,Any})
    outpath = ensure_json_path(outpath)

    open(outpath, "w") do io
        JSON.print(io, summary, 4)
        println(io)
    end

    return outpath
end

function clade_label(clade::Vector{String})
    return join(sort(clade), ",")
end

function hybrid_clade_labels(r)
    return sort([clade_label(hc) for hc in r.hybrid_clades])
end

function boundary_key(r)
    return sort([clade_label(c) for c in r.boundary_clades])
end

function compare_two(
    nwk1::String,
    nwk2::String;
    outjson::Union{Nothing,String}=nothing,
    min_boundary_clades::Int=5
)
    _, res1_all = analyze_network(nwk1)
    _, res2_all = analyze_network(nwk2)

    res1 = [r for r in res1_all if r.boundary_count >= min_boundary_clades]
    res2 = [r for r in res2_all if r.boundary_count >= min_boundary_clades]

    boundary_counts_net1_all = [r.boundary_count for r in res1_all]
    boundary_counts_net2_all = [r.boundary_count for r in res2_all]
    boundary_counts_net1 = [r.boundary_count for r in res1]
    boundary_counts_net2 = [r.boundary_count for r in res2]

    # ---- global hybrid clade sets ----
    function all_hybrid_clade_labels(res)
        labs = String[]
        for r in res
            append!(labs, hybrid_clade_labels(r))
        end
        return labs
    end

    set1 = Set(all_hybrid_clade_labels(res1))
    set2 = Set(all_hybrid_clade_labels(res2))

    hybrid_unique_net1_only = length(setdiff(set1, set2))
    hybrid_unique_net2_only = length(setdiff(set2, set1))
    hybrid_shared_count = length(intersect(set1, set2))

    # ---- match blobs by multiset of boundary clades ----
    dict2 = Dict{String, Vector{Any}}()

    for r in res2
        k = join(boundary_key(r), "||")
        push!(get!(dict2, k, Any[]), r)
    end

    raw_ds = Float64[]
    norm_ds = Float64[]

    matched_cycles = 0
    missing_cycles = 0
    matched_blob_hybrid_info = Vector{Dict{String,Any}}()

    for r1 in res1
        k = join(boundary_key(r1), "||")

        if !haskey(dict2, k) || isempty(dict2[k])
            missing_cycles += 1

            h1 = hybrid_clade_labels(r1)

            push!(matched_blob_hybrid_info, Dict{String,Any}(
                "matched" => false,
                "net1_blob_index" => r1.blob_index,
                "net2_blob_index" => nothing,
                "boundary_clades" => boundary_key(r1),

                "net1_hybrid_nodes" => r1.hybrid_nodes,
                "net2_hybrid_nodes" => Any[],

                "net1_hybrid_clades" => h1,
                "net2_hybrid_clades" => Any[],

                "correct_hybrid_clades" => Any[],
                "net1_only_hybrid_clades" => h1,
                "net2_only_hybrid_clades" => Any[],

                "has_correct_hybrid" => false,
                "num_correct_hybrids" => 0,
                "num_net1_only_hybrids" => length(h1),
                "num_net2_only_hybrids" => 0,

                "kendall_raw" => nothing,
                "kendall_norm" => nothing
            ))

            continue
        end

        r2 = popfirst!(dict2[k])
        matched_cycles += 1

        p = r1.boundary_order_labels
        q = r2.boundary_order_labels

        d = 0
        nused = 0

        if Set(p) == Set(q) && length(p) == length(q)
            d = circular_kendall_tau(p, q)
            nused = length(p)
        else
            d, nused = circular_kendall_tau_common(p, q)
        end

        dn = normalize_kendall(d, nused)

        push!(raw_ds, float(d))
        push!(norm_ds, dn)

        h1 = Set(hybrid_clade_labels(r1))
        h2 = Set(hybrid_clade_labels(r2))

        correct_hybrids = sort(collect(intersect(h1, h2)))
        net1_only_hybrids = sort(collect(setdiff(h1, h2)))
        net2_only_hybrids = sort(collect(setdiff(h2, h1)))

        push!(matched_blob_hybrid_info, Dict{String,Any}(
            "matched" => true,
            "net1_blob_index" => r1.blob_index,
            "net2_blob_index" => r2.blob_index,
            "boundary_clades" => boundary_key(r1),

            "net1_hybrid_nodes" => r1.hybrid_nodes,
            "net2_hybrid_nodes" => r2.hybrid_nodes,

            "net1_hybrid_clades" => sort(collect(h1)),
            "net2_hybrid_clades" => sort(collect(h2)),

            "correct_hybrid_clades" => correct_hybrids,
            "net1_only_hybrid_clades" => net1_only_hybrids,
            "net2_only_hybrid_clades" => net2_only_hybrids,

            "has_correct_hybrid" => !isempty(correct_hybrids),
            "num_correct_hybrids" => length(correct_hybrids),
            "num_net1_only_hybrids" => length(net1_only_hybrids),
            "num_net2_only_hybrids" => length(net2_only_hybrids),

            "kendall_raw" => float(d),
            "kendall_norm" => dn
        ))
    end

    # Remaining res2 blobs were not matched by any res1 blob.
    unmatched_net2_blobs = Vector{Dict{String,Any}}()
    for remaining in values(dict2)
        for r2 in remaining
            h2 = hybrid_clade_labels(r2)
            push!(unmatched_net2_blobs, Dict{String,Any}(
                "net2_blob_index" => r2.blob_index,
                "boundary_clades" => boundary_key(r2),
                "net2_hybrid_nodes" => r2.hybrid_nodes,
                "net2_hybrid_clades" => h2,
                "num_net2_only_hybrids" => length(h2)
            ))
        end
    end

    avg_norm = isempty(norm_ds) ? NaN : sum(norm_ds) / length(norm_ds)

    # ---- hybrid error rate, treating net1 as the ground truth ----
    #
    # hybrid_unique_net1_only counts hybrid clades present in the
    # true network (net1) but absent from the estimated network (net2).
    #
    # Use the same denominator as the downstream analysis:
    # the number of net1 blobs remaining after the boundary filter.
    #
    #     hybrid_error_rate =
    #         hybrid_unique_net1_only /
    #         net1_blobs_minV5_boundaryFiltered
    #
    # If there are no retained true blobs, record JSON null.
    true_blob_count = length(res1)
    hybrid_error_rate = true_blob_count > 0 ?
        hybrid_unique_net1_only / true_blob_count :
        nothing

    matched_infos = [x for x in matched_blob_hybrid_info if x["matched"] == true]

    summary = Dict{String,Any}(
        "net1_path" => nwk1,
        "net2_path" => nwk2,

        "min_boundary_clades" => min_boundary_clades,

        "net1_blobs_minV5_all" => length(res1_all),
        "net2_blobs_minV5_all" => length(res2_all),
        "net1_blobs_minV5_boundaryFiltered" => length(res1),
        "net2_blobs_minV5_boundaryFiltered" => length(res2),

        "matched_cycles" => matched_cycles,
        "missing_cycles" => missing_cycles,
        "unmatched_net2_cycles" => length(unmatched_net2_blobs),

        "boundary_counts_net1_all" => boundary_counts_net1_all,
        "boundary_counts_net2_all" => boundary_counts_net2_all,
        "boundary_counts_net1" => boundary_counts_net1,
        "boundary_counts_net2" => boundary_counts_net2,

        "boundary_orders_net1" => [r.boundary_order_labels for r in res1],
        "boundary_orders_net2" => [r.boundary_order_labels for r in res2],

        "hybrid_unique_net1_only" => hybrid_unique_net1_only,
        "hybrid_unique_net2_only" => hybrid_unique_net2_only,
        "hybrid_shared_count" => hybrid_shared_count,
        "hybrid_error_rate" => hybrid_error_rate,

        "matched_blob_hybrid_info" => matched_blob_hybrid_info,
        "unmatched_net2_blob_hybrid_info" => unmatched_net2_blobs,

        "matched_blob_has_correct_hybrid" =>
            [x["has_correct_hybrid"] for x in matched_infos],

        "matched_blob_num_correct_hybrids" =>
            [x["num_correct_hybrids"] for x in matched_infos],

        "num_matched_blobs_with_correct_hybrid" =>
            count(x -> x["has_correct_hybrid"], matched_infos),

        "num_matched_blobs_without_correct_hybrid" =>
            count(x -> !x["has_correct_hybrid"], matched_infos),

        "kendall_raw_list" => raw_ds,
        "kendall_norm_list" => norm_ds,
        "avg_kendall_norm" => avg_norm
    )

    if outjson !== nothing
        outpath = write_summary_json(outjson, summary)
        println("Wrote JSON summary to: ", outpath)
    end

    println(
        "SUMMARY: blobs(net1,after_boundary_filter)=",
        summary["net1_blobs_minV5_boundaryFiltered"],
        ", blobs(net2,after_boundary_filter)=",
        summary["net2_blobs_minV5_boundaryFiltered"],
        ", matched_cycles=",
        summary["matched_cycles"],
        ", missing_cycles=",
        summary["missing_cycles"],
        ", matched_blobs_with_correct_hybrid=",
        summary["num_matched_blobs_with_correct_hybrid"],
        ", hybrid_error_rate=",
        summary["hybrid_error_rate"],
        ", avg_kendall_norm=",
        summary["avg_kendall_norm"]
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) == 2 || length(ARGS) == 3)
        error("usage: julia compare_two_networks_by_hybrid_and_circle_ordering.jl net1.nwk net2.nwk [out.json]")
    end

    outjson = (length(ARGS) == 3) ? ARGS[3] : nothing
    compare_two(ARGS[1], ARGS[2]; outjson=outjson, min_boundary_clades=5)
end
