#!/usr/bin/env julia

using PhyloNetworks

function edge_nodes(e)
    return e.node[1], e.node[2]
end

function node_id(n)
    return n.number
end

function blob_degree(blob_edges, all_edges)
    blob_nodes = Set{Int}()

    for e in blob_edges
        n1, n2 = edge_nodes(e)
        push!(blob_nodes, node_id(n1))
        push!(blob_nodes, node_id(n2))
    end

    deg = 0
    incident_cut_edges = []

    for e in all_edges
        n1, n2 = edge_nodes(e)
        id1 = node_id(n1)
        id2 = node_id(n2)

        in1 = id1 in blob_nodes
        in2 = id2 in blob_nodes

        # exactly one endpoint is inside the blob
        if xor(in1, in2)
            deg += 1
            push!(incident_cut_edges, e)
        end
    end

    return deg, blob_nodes, incident_cut_edges
end

function find_degree4_blobs(net)
    bccs = biconnectedcomponents(net)

    degree4 = []

    for (i, blob_edges) in enumerate(bccs)

        # skip trivial blobs: a single cut edge
        if length(blob_edges) <= 1
            continue
        end

        deg, blob_nodes, cut_edges = blob_degree(blob_edges, net.edge)

        if deg == 4
            push!(degree4, (i, blob_nodes, cut_edges))
        end
    end

    return degree4
end

function main()
    if length(ARGS) != 1
        println("Usage: julia check_four_blob.jl network.nwk")
        exit(1)
    end

    netfile = ARGS[1]

    net = readnewick(strip(read(netfile, String)))
    directedges!(net)

    blobs = find_degree4_blobs(net)

    if isempty(blobs)
        println("No degree-4 blob found.")
    else
        println("Contains degree-4 blob(s): YES")

        for (idx, nodes, cut_edges) in blobs
            println("Blob component index: ", idx)
            println("  nodes: ", sort(collect(nodes)))
            println("  degree: ", length(cut_edges))
        end
    end
end

main()