#!/usr/bin/env julia
# Compare two networks using:
#  - Hardwired Cluster Distance (HWCD)
#  - Normalized HWCD 
#
# Output (CSV, one line, no header):
#   n_leaves,edges_net1,edges_net2,hwcd,hwcd_norm
#
# Usage:
#   julia compare_hwcd_paper_norm.jl net1.nwk net2.nwk [rooted=true|false]

using PhyloNetworks

# ---------------- utilities ----------------

function parse_bool(s::AbstractString)
    ls = lowercase(strip(s))
    if ls in ("true","t","1","yes","y")
        return true
    elseif ls in ("false","f","0","no","n")
        return false
    else
        error("Could not parse boolean from '$s' (expected true/false)")
    end
end

# number of leaves |V_L|
nleaves(net::HybridNetwork) =
    count(n -> getfield(n, :leaf), net.node)

# denominator term: |E| - |V_L|
denom_term(net::HybridNetwork) =
    length(net.edge) - nleaves(net)

# ---------------- main ----------------

function main()
    if length(ARGS) < 2 || length(ARGS) > 3
        println(stderr,
            "Usage: julia compare_hwcd_paper_norm.jl net1.nwk net2.nwk [rooted=true|false]")
        exit(2)
    end

    net1_path = ARGS[1]
    net2_path = ARGS[2]
    rooted = (length(ARGS) == 3) ? parse_bool(ARGS[3]) : false

    net1 = readnewick(net1_path)
    net2 = readnewick(net2_path)

    # ensure consistent edge directions
    directEdges!(net1)
    directEdges!(net2)

    # counts
    nL = nleaves(net1)              # assumed same taxa set
    e1 = length(net1.edge)
    e2 = length(net2.edge)

    # HWCD
    hwcd = hardwiredclusterdistance(net1, net2, rooted)

    # normalized HWCD (paper definition)
    denom = denom_term(net1) + denom_term(net2)
    hwcd_norm = denom > 0 ? hwcd / denom : 0.0

    # CSV output, no header (bash-friendly)
    println(string(nL), ",",
            string(e1), ",",
            string(e2), ",",
            string(hwcd), ",",
            string(hwcd_norm))
end

main()
