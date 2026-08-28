#!/usr/bin/env julia

using PhyloNetworks
using ArgParse
using CSV
using DataFrames

function read_network(path::String)
    s = strip(read(path, String))
    return readnewick(s) |> directedges!
end

function leaf_nodes(net)
    return [n for n in net.node if n.leaf]
end

function leaf_names(net)
    return sort([n.name for n in leaf_nodes(net)])
end

function hybrid_above_bits(net, taxa_order::Vector{String})
    leaves_by_name = Dict(n.name => n for n in leaf_nodes(net))
    hybrid_nodes = [n for n in net.node if n.hybrid]

    bits = Int[]
    for taxon in taxa_order
        leaf = leaves_by_name[taxon]
        has_hybrid_ancestor = any(h -> PhyloNetworks.isdescendant(leaf, h), hybrid_nodes)
        push!(bits, has_hybrid_ancestor ? 1 : 0)
    end
    return bits
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--net1"
            help = "First/reference network file in extended Newick format"
            required = true
        "--net2"
            help = "Second/estimated network file in extended Newick format"
            required = true
        "--out"
            help = "Output CSV file"
            default = "hybrid_leaf_bitvectors.csv"
    end

    args = parse_args(s)

    net1 = read_network(args["net1"])
    net2 = read_network(args["net2"])

    taxa1 = leaf_names(net1)
    taxa2 = leaf_names(net2)

    if taxa1 != taxa2
        error("The two networks do not have the same sorted leaf set.")
    end

    taxa = taxa1
    bits1 = hybrid_above_bits(net1, taxa)
    bits2 = hybrid_above_bits(net2, taxa)

    fp = sum((bits1 .== 0) .& (bits2 .== 1))
    fn = sum((bits1 .== 1) .& (bits2 .== 0))

    n1_ones = sum(bits1)
    n2_ones = sum(bits2)

    fn_rate = n1_ones == 0 ? 0.0 : fn / n1_ones
    fp_rate = n2_ones == 0 ? 0.0 : fp / n2_ones

    df = DataFrame(
        taxon = taxa,
        net1_bit = bits1,
        net2_bit = bits2,
        status = [
            bits1[i] == 0 && bits2[i] == 1 ? "FP" :
            bits1[i] == 1 && bits2[i] == 0 ? "FN" :
            "match"
            for i in eachindex(taxa)
        ],
    )

    summary_row = DataFrame(
        taxon = ["SUMMARY"],
        net1_bit = [n1_ones],
        net2_bit = [n2_ones],
        status = ["FN=$(fn),FP=$(fp),FN_rate=$(fn_rate),FP_rate=$(fp_rate)"]
    )

    df_out = vcat(df, summary_row)

    CSV.write(args["out"], df_out)

    println("Taxa order:")
    println(join(taxa, ","))

    println("\nNetwork 1 bit vector:")
    println(join(bits1, ""))

    println("\nNetwork 2 bit vector:")
    println(join(bits2, ""))

    println("\n# leaves with 1 in network 1 = ", n1_ones)
    println("# leaves with 1 in network 2 = ", n2_ones)

    println("\nFP = ", fp)
    println("FN = ", fn)
    println("FP rate = ", fp_rate)
    println("FN rate = ", fn_rate)

    println("\nWrote: ", args["out"])
end

main()