#!/usr/bin/env julia

using CSV
using DataFrames
using PhyloNetworks

function read_network_from_file(path::String)
    s = strip(read(path, String))
    return PhyloNetworks.readnewick(s)
end

function count_reticulations_from_file(network_file::String)
    net = read_network_from_file(network_file)

    # Try PhyloNetworks-based count first
    try
        return size(PhyloNetworks.minorreticulationmatrix(net), 1)
    catch
        # Fallback: count unique #H labels in the text
        s = strip(read(network_file, String))
        labels = Set{String}()
        for m in eachmatch(r"#H\d+", s)
            push!(labels, m.match)
        end
        return length(labels)
    end
end

# Clean a header string: remove BOM and surrounding whitespace
function clean_header(x)
    s = String(x)
    s = replace(s, '\ufeff' => "")   # remove UTF-8 BOM if present
    s = strip(s)
    return s
end

function find_column(df::DataFrame, wanted::String)
    cleaned = [clean_header(string(n)) for n in names(df)]
    idx = findfirst(==(wanted), cleaned)
    return idx
end

function choose_row(df::DataFrame, target::Int)
    idx_num = find_column(df, "Number of Branches")
    idx_nwk = find_column(df, "Extended Newick")

    println("Detected columns: ", [repr(clean_header(string(n))) for n in names(df)])

    if idx_num === nothing
        error("Missing required column: Number of Branches")
    end
    if idx_nwk === nothing
        error("Missing required column: Extended Newick")
    end

    numcol = names(df)[idx_num]
    nwkcol = names(df)[idx_nwk]

    # rename to stable internal names
    rename!(df, numcol => :Number_of_Branches, nwkcol => :Extended_Newick)

    sort!(df, :Number_of_Branches)

    vals = df[!, :Number_of_Branches]
    idx = findfirst(==(target), vals)

    if idx !== nothing
        return df[idx, :], "exact_match"
    end

    # your requested rule: if no exact match, use the last line
    return df[nrow(df), :], maximum(vals) < target ?
        "fallback_last_line_max_less_than_target" :
        "fallback_last_line_no_exact_match"
end

function main()
    if length(ARGS) != 3
        println(stderr, "Usage: julia select_network_by_hybrid_num.jl input_network.nwk input.csv output.nwk")
        exit(1)
    end

    network_file = ARGS[1]
    csv_file = ARGS[2]
    output_file = ARGS[3]

    target = count_reticulations_from_file(network_file)
    println("Reticulation count in input network: ", target)

    df = CSV.read(csv_file, DataFrame)
    row, mode = choose_row(df, target)

    chosen_newick = String(row[:Extended_Newick])
    chosen_branches = row[:Number_of_Branches]

    open(output_file, "w") do io
        println(io, chosen_newick)
    end

    println("Selection mode: ", mode)
    println("Chosen Number of Branches: ", chosen_branches)
    println("Wrote selected network to: ", output_file)
end

main()