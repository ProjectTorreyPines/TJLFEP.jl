#!/usr/bin/env julia
# Harvest the UCP_complete SCAN_N=20 timing sweep into timing_runs/ucp_scan20_timing.csv.
# Keys each series by its ucp-tagged log filename (time_scan20_ucp_nb${nb}_<label>_<jobid>.out)
# and reads the phase=total_job (fallback scan/compute) `seconds=` and `nodes=` tokens.
# node-hours = nodes * seconds / 3600. Independent of the DIII-D collector (no fallbacks).
#   cd dev && julia --project=.. timing/collect_ucp_timing.jl

using Pkg
Pkg.activate(normpath(@__DIR__, "..", ".."))

const BUILD = normpath(@__DIR__, "..")
const OUT_CSV = joinpath(BUILD, "timing_runs", "ucp_scan20_timing.csv")
const BASIS = [6, 8, 16, 32, 40, 48]

# label => default node count (used only if a log lacks a nodes= token)
const SERIES = [
    ("fortran",             "fortran",             10),
    ("julia_gpu",           "julia_gpu",            5),
    ("julia_gpu_ad_only",   "julia_gpu_ad_only",    5),
    ("julia_gpu_ad_locate", "julia_gpu_ad_locate",  1),
    ("julia_gpu_ad_wide",   "julia_gpu_ad_wide",    1),
    ("julia_cpu",           "julia_cpu",           10),
]

_jobid(f) = (m = match(r"_(\d+)\.out$", basename(f)); m === nothing ? -1 : parse(Int, m.captures[1]))

# newest log for a given nb + label (highest job id wins).
function newest(nb::Int, label::String)
    rx = Regex("^time_scan20_ucp_nb$(nb)_$(label)_\\d+\\.out\$")
    files = [joinpath(BUILD, f) for f in readdir(BUILD) if occursin(rx, f)]
    isempty(files) && return nothing
    return first(sort(files; by=f -> (_jobid(f), mtime(f)), rev=true))
end

# (seconds, nodes) from the best TIMING_RESULT phase in a log, or (nothing, nothing).
function parse_log(path)
    path === nothing && return (nothing, nothing)
    isfile(path) || return (nothing, nothing)
    best = Dict{String,Tuple{Float64,Float64}}()   # phase => (sec, nodes)
    for line in readlines(path)
        occursin("TIMING_RESULT", line) || continue
        s = replace(line, " " => "")
        m = match(r"seconds=([0-9.]+)", s); m === nothing && continue
        sec = parse(Float64, m.captures[1])
        mn = match(r"nodes=([0-9]+)", s); nn = mn === nothing ? NaN : parse(Float64, mn.captures[1])
        ph = occursin("phase=total_job", s) ? "total" :
             occursin("phase=scan", s)      ? "scan"  :
             occursin("phase=compute", s)   ? "compute" : ""
        isempty(ph) && continue
        best[ph] = (sec, nn)
    end
    for ph in ("total", "scan", "compute")
        haskey(best, ph) && return best[ph]
    end
    return (nothing, nothing)
end

function main()
    hdr = String["n_basis"]
    for (label, _, _) in SERIES
        push!(hdr, "$(label)_s")
    end
    for (label, _, _) in SERIES
        push!(hdr, "$(label)_nh")
    end
    push!(hdr, "notes")
    rows = [join(hdr, ",")]

    println(rpad("nb", 5), join([rpad(l, 22) for (l, _, _) in SERIES]))
    for nb in BASIS
        secs = Union{Nothing,Float64}[]
        nhs  = Union{Nothing,Float64}[]
        missinglist = String[]
        for (label, glob, defn) in SERIES
            sec, nodes = parse_log(newest(nb, glob))
            n = (nodes === nothing || (nodes isa Float64 && isnan(nodes))) ? Float64(defn) : nodes
            push!(secs, sec)
            push!(nhs, sec === nothing ? nothing : sec * n / 3600.0)
            sec === nothing && push!(missinglist, label)
        end
        cell(x) = x === nothing ? "" : string(round(x; digits=3))
        note = isempty(missinglist) ? "ok" : "missing:" * join(missinglist, ";")
        push!(rows, join(vcat(string(nb), cell.(secs), cell.(nhs), note), ","))
        println(rpad(nb, 5), join([rpad(secs[i] === nothing ? "-" : string(round(secs[i]; digits=1)), 22) for i in eachindex(SERIES)]))
    end

    mkpath(dirname(OUT_CSV))
    write(OUT_CSV, join(rows, "\n") * "\n")
    println("\nWrote ", OUT_CSV)
end

main()
