# Harvest SCAN_N=20 compute times (seconds) into timing_runs/scan20_timing.csv

const BUILD = @__DIR__
const OUT_CSV = joinpath(BUILD, "timing_runs", "scan20_timing.csv")
const BASIS = [6, 8, 16, 32]

const NB6_FALLBACK = (
    fortran = 62.506,
    julia_cpu = 154.719,
    julia_gpu = 199.813,
)

function parse_timing_log(path::String)
    fort = cpu = gpu = nothing
    isfile(path) || return fort, cpu, gpu
    for line in readlines(path)
        occursin("TIMING_RESULT", line) || continue
        m = match(r"seconds=([0-9.]+)", replace(line, " " => ""))
        m === nothing && continue
        t = parse(Float64, m.captures[1])
        if occursin("backend=fortran", line) && occursin("phase=scan", line)
            fort = t
        elseif occursin("backend=julia", line) && occursin("device=cpu", line) &&
               (occursin("phase=compute", line) || occursin("phase=scan", line))
            cpu = t
        elseif occursin("backend=julia", line) && occursin("device=gpu", line) &&
               occursin("phase=scan", line)
            gpu = t
        end
    end
    return fort, cpu, gpu
end

function parse_ok_in(path::String)
    isfile(path) || return nothing
    for line in reverse(readlines(path))
        m = match(r"OK in ([0-9.]+) s", line)
        m !== nothing && return parse(Float64, m.captures[1])
    end
    return nothing
end

function newest(pattern::String)
    rx = Regex("^" * replace(pattern, "*" => ".*") * "\$")
    files = String[]
    for f in readdir(BUILD)
        occursin(rx, f) || continue
        push!(files, joinpath(BUILD, f))
    end
    isempty(files) && return nothing
    return first(sort(files, by=f -> mtime(f), rev=true))
end

function collect_nb(nb::Int)
    fort = cpu = gpu = nothing

    if nb == 6
        f = newest("time_scan20_fortran_*.out")
        fort, _, _ = f === nothing ? (nothing, nothing, nothing) : parse_timing_log(f)
        cpath = newest("time_scan20_julia_cpu_*.out")
        _, c, _ = cpath === nothing ? (nothing, nothing, nothing) : parse_timing_log(cpath)
        gpath = newest("time_scan20_julia_gpu_*.out")
        _, _, g = gpath === nothing ? (nothing, nothing, nothing) : parse_timing_log(gpath)
        fort = something(fort, NB6_FALLBACK.fortran)
        cpu = something(c, NB6_FALLBACK.julia_cpu)
        gpu = something(g, NB6_FALLBACK.julia_gpu)
    end

    # Timing harness logs (preferred for timing runs)
    j_time = newest("time_scan20_nb$(nb)_julia_cpu_*.out")
    g_time = newest("time_scan20_nb$(nb)_julia_gpu_*.out")
    if cpu === nothing && j_time !== nothing
        _, ct, _ = parse_timing_log(j_time)
        cpu = ct
    end
    if gpu === nothing && g_time !== nothing
        _, _, gt = parse_timing_log(g_time)
        gpu = gt
    end

    f_dbg = newest("debug_nb$(nb)_fortran20_10n_*.out")
    j_dbg = newest("debug_nb$(nb)_julia20_10n_*.out")
    g_dbg = newest("gacode_nb$(nb)_scan20_gpu5_*.out")

    if fort === nothing && f_dbg !== nothing
        ft, _, _ = parse_timing_log(f_dbg)
        fort = ft
    end
    if cpu === nothing && j_dbg !== nothing
        _, ct, _ = parse_timing_log(j_dbg)
        ok = parse_ok_in(j_dbg)
        cpu = ct !== nothing ? ct : ok
    end
    if gpu === nothing && g_dbg !== nothing
        _, _, gt = parse_timing_log(g_dbg)
        gpu = gt
    end

    return fort, cpu, gpu
end

function collect_scan20_timing!()
    rows = ["n_basis,fortran_s,julia_cpu_s,julia_gpu_s,notes"]
    for nb in BASIS
        fort, cpu, gpu = collect_nb(nb)
        notes = String[]
        fort === nothing && push!(notes, "fortran_missing")
        cpu === nothing && push!(notes, "julia_cpu_missing")
        gpu === nothing && push!(notes, "julia_gpu_missing")
        note = isempty(notes) ? "ok" : join(notes, ";")
        fs = fort === nothing ? "" : string(fort)
        cs = cpu === nothing ? "" : string(cpu)
        gs = gpu === nothing ? "" : string(gpu)
        push!(rows, "$nb,$fs,$cs,$gs,$note")
    end
    mkpath(dirname(OUT_CSV))
    write(OUT_CSV, join(rows, "\n") * "\n")
    return OUT_CSV, rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(normpath(@__DIR__, ".."))
    path, rows = collect_scan20_timing!()
    println("Wrote ", path)
    for line in rows[2:end]
        println("  ", line)
    end
end
