# Harvest SCAN_N=20 compute times (seconds) into timing_runs/scan20_timing.csv

const BUILD = normpath(@__DIR__, "..")
const OUT_CSV = joinpath(BUILD, "timing_runs", "scan20_timing.csv")
const BASIS = [6, 8, 16, 32]

const NB6_FALLBACK = (
    fortran = 62.506,
    julia_cpu = 154.719,
    julia_gpu = 199.813,
)

# Completed comparison/debug jobs without TIMING_RESULT in the log (sacct srun step, seconds).
const KNOWN_SCAN20 = Dict(
    16 => (fortran = 347.0,  julia_cpu = nothing,  julia_gpu = 294.0),   # 53260814.0 (5:47), 53260865 @time, 53260819.0 (4:54)
    32 => (fortran = 1546.0, julia_cpu = 13200.76, julia_gpu = 581.0),   # 53260741.0 (25:46), 53260866 @time, 53260659.0 (9:41)
)

# Report TOTAL WALLCLOCK of the run (phase=total_job): end-to-end wall including process
# spawn, sysimage load, compute, and merge. The Fortran binary compile and the Julia sysimage
# build are the equivalent one-time costs and are EXCLUDED from both (neither is in total_job).
# Falls back to scan/compute for older logs that lack a total_job line.
# Returns (fort, cpu, gpu, cpu_ad, gpu_ad, gpu_truth, gpu_robust) seconds, FOLLOWED BY the matching
# node counts (fort_n, cpu_n, gpu_n, cpu_ad_n, gpu_ad_n, gpu_truth_n, gpu_robust_n) parsed from the
# `nodes=` token of the SAME picked line — so node-hours = nodes × seconds adapts automatically to a
# 5-node wave run vs a 1-node backfill run. The seconds occupy positions 1..7 (existing callers index
# those), node counts occupy 8..14. The `solver` token (grid|ad|robust_ad|truth) distinguishes the
# series; lines without it are treated as solver=grid (back-compat).
function parse_timing_log(path::String)
    if !isfile(path)
        return ntuple(_ -> nothing, 14)
    end
    vals  = Dict{String,Float64}()   # "backend|device|solver|phase" => seconds
    nodes = Dict{String,Float64}()   # same key => node count (from nodes= token)
    for line in readlines(path)
        occursin("TIMING_RESULT", line) || continue
        s = replace(line, " " => "")
        m = match(r"seconds=([0-9.]+)", s)
        m === nothing && continue
        t = parse(Float64, m.captures[1])
        mn = match(r"nodes=([0-9]+)", s)
        nn = mn === nothing ? NaN : parse(Float64, mn.captures[1])
        be  = occursin("backend=fortran", s) ? "fortran" :
              occursin("backend=julia", s)   ? "julia"   : ""
        dev = occursin("device=gpu", s) ? "gpu" :
              occursin("device=cpu", s) ? "cpu" : ""
        # NB: check solver=robust_ad BEFORE solver=ad ("solver=ad" is not a substring of
        # "solver=robust_ad", but keep the explicit ordering for clarity).
        sol = occursin("solver=robust_ad", s) ? "robust_ad" :
              occursin("solver=truth", s)     ? "truth"     :
              occursin("solver=ad", s)        ? "ad"        : "grid"
        ph  = occursin("phase=total_job", s) ? "total"   :
              occursin("phase=scan", s)      ? "scan"    :
              occursin("phase=compute", s)   ? "compute" : ""
        (isempty(be) || isempty(ph)) && continue
        vals["$be|$dev|$sol|$ph"]  = t
        nodes["$be|$dev|$sol|$ph"] = nn
    end
    # pick returns (seconds, nodes) for the first present key, preserving fallback order.
    pick(ks...) = (for k in ks; haskey(vals, k) && return (vals[k], get(nodes, k, NaN)); end; (nothing, nothing))
    fort,    fort_n    = pick("fortran|cpu|grid|total", "fortran|cpu|grid|scan")
    cpu,     cpu_n     = pick("julia|cpu|grid|total", "julia|cpu|grid|compute", "julia|cpu|grid|scan")
    gpu,     gpu_n     = pick("julia|gpu|grid|total", "julia|gpu|grid|scan")
    cpu_ad,  cpu_ad_n  = pick("julia|cpu|ad|total", "julia|cpu|ad|compute", "julia|cpu|ad|scan")
    gpu_ad,  gpu_ad_n  = pick("julia|gpu|ad|total", "julia|gpu|ad|scan")
    gpu_truth, gpu_truth_n = pick("julia|gpu|truth|total", "julia|gpu|truth|scan")
    gpu_robust, gpu_robust_n = pick("julia|gpu|robust_ad|total", "julia|gpu|robust_ad|scan")
    return fort, cpu, gpu, cpu_ad, gpu_ad, gpu_truth, gpu_robust,
           fort_n, cpu_n, gpu_n, cpu_ad_n, gpu_ad_n, gpu_truth_n, gpu_robust_n
end

function parse_ok_in(path::String)
    isfile(path) || return nothing
    for line in reverse(readlines(path))
        m = match(r"OK in ([0-9.]+) s", line)
        m !== nothing && return parse(Float64, m.captures[1])
    end
    return nothing
end

"""@time runTHD pmap wall from debug_compare_*_distributed.jl logs."""
function parse_julia_pmap_compute(path::String)
    isfile(path) || return nothing
    for line in reverse(readlines(path))
        m = match(r"^([0-9.]+) seconds \(", line)
        m !== nothing && return parse(Float64, m.captures[1])
    end
    return nothing
end

"""Max per-task 'OK scan_index=… in X s' from gacode GPU comparison logs."""
function parse_gpu_max_task(path::String)
    isfile(path) || return nothing
    times = Float64[]
    for line in readlines(path)
        m = match(r"OK scan_index=\d+ ir=\d+ sfmin=[^\n]+ in ([0-9.]+) s", line)
        m === nothing && continue
        push!(times, parse(Float64, m.captures[1]))
    end
    return isempty(times) ? nothing : maximum(times)
end

# Trailing SLURM job id from a log filename (…_<jobid>.out), or -1 if absent.
_jobid(f::String) = (m = match(r"_(\d+)\.out$", basename(f)); m === nothing ? -1 : parse(Int, m.captures[1]))

function newest(pattern::String)
    rx = Regex("^" * replace(pattern, "*" => ".*") * "\$")
    files = String[]
    for f in readdir(BUILD)
        occursin(rx, f) || continue
        push!(files, joinpath(BUILD, f))
    end
    isempty(files) && return nothing
    # Rank by (jobid, mtime) descending: the highest SLURM job id is the latest submission,
    # so the most recent re-run always wins. This is robust to mtime corruption from a bulk
    # checkout/touch (which can leave an older contended run with a newer mtime). mtime only
    # breaks ties among logs with no parseable job id.
    return first(sort(files, by=f -> (_jobid(f), mtime(f)), rev=true))
end

function collect_nb(nb::Int)
    fort = cpu = gpu = cpu_ad = gpu_ad = gpu_ad_mps = gpu_truth = gpu_robust = nothing

    # Timing harness logs (preferred -- these are the fresh per-nbasis runs). Take these
    # FIRST so a fresh nb6 run wins over the legacy no-nbasis logs / NB6_FALLBACK below.
    # Grid globs use _[0-9]* so they do NOT match the autodiff _ad_ logs below.
    f_time = newest("time_scan20_nb$(nb)_fortran_[0-9]*.out")
    j_time = newest("time_scan20_nb$(nb)_julia_cpu_[0-9]*.out")
    g_time = newest("time_scan20_nb$(nb)_julia_gpu_[0-9]*.out")
    if fort === nothing && f_time !== nothing
        ft, _, _ = parse_timing_log(f_time)
        fort = ft
    end
    if cpu === nothing && j_time !== nothing
        _, ct, _ = parse_timing_log(j_time)
        cpu = ct
    end
    if gpu === nothing && g_time !== nothing
        _, _, gt = parse_timing_log(g_time)
        gpu = gt
    end

    # Autodiff (solver=:ad) harness logs. The MPS-AD log (_gpu_ad_mps_) carries the SAME
    # device=gpu solver=ad TIMING_RESULT tokens as the baseline AD-GPU log, so they are
    # disambiguated only by filename: the baseline globs use _[0-9]* to exclude _mps_.
    ja_time  = newest("time_scan20_nb$(nb)_julia_cpu_ad_[0-9]*.out")
    ga_time  = newest("time_scan20_nb$(nb)_julia_gpu_ad_[0-9]*.out")
    gam_time = newest("time_scan20_nb$(nb)_julia_gpu_ad_mps_*.out")
    if ja_time !== nothing
        cpu_ad = parse_timing_log(ja_time)[4]
    end
    if ga_time !== nothing
        gpu_ad = parse_timing_log(ga_time)[5]
    end
    if gam_time !== nothing
        gpu_ad_mps = parse_timing_log(gam_time)[5]
    end

    # Physical-truth (solver=:truth) GPU+MPS harness logs (device=gpu solver=truth tokens).
    gt_time = newest("time_scan20_nb$(nb)_julia_gpu_truth_[0-9]*.out")
    if gt_time !== nothing
        gpu_truth = parse_timing_log(gt_time)[6]
    end

    # Robust_ad (solver=:robust_ad, width-extended) GPU+MPS harness logs (device=gpu
    # solver=robust_ad tokens) -- the WIDTH tier (no nbasis ladder).
    gr_time = newest("time_scan20_nb$(nb)_julia_gpu_robust_ad_[0-9]*.out")
    if gr_time !== nothing
        gpu_robust = parse_timing_log(gr_time)[7]
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
        pmap_t = parse_julia_pmap_compute(j_dbg)
        ok = parse_ok_in(j_dbg)
        if ct !== nothing
            cpu = ct
        elseif pmap_t !== nothing
            cpu = pmap_t
        elseif ok !== nothing
            cpu = ok
        end
    end
    if gpu === nothing && g_dbg !== nothing
        _, _, gt = parse_timing_log(g_dbg)
        gpu = gt
    end

    # KNOWN_SCAN20 are OLD pre-optimization/pre-sysimage numbers: use them only as a
    # FALLBACK when no fresh log was found, never to override a freshly measured run.
    if haskey(KNOWN_SCAN20, nb)
        k = KNOWN_SCAN20[nb]
        if k.fortran !== nothing
            fort = fort === nothing ? k.fortran : fort
        end
        if k.julia_cpu !== nothing
            cpu = cpu === nothing ? k.julia_cpu : cpu
        end
        if k.julia_gpu !== nothing
            gpu = gpu === nothing ? k.julia_gpu : gpu
        end
    end
    if gpu === nothing && g_dbg !== nothing
        gpu = parse_gpu_max_task(g_dbg)
    end

    # nb6 legacy fallbacks: older runs used logs without the "nb6" tag. Only consulted if
    # no fresh nb6-tagged log was found above.
    if nb == 6
        if fort === nothing
            lf = newest("time_scan20_fortran_*.out")
            lf !== nothing && (fort = parse_timing_log(lf)[1])
        end
        if cpu === nothing
            lc = newest("time_scan20_julia_cpu_*.out")
            lc !== nothing && (cpu = parse_timing_log(lc)[2])
        end
        if gpu === nothing
            lg = newest("time_scan20_julia_gpu_*.out")
            lg !== nothing && (gpu = parse_timing_log(lg)[3])
        end
        fort = fort === nothing ? NB6_FALLBACK.fortran : fort
        cpu = cpu === nothing ? NB6_FALLBACK.julia_cpu : cpu
        gpu = gpu === nothing ? NB6_FALLBACK.julia_gpu : gpu
    end

    # ── Node counts (for node-hours = nodes × seconds) ───────────────────────────
    # Parse the nodes= token of the SAME chosen log so a 1-node backfill run reports nodes=1 and a
    # 5/10-node wave run reports its real count. Fall back to the fixed batch layout when the value
    # came from a constant/legacy/debug log without a nodes= token. node idx in parse_timing_log:
    # 8=fort 9=cpu 10=gpu 11=cpu_ad 12=gpu_ad 13=gpu_truth 14=gpu_robust (gpu_ad_mps shares idx 12).
    DEFAULT_N = (fortran=10.0, julia_cpu=10.0, julia_gpu=5.0,
                 cpu_ad=10.0, gpu_ad=5.0, gpu_ad_mps=5.0, gpu_truth=5.0, gpu_robust=5.0)
    node_or(default, log, idx) = begin
        log === nothing && return default
        v = parse_timing_log(log)[idx]
        (v === nothing || (v isa Float64 && isnan(v))) ? default : v
    end
    fort_n       = node_or(DEFAULT_N.fortran,   f_time !== nothing ? f_time : f_dbg, 8)
    cpu_n        = node_or(DEFAULT_N.julia_cpu,  j_time !== nothing ? j_time : j_dbg, 9)
    gpu_n        = node_or(DEFAULT_N.julia_gpu,  g_time !== nothing ? g_time : g_dbg, 10)
    cpu_ad_n     = node_or(DEFAULT_N.cpu_ad,     ja_time,  11)
    gpu_ad_n     = node_or(DEFAULT_N.gpu_ad,     ga_time,  12)
    gpu_ad_mps_n = node_or(DEFAULT_N.gpu_ad_mps, gam_time, 12)
    gpu_truth_n  = node_or(DEFAULT_N.gpu_truth,  gt_time,  13)
    gpu_robust_n = node_or(DEFAULT_N.gpu_robust, gr_time,  14)

    return fort, cpu, gpu, cpu_ad, gpu_ad, gpu_ad_mps, gpu_truth, gpu_robust,
           fort_n, cpu_n, gpu_n, cpu_ad_n, gpu_ad_n, gpu_ad_mps_n, gpu_truth_n, gpu_robust_n
end

function collect_scan20_timing!()
    # Schema: 8 wallclock-seconds columns, then 8 node-hours columns (nodes × seconds / 3600), then
    # notes (always last). node-hours is the resource-cost metric for the node-hours-vs-nbasis plot.
    rows = ["n_basis," *
            "fortran_s,julia_cpu_s,julia_gpu_s,julia_cpu_ad_s,julia_gpu_ad_s,julia_gpu_ad_mps_s,julia_gpu_truth_s,julia_gpu_robust_ad_s," *
            "fortran_nh,julia_cpu_nh,julia_gpu_nh,julia_cpu_ad_nh,julia_gpu_ad_nh,julia_gpu_ad_mps_nh,julia_gpu_truth_nh,julia_gpu_robust_ad_nh," *
            "notes"]
    for nb in BASIS
        fort, cpu, gpu, cpu_ad, gpu_ad, gpu_ad_mps, gpu_truth, gpu_robust,
            fort_n, cpu_n, gpu_n, cpu_ad_n, gpu_ad_n, gpu_ad_mps_n, gpu_truth_n, gpu_robust_n = collect_nb(nb)
        notes = String[]
        fort === nothing && push!(notes, "fortran_missing")
        cpu === nothing && push!(notes, "julia_cpu_missing")
        gpu === nothing && push!(notes, "julia_gpu_missing")
        cpu_ad === nothing && push!(notes, "julia_cpu_ad_missing")
        gpu_ad === nothing && push!(notes, "julia_gpu_ad_missing")
        # gpu_ad_mps is Phase B (pending the GPU/MPS AD kernel); absence is expected, not flagged.
        note = isempty(notes) ? "ok" : join(notes, ";")
        cell(x) = x === nothing ? "" : string(x)
        # node-hours = nodes × seconds / 3600 (blank when seconds missing).
        nh(sec, n) = (sec === nothing || n === nothing) ? "" : string(sec * n / 3600.0)
        push!(rows, join((nb,
            cell(fort), cell(cpu), cell(gpu), cell(cpu_ad), cell(gpu_ad), cell(gpu_ad_mps), cell(gpu_truth), cell(gpu_robust),
            nh(fort, fort_n), nh(cpu, cpu_n), nh(gpu, gpu_n), nh(cpu_ad, cpu_ad_n), nh(gpu_ad, gpu_ad_n),
            nh(gpu_ad_mps, gpu_ad_mps_n), nh(gpu_truth, gpu_truth_n), nh(gpu_robust, gpu_robust_n),
            note), ","))
    end
    mkpath(dirname(OUT_CSV))
    write(OUT_CSV, join(rows, "\n") * "\n")
    return OUT_CSV, rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Pkg
    Pkg.activate(normpath(@__DIR__, "..", ".."))
    path, rows = collect_scan20_timing!()
    println("Wrote ", path)
    for line in rows[2:end]
        println("  ", line)
    end
end
