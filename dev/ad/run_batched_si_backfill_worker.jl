using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLFEP, TJLF, Printf
using TJLFEP: preprocess_gacode_inputs, kwscale_scan

# One backfill worker PROCESS pinned to a single GPU (via CUDA_VISIBLE_DEVICES). Workers drain a
# shared radius queue by atomic-mkdir work-stealing on QDIR: worker that first creates QDIR/c<i>
# owns radius i. This is process-based (unlike the earlier in-process Threads.@spawn version, which
# crashed because kwscale_scan's dense endpoints call Threads.@threads and nesting @threads inside
# @spawn on the shared pool is illegal). Each process therefore owns its own -t thread pool and its
# own GPU — identical resourcing to the working fixed-shard runner (run_batched_si_sfmin.jl), so the
# node-hours are directly comparable; only the radius->GPU assignment differs (dynamic vs static).
CASE = get(ENV, "CASE_DIR", normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1"))
NB   = parse(Int, get(ENV, "NB", "32"))
GAC  = joinpath(CASE, "input.gacode")
TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
QDIR = ENV["QDIR"]                          # shared claim dir (must already exist)
WID  = parse(Int, get(ENV, "WID", "0"))     # worker id, for logging
OUT  = ENV["OUT"]                           # this worker's result shard: "idx IR sfmin secs wid"
grid = (nfactor = parse(Int, get(ENV,"NFACTOR","8")), nefwid = parse(Int, get(ENV,"NEFWID","8")),
        nkyhat  = parse(Int, get(ENV,"NKYHAT","4")), k_max  = parse(Int, get(ENV,"KMAX","4")))

base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
scan_n = Int(base_ep.SCAN_N)
dev = CUDA.device()
@printf("worker %d: gpu=%s nb=%d radii=%d grid=%s\n", WID, repr(dev), NB, scan_n, grid); flush(stdout)

open(OUT, "w") do io
    for i in 1:scan_n
        claimed = try
            mkdir(joinpath(QDIR, "c$(i)")); true      # atomic test-and-set on the shared FS
        catch
            false
        end
        claimed || continue
        ep = deepcopy(base_ep); ep.IR = base_ep.IR_EXP[i]
        ep.WIDTH_IN_FLAG = false; ep.MODE_IN = 2; ep.KY_MODEL = 3; ep.PROCESS_IN = 5
        ep.FACTOR_IN = Float64(base_ep.FACTOR[i])
        t = @elapsed begin
            _g, epo, = kwscale_scan(ep, prof, false; use_gpu=true, inner=:batched_si, grid...)
        end
        sf = Float64(epo.FACTOR_IN)
        @printf(io, "%d %d %.16g %.1f %d\n", i, base_ep.IR_EXP[i], sf, t, WID); flush(io)
        @printf("  [w%d gpu=%s] ir=%3d sfmin=%9.4g %5.0fs\n", WID, repr(dev), base_ep.IR_EXP[i], sf, t)
        flush(stdout)
    end
end
@printf("worker %d done\n", WID); flush(stdout)
