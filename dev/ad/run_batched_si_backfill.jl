using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
import CUDA
using TJLFEP, TJLF, Printf
using TJLFEP: preprocess_gacode_inputs, kwscale_scan

# Node-hours of the 20-radius inner=:batched_si scan under a DYNAMIC BACKFILL layout: NGPU workers
# on one node, each pinned to a GPU, pull the next radius index from a shared queue until drained
# (work-stealing). This keeps every GPU busy despite the non-uniform per-radius cost (the dense
# fallback / edge-radius stragglers), unlike a fixed round-robin shard whose wall = slowest shard.
# node-hours = (#nodes = 1) * wall_seconds / 3600. Compare to the fixed-shard number and to
# grid/ad in docs/plots/scan20_timing.csv.
CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
NB   = parse(Int, get(ENV, "NB", "32"))
GAC  = joinpath(CASE, "input.gacode")
TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
NGPU = min(parse(Int, get(ENV, "NGPU", "4")), length(CUDA.devices()))
grid = (nfactor = parse(Int, get(ENV,"NFACTOR","8")), nefwid = parse(Int, get(ENV,"NEFWID","8")),
        nkyhat  = parse(Int, get(ENV,"NKYHAT","4")), k_max  = parse(Int, get(ENV,"KMAX","4")))

base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
scan_n = Int(base_ep.SCAN_N)
@printf("backfill batched_si: nb=%d ngpu=%d radii=%d grid=%s\n", NB, NGPU, scan_n, grid)

sfmin = Vector{Float64}(undef, scan_n)
rtime = Vector{Float64}(undef, scan_n)
rgpu  = Vector{Int}(undef, scan_n)
q = Channel{Int}(scan_n); for i in 1:scan_n; put!(q, i); end; close(q)

twall = @elapsed begin
    @sync for g in 1:NGPU
        Threads.@spawn begin
            CUDA.device!(g - 1)
            for i in q                      # work-stealing: drains until channel empty
                ep = deepcopy(base_ep); ep.IR = base_ep.IR_EXP[i]
                ep.WIDTH_IN_FLAG = false; ep.MODE_IN = 2; ep.KY_MODEL = 3; ep.PROCESS_IN = 5
                ep.FACTOR_IN = Float64(base_ep.FACTOR[i])
                t = @elapsed begin
                    _gg, epo, = kwscale_scan(ep, prof, false; use_gpu=true, inner=:batched_si, grid...)
                end
                sfmin[i] = Float64(epo.FACTOR_IN); rtime[i] = t; rgpu[i] = g - 1
                @printf("  [gpu%d] ir=%3d sfmin=%9.4g %5.0fs\n", g-1, base_ep.IR_EXP[i], sfmin[i], t); flush(stdout)
            end
        end
    end
end

OUT = get(ENV, "OUT", joinpath(@__DIR__, "batched_si_sfmin_nb$(NB).txt"))
open(OUT, "w") do io
    for i in 1:scan_n; @printf(io, "%d %d %.16g\n", i, base_ep.IR_EXP[i], sfmin[i]); end
end
srt = sort(rtime; rev=true)
@printf("\nnb=%d BACKFILL: wall=%.1fs  node_hours=%.4f  (ngpu=%d)\n", NB, twall, twall/3600, NGPU)
@printf("  per-radius time: min=%.0fs median=%.0fs max=%.0fs  (top3 stragglers: %.0f/%.0f/%.0f s)\n",
        minimum(rtime), srt[cld(end,2)], maximum(rtime), srt[1], srt[2], srt[3])
@printf("  sum(per-radius)=%.0fs  ideal_1gpu_wall=%.0fs  backfill_wall=%.0fs  parallel_eff=%.2f\n",
        sum(rtime), sum(rtime)/NGPU, twall, (sum(rtime)/NGPU)/twall)
csv = get(ENV, "CSV_OUT", "")
isempty(csv) || open(csv, "a") do io
    @printf(io, "%d,backfill,%d,%.1f,%.4f\n", NB, NGPU, twall, twall/3600)
end
