using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
using TJLFEP, TJLF
using TJLFEP: preprocess_gacode_inputs, kwscale_scan

# Harvest the exact (A,B) grid pencils for ONE radius by running its kwscale grid with the dense
# eigensolver (inner=:threads, use_gpu=false) and TJLF_DUMP_PENCILS set. These are the same pencils
# inner=:batched_si would route through its GPU solver, so running geev-vs-SI on them isolates
# eigensolver accuracy from the kwscale collect/replay wiring.
CASE = normpath(@__DIR__, "..", "..", "examples", "DIIID_202017C42_500ms_v3.1")
NB   = parse(Int, get(ENV, "NB", "16"))
IR   = parse(Int, get(ENV, "IR", "101"))      # experimental radius label (base_ep.IR_EXP value)
GAC  = joinpath(CASE, "input.gacode")
TGL  = joinpath(CASE, "input_scan20_nb$(NB).TGLFEP")
grid = (nfactor = parse(Int, get(ENV,"NFACTOR","8")),
        nefwid  = parse(Int, get(ENV,"NEFWID","8")),
        nkyhat  = parse(Int, get(ENV,"NKYHAT","4")),
        k_max   = parse(Int, get(ENV,"KMAX","4")))

base_ep, prof, _ = preprocess_gacode_inputs(GAC, TGL)
i = findfirst(==(IR), Int.(base_ep.IR_EXP))
i === nothing && error("IR=$IR not in scan (IR_EXP=$(Int.(base_ep.IR_EXP)))")

ep = deepcopy(base_ep); ep.IR = base_ep.IR_EXP[i]
ep.WIDTH_IN_FLAG = false; ep.MODE_IN = 2; ep.KY_MODEL = 3; ep.PROCESS_IN = 5
ep.FACTOR_IN = Float64(base_ep.FACTOR[i])
@info "harvesting" IR NB grid factor_in=ep.FACTOR_IN dump=get(ENV,"TJLF_DUMP_PENCILS","")

g, epo, = kwscale_scan(ep, prof, false; use_gpu=false, inner=:threads, grid...)
@info "harvest done" sfmin=Float64(epo.FACTOR_IN) ky=Float64(epo.KYMARK) w=Float64(epo.WIDTH_IN)
