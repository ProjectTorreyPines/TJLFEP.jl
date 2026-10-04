#!/usr/bin/env julia
# Exactness gate for the speedup work: compare a candidate benchmark CSV against a golden one,
# per radius, for one solver column (default :grid). Asserts sfmin/ky/width agree within TOL
# (default rtol 1e-9 == "bitwise" at the CSV's 6-sig-fig precision). Nonzero exit on any mismatch,
# so batch jobs and CI can gate on it.
#
#   julia --project=. dev/ad/check_grid_exact.jl GOLDEN.csv CANDIDATE.csv
# or via env: GOLDEN=... CANDIDATE=... SOLVER=grid TOL=1e-9 FIELDS=sfmin,ky,width
#
# CSV schema (written by benchmark_nls_solvers.jl): ir,solver,sfmin,ky,width,evals,wall_s
using Printf

_arg(i, envkey, default) = length(ARGS) >= i ? ARGS[i] : get(ENV, envkey, default)

const GOLDEN    = _arg(1, "GOLDEN", "")
const CANDIDATE = _arg(2, "CANDIDATE", "")
const SOLVER    = get(ENV, "SOLVER", "grid")
const TOL       = parse(Float64, get(ENV, "TOL", "1e-9"))
const FIELDS    = split(get(ENV, "FIELDS", "sfmin,ky,width"), ",")

if isempty(GOLDEN) || isempty(CANDIDATE)
    error("usage: check_grid_exact.jl GOLDEN.csv CANDIDATE.csv  (or GOLDEN=/CANDIDATE= env)")
end

# Parse a benchmark CSV into Dict{ir => (sfmin,ky,width)} for `solver` rows only.
function _load(path::String, solver::String)
    isfile(path) || error("missing CSV: $path")
    rows = Dict{Int,NamedTuple}()
    open(path) do io
        readline(io)  # header
        for ln in eachline(io)
            isempty(strip(ln)) && continue
            f = split(ln, ",")
            f[2] == solver || continue
            rows[parse(Int, f[1])] = (sfmin=parse(Float64, f[3]), ky=parse(Float64, f[4]),
                                      width=parse(Float64, f[5]))
        end
    end
    return rows
end

_reldiff(a, b) = (a == b) ? 0.0 : abs(a - b) / max(abs(a), abs(b), eps())

function main()
    g = _load(GOLDEN, SOLVER)
    c = _load(CANDIDATE, SOLVER)
    @printf("Exactness check: solver=%s  tol(rtol)=%g\n", SOLVER, TOL)
    @printf("  golden=%s (%d radii)\n  candidate=%s (%d radii)\n", GOLDEN, length(g), CANDIDATE, length(c))

    irs = sort(collect(keys(g)))
    missing_irs = [ir for ir in irs if !haskey(c, ir)]
    extra_irs   = [ir for ir in keys(c) if !haskey(g, ir)]
    nfail = 0
    worst = Dict(f => (0.0, 0) for f in FIELDS)
    for ir in irs
        haskey(c, ir) || continue
        for f in FIELDS
            sym = Symbol(f)
            ga = getproperty(g[ir], sym); ca = getproperty(c[ir], sym)
            rd = _reldiff(ga, ca)
            if rd > worst[f][1]; worst[f] = (rd, ir); end
            if rd > TOL
                nfail += 1
                @printf("  MISMATCH ir=%3d %-6s golden=%.8g candidate=%.8g  reldiff=%.3g\n",
                        ir, f, ga, ca, rd)
            end
        end
    end

    println("\n  worst reldiff per field:")
    for f in FIELDS
        @printf("    %-6s %.3g (ir=%d)\n", f, worst[f][1], worst[f][2])
    end
    if !isempty(missing_irs); @printf("  MISSING radii in candidate: %s\n", join(missing_irs, ",")); end
    if !isempty(extra_irs);   @printf("  EXTRA radii in candidate: %s\n", join(extra_irs, ",")); end

    ok = nfail == 0 && isempty(missing_irs)
    println(ok ? "\nPASS: candidate matches golden within tol." :
                 "\nFAIL: $nfail field mismatch(es)$(isempty(missing_irs) ? "" : " + missing radii").")
    exit(ok ? 0 : 1)
end

main()
