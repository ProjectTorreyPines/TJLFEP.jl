# Per-radius parity check between two scan logs that each contain a `SFmin = [...]` line.
# Usage:  julia --project=. timing/check_sfmin_parity.jl <baseline.out> <new.out> [tol]
# Exits nonzero if the max per-radius relative difference exceeds `tol` (default 1e-6).

using Printf

function read_sfmin(path::AbstractString)
    line = nothing
    for l in eachline(path)
        if occursin("SFmin", l) && occursin('[', l)
            line = l  # keep the LAST occurrence (final reported array)
        end
    end
    line === nothing && error("no 'SFmin = [...]' line found in $path")
    inner = line[(findfirst('[', line)+1):(findlast(']', line)-1)]
    return parse.(Float64, split(inner, ','))
end

function main()
    length(ARGS) >= 2 || error("usage: check_sfmin_parity.jl <baseline.out> <new.out> [tol]")
    base = read_sfmin(ARGS[1])
    new  = read_sfmin(ARGS[2])
    tol  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-6
    length(base) == length(new) || error("length mismatch: baseline=$(length(base)) new=$(length(new))")
    @printf("parity: baseline=%s\n        new     =%s\n", basename(ARGS[1]), basename(ARGS[2]))
    @printf("  %-4s %-22s %-22s %-12s\n", "IR", "baseline", "new", "rel_diff")
    maxrel = 0.0; imax = 0
    for i in eachindex(base)
        b = base[i]; n = new[i]
        rel = (b == n) ? 0.0 : abs(b - n) / max(abs(b), abs(n), eps())
        rel > maxrel && (maxrel = rel; imax = i)
        flag = rel <= tol ? "" : "  <-- EXCEEDS"
        @printf("  %-4d %-22.16g %-22.16g %-12.3e%s\n", i, b, n, rel, flag)
    end
    @printf("\n  max rel diff = %.3e at IR index %d  (tol=%.1e)\n", maxrel, imax, tol)
    if maxrel <= tol
        println("  PARITY OK")
    else
        println("  PARITY FAILED")
        exit(1)
    end
end

main()
