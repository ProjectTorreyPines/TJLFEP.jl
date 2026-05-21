#!/usr/bin/env julia
# Compare preprocessing outputs: Fortran/file reference vs IMAS-generated inputs.
#
# Usage:
#   REF_DIR=build/imas_ref_local TEST_DIR=build/imas_test_local \
#     julia --startup-file=no --project=.. build/compare_preprocess_inputs.jl

using Printf
using TJLFEP

const REF_DIR = abspath(get(ENV, "REF_DIR", "imas_ref_local"))
const TEST_DIR = abspath(get(ENV, "TEST_DIR", "imas_test_local"))

"""Build species index map: ref (Fortran) is -> test (IMAS) is."""
function species_map(ns_ref::Int, ns_test::Int)
    if ns_ref == ns_test
        return Dict(is => is for is in 1:ns_ref)
    elseif ns_ref == 3 && ns_test == 4
        return Dict(1 => 1, 2 => 2, 3 => 4)
    else
        return Dict(is => is for is in 1:min(ns_ref, ns_test))
    end
end

"""Parse `key=value` lines from MTGLF / EXPRO into Dict{String,Float64}."""
function parse_kv_file(path::AbstractString)
    d = Dict{String, Float64}()
    isfile(path) || return d
    for line in readlines(path)
        s = strip(line)
        isempty(s) || !occursin('=', s) && continue
        parts = split(s, '=', limit=2)
        length(parts) < 2 && continue
        key = strip(parts[1])
        val = tryparse(Float64, strip(parts[2]))
        val === nothing && continue
        d[key] = val
    end
    return d
end

function rel_err(a::Float64, b::Float64)
    denom = max(abs(a), 1e-30)
    return abs(b - a) / denom
end

function compare_kv_files(name::AbstractString, ref_path::AbstractString, test_path::AbstractString;
        key_map::Function=identity)
    ref = parse_kv_file(ref_path)
    test = parse_kv_file(test_path)
    ref_keys = Set(keys(ref))
    test_keys = Set(keys(test))

    only_ref = setdiff(ref_keys, test_keys)
    only_test = setdiff(test_keys, ref_keys)

    shared = intersect(ref_keys, test_keys)
    errs = Tuple{String, Float64, Float64, Float64}[]
    for k in shared
        rk = key_map(k)
        rk === nothing && continue
        haskey(test, rk) || continue
        push!(errs, (k, ref[k], test[rk], rel_err(ref[k], test[rk])))
    end
    sort!(errs, by=x -> x[4], rev=true)

    println("\n=== $name ===")
    println("  ref keys: $(length(ref))  test keys: $(length(test))  shared (mapped): $(length(errs))")
    if !isempty(only_ref)
        println("  only in ref ($(length(only_ref))): ", join(sort(collect(only_ref))[1:min(5, end)], ", "), length(only_ref) > 5 ? "..." : "")
    end
    if !isempty(only_test)
        println("  only in test ($(length(only_test))): ", join(sort(collect(only_test))[1:min(5, end)], ", "), length(only_test) > 5 ? "..." : "")
    end
    if isempty(errs)
        println("  no shared keys to compare")
        return NaN
    end
    rels = [e[4] for e in errs]
    @printf("  max rel err = %.4g  mean rel err = %.4g\n", maximum(rels), sum(rels) / length(rels))
    println("  worst 10:")
    @printf("    %-28s %14s %14s %10s\n", "key", "ref", "test", "rel_err")
    for (k, rv, tv, re) in errs[1:min(10, length(errs))]
        @printf("    %-28s %14.6g %14.6g %10.4g\n", k, rv, tv, re)
    end
    return maximum(rels)
end

"""Map MTGLF matrix keys from ref species index to test species index."""
function mtglf_key_map(ref_key::AbstractString, f2j::Dict{Int, Int})
    m = match(r"^(.+)_(\d+)_(\d+)$", ref_key)
    if m !== nothing
        field, ir_str, is_f = m.captures
        is_f = parse(Int, is_f)
        haskey(f2j, is_f) || return nothing
        is_j = f2j[is_f]
        return "$(field)_$(ir_str)_$(is_j)"
    end
    m = match(r"^(.+)_(\d+)$", ref_key)
    if m !== nothing
        field, is_f = m.captures
        is_f = parse(Int, is_f)
        haskey(f2j, is_f) || return nothing
        is_j = f2j[is_f]
        return "$(field)_$(is_j)"
    end
    return ref_key
end

"""Map EXPRO keys (species index in middle)."""
function expro_key_map(ref_key::AbstractString, f2j::Dict{Int, Int})
    m = match(r"^EXPRO_(\w+)_(\d+)_(\d+)$", ref_key)
    if m !== nothing
        field, is_f, ir = m.captures
        is_f = parse(Int, is_f)
        haskey(f2j, is_f) || return nothing
        is_j = f2j[is_f]
        return "EXPRO_$(field)_$(is_j)_$(ir)"
    end
    return ref_key
end

function compare_mtglf_struct(ref_path::AbstractString, test_path::AbstractString; ir_focus::Int=2)
    prof_r, ir_r = readMTGLF(ref_path)
    prof_t, ir_t = readMTGLF(test_path)
    f2j = species_map(prof_r.NS, prof_t.NS)

    println("\n=== MTGLF struct (species-mapped) ===")
    println("  ref NR/NS: ", prof_r.NR, "/", prof_r.NS, "  IR_EXP: ", ir_r)
    println("  test NR/NS: ", prof_t.NR, "/", prof_t.NS, "  IR_EXP: ", ir_t)
    println("  species map (ref→test): ", f2j)

    matrix_fields = (:AS, :TAUS, :RLNS, :RLTS)
    max_rel = 0.0
    for field in matrix_fields
        mat_r = getfield(prof_r, field)
        mat_t = getfield(prof_t, field)
        for is_f in 1:prof_r.NS
            is_j = get(f2j, is_f, nothing)
            is_j === nothing && continue
            is_j > prof_t.NS && continue
            for ir in 1:prof_r.NR
                a = mat_r[ir, is_f]
                b = mat_t[ir, is_j]
                (ismissing(a) || isnan(a) || isnan(b)) && continue
                max_rel = max(max_rel, rel_err(a, b))
            end
        end
    end
    @printf("  matrix fields max rel (all radii): %.4g\n", max_rel)

    println("  at IR_EXP=$ir_focus:")
    @printf("    %-8s %12s %12s %10s\n", "field", "ref(F→J)", "test(J)", "rel_err")
    for field in matrix_fields
        mat_r = getfield(prof_r, field)
        mat_t = getfield(prof_t, field)
        for is_f in 1:prof_r.NS
            is_j = get(f2j, is_f, nothing)
            is_j === nothing && continue
            is_j > prof_t.NS && continue
            a = mat_r[ir_focus, is_f]
            b = mat_t[ir_focus, is_j]
            (ismissing(a) || isnan(a) || isnan(b)) && continue
            re = rel_err(a, b)
            @printf("    %-8s F%d→J%d %12.6g %12.6g %10.4g\n", field, is_f, is_j, a, b, re)
        end
    end

    for (fname, fr, ft) in (("RMIN", prof_r.RMIN, prof_t.RMIN),
                            ("BETAE", prof_r.BETAE, prof_t.BETAE),
                            ("Q", prof_r.Q, prof_t.Q))
        if !ismissing(fr) && !ismissing(ft) && length(fr) >= ir_focus && length(ft) >= ir_focus
            re = rel_err(fr[ir_focus], ft[ir_focus])
            @printf("    %-8s ir=%d  ref=%.6g test=%.6g rel=%.4g\n", fname, ir_focus, fr[ir_focus], ft[ir_focus], re)
        end
    end
    return max_rel
end

function compare_tglfep(ref_path::AbstractString, test_path::AbstractString)
    _, ir_r = readMTGLF(joinpath(REF_DIR, "input.MTGLF"))
    opts_r = readTGLFEP(ref_path, ir_r)
    opts_t = readTGLFEP(test_path, ir_r)

    println("\n=== input.TGLFEP (readTGLFEP) ===")
    skip = Set([:IS_EP])  # Fortran 2 vs IMAS ep_slot-1 (=3)
    max_rel = 0.0
    mismatches = String[]
    for key in fieldnames(typeof(opts_r))
        vr = getfield(opts_r, key)
        vt = getfield(opts_t, key)
        key in skip && continue
        if vr isa Bool
            vr == vt || push!(mismatches, "$key: ref=$vr test=$vt")
        elseif vr isa Number && vt isa Number
            if vr isa AbstractFloat || vt isa AbstractFloat
                re = rel_err(Float64(vr), Float64(vt))
                max_rel = max(max_rel, re)
                re > 1e-10 && push!(mismatches, @sprintf("%s: ref=%g test=%g rel=%.4g", key, vr, vt, re))
            elseif vr != vt
                push!(mismatches, "$key: ref=$vr test=$vt")
            end
        elseif vr isa AbstractVector && vt isa AbstractVector && length(vr) == length(vt)
            for i in eachindex(vr)
                if vr[i] isa Number && vt[i] isa Number
                    re = rel_err(Float64(vr[i]), Float64(vt[i]))
                    max_rel = max(max_rel, re)
                end
            end
        end
    end
    println("  IS_EP (semantic): ref=$(opts_r.IS_EP) test=$(opts_t.IS_EP)  [Fortran 2 ↔ IMAS ep_slot-1]")
    println("  IR_EXP: ref=$(opts_r.IR_EXP) test=$(opts_t.IR_EXP)")
    if isempty(mismatches)
        println("  scalar fields match (within tol)")
    else
        println("  mismatches / large rel err:")
        for m in mismatches[1:min(15, end)]
            println("    ", m)
        end
    end
    return max_rel
end

function main()
    println("Reference: ", REF_DIR)
    println("Test:      ", TEST_DIR)
    for d in (REF_DIR, TEST_DIR)
        isdir(d) || error("missing directory: $d")
    end

    ir_focus = 2
    mtglf_path_r = joinpath(REF_DIR, "input.MTGLF")
    mtglf_path_t = joinpath(TEST_DIR, "input.MTGLF")
    prof_r, ir_r = readMTGLF(mtglf_path_r)
    prof_t, _ = readMTGLF(mtglf_path_t)
    f2j = species_map(prof_r.NS, prof_t.NS)
    if !isempty(ir_r)
        ir_focus = ir_r[1]
    end

    compare_kv_files("input.MTGLF (raw keys, species-mapped)", mtglf_path_r, mtglf_path_t;
        key_map=k -> mtglf_key_map(k, f2j))
    compare_mtglf_struct(mtglf_path_r, mtglf_path_t; ir_focus=ir_focus)

    compare_kv_files("input.EXPRO (species-mapped)", joinpath(REF_DIR, "input.EXPRO"),
        joinpath(TEST_DIR, "input.EXPRO"); key_map=k -> expro_key_map(k, f2j))

    compare_tglfep(joinpath(REF_DIR, "input.TGLFEP"), joinpath(TEST_DIR, "input.TGLFEP"))
end

main()
