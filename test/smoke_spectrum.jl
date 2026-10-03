# Manual smoke test for process_in=3 (spectrum mode). Not part of runtests.jl.
using TJLFEP

const CASE = joinpath(@__DIR__, "..", "examples", "DIIID_202017C42_500ms_v3.1")
const GACODE = joinpath(CASE, "input.gacode")
const TGLFEP = joinpath(CASE, "input_spectrum.TGLFEP")

out_dir = mktempdir()
@info "running process_in=3 spectrum" GACODE TGLFEP out_dir

# diagnostics on the profile grid
let
    Opts, prof, expro = TJLFEP.preprocess_gacode_inputs(GACODE, TGLFEP)
    TJLFEP._apply_runthd_expro_setup!(Opts, prof, expro)
    @info "profile diag" NR=prof.NR IRS=prof.IRS IR_EXP=Opts.IR_EXP SCAN_N=Opts.SCAN_N
    ir = Opts.IR_EXP[1]
    @info "Q diag" ir Q_ir=prof.Q[ir] nan_in_Q=count(isnan, prof.Q) len_Q=length(prof.Q) RMIN_end=prof.RMIN[end] RMAJ_ir=prof.RMAJ[ir] KAPPA_ir=prof.KAPPA[ir]
end

r = run_gacode_scan_task(GACODE, TGLFEP, 1; out_dir=out_dir, use_gpu=false, printout=true)

@info "result" ir=r.ir width=r.width kymark=r.kymark
@assert r.spectra !== nothing "spectra missing from result"
for mode in (1, 2, 4)
    s = r.spectra[mode]
    @info "mode $mode" nky=length(s.ky) size_gamma=size(s.gamma) any_finite=any(isfinite, s.gamma)
    @assert length(s.ky) == 30 "expected nky=30, got $(length(s.ky))"
    @assert size(s.gamma, 1) == 30 "gamma rows != 30"
    @assert all(isfinite, s.ky) "non-finite ky"
    @assert any(isfinite, s.gamma) "no finite gamma in mode $mode"
end

println("\nfiles written to out_dir:")
for f in readdir(out_dir)
    println("  ", f)
end

# Save Julia spectra for offline comparison with the Fortran golden.
open(joinpath(@__DIR__, "julia_spectrum_dump.txt"), "w") do io
    for mode in (1, 2, 4)
        s = r.spectra[mode]
        println(io, "# mode_in=", mode)
        for i in 1:length(s.ky)
            print(io, lpad(string(round(s.ky[i]; digits=4)), 9))
            for n in 1:size(s.gamma, 2)
                print(io, "  ", round(s.gamma[i,n]; digits=7), "  ", round(s.freq[i,n]; digits=7))
            end
            println(io)
        end
    end
end

# Print mode_in=1 spectrum head (background+EP, thermal gradients retained)
s1 = r.spectra[1]
println("\nmode_in=1 spectrum (ky, gamma_1, freq_1):")
for i in 1:5
    println("  ", round(s1.ky[i]; digits=4), "  ", round(s1.gamma[i,1]; digits=6), "  ", round(s1.freq[i,1]; digits=6))
end
println("SMOKE OK")
