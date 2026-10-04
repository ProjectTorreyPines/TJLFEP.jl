# Demonstrate + validate AD-exact growth-rate sensitivities ∂γ/∂(plasma inputs)
# at a fixed TGLF-EP operating point. This is the differentiable, grid-independent
# quantity the Fortran code can't provide (the critical factor sfmin sits at a
# discrete keep-flag transition and isn't smoothly differentiable).
#
# Case: ITER IR=83, the marked AE point (kyhat=0.151, width=1.494) at FACTOR_IN=20
# where the Alfvén eigenmode is kept and unstable (γ≈0.017, freq≈-0.33, in-band).
#
# Validation: AD ∂γ/∂knob vs central finite difference on the same configured
# InputTJLF, for every knob.
#
#   module load julia/1.11.7 && export JULIA_DEPOT_PATH=$PSCRATCH/.julia
#   julia -t 8 --project=. dev/ad/demo_gamma_sensitivities.jl

using TJLFEP
using Printf
import ForwardDiff

const DIR = joinpath(@__DIR__, "..", "..", "examples", "ITER")

function build_inputs()
    prof = TJLFEP.readMTGLF(joinpath(DIR, "input.MTGLF"))
    profile = prof[1]; ir_exp = prof[2]
    Options = TJLFEP.readTGLFEP(joinpath(DIR, "input.TGLFEP"), ir_exp)
    ni, Ti, dlnnidr, dlntidr, cs, rmin_ex, gammaE, gammap, omegaGAM =
        TJLFEP.read_expro_for_alpha(joinpath(DIR, "input.EXPRO"), profile, Options.IS_EP; gacode_file = nothing)
    expro = (ni=ni, Ti=Ti, dlnnidr=dlnnidr, dlntidr=dlntidr, cs=cs, rmin_ex=rmin_ex,
             gammaE=gammaE, gammap=gammap, omegaGAM=omegaGAM)
    TJLFEP._apply_runthd_expro_setup!(Options, profile, expro)
    return Options, profile
end

# Float64 γ-vector at a configured InputTJLF (same path the AD pass uses).
gamma_vec(inputF) = ForwardDiff.value.(TJLFEP._tjlf_run_dual(inputF, Float64).eigenvalue[:, 1, 1])

function main()
    @printf("threads = %d\n", Threads.nthreads())
    opts, prof = build_inputs()
    ep = deepcopy(opts); ep.IR = opts.IR_EXP[3]
    ep.KYHAT_IN = 0.151; ep.WIDTH_IN = 1.494; ep.FACTOR_IN = 20.0

    s = gamma_input_sensitivities(ep, prof)
    fu = -abs(prof.omegaGAM[ep.IR])

    # pick the AE mode: in-band (freq<fu) with the largest γ
    cand = findall(<(fu), s.freq)
    m = isempty(cand) ? argmax(s.gamma) : cand[argmax(@view s.gamma[cand])]
    @printf("\nITER IR=%d kyhat=%.3f width=%.3f FACTOR_IN=%.1f  (IS_EP=%d, FREQ_AE_UPPER=%.4f)\n",
            ep.IR, ep.KYHAT_IN, ep.WIDTH_IN, ep.FACTOR_IN, ep.IS_EP, fu)
    @printf("AE mode m=%d: γ=%.6e  freq=%.6e  (in-band=%s)\n",
            m, s.gamma[m], s.freq[m], s.freq[m] < fu ? "yes" : "NO")

    # central-FD reference on the same configured input
    inputF = TJLFEP.TJLF_map(ep, prof)
    TJLFEP._configure_inputTJLF_for_ky!(inputF, ep)

    @printf("\n%-12s  %12s  %12s  %12s  %10s  %12s\n",
            "knob", "base", "dγ/dknob(AD)", "dγ/dknob(FD)", "rel.err", "x·dγ/dx")
    for (i, (fld, idx)) in enumerate(s.knobs)
        x0 = s.base[i]
        δ = max(1e-4, 1e-4 * abs(x0))
        function setval!(inp, val)
            if idx === nothing; setfield!(inp, fld, val); else getfield(inp, fld)[idx] = val; end
        end
        ip = deepcopy(inputF); setval!(ip, x0 + δ); gp = gamma_vec(ip)[m]
        im = deepcopy(inputF); setval!(im, x0 - δ); gm = gamma_vec(im)[m]
        fd = (gp - gm) / (2δ)
        ad = s.dgamma[m, i]
        rel = abs(ad - fd) / (abs(fd) + 1e-12)
        @printf("%-12s  %12.5e  %12.5e  %12.5e  %10.2e  %12.5e\n",
                s.labels[i], x0, ad, fd, rel, s.logsens[m, i])
    end

    # ranked drive/damping by |logarithmic sensitivity|
    order = sortperm([abs(s.logsens[m, i]) for i in 1:length(s.knobs)]; rev = true)
    println("\n--- ranked drive(+)/damping(-) of AE γ by |x·dγ/dx| ---")
    for i in order
        @printf("  %-12s  x·dγ/dx = %+.5e\n", s.labels[i], s.logsens[m, i])
    end
end

main()
