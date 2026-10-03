using Pkg; Pkg.activate(normpath(@__DIR__, "..", ".."))
using Plots, Printf
gr()

# Parse the 20-radii e2e log: lines like
#   "   2 | sf=   0.9374 ky=... w=...  132s | sf=   0.8203 ky=... w=...   40s | 5.09e-06 <-- MISMATCH"
# -> (IR, golden_sfmin, batched_si_sfmin).
LOG = get(ENV, "LOG", normpath(@__DIR__, "si_full20_55806342.out"))
rowre = r"^\s*(\d+)\s*\|\s*sf=\s*([\d.eE+-]+).*?\|\s*sf=\s*([\d.eE+-]+).*?\|\s*([\d.eE+-]+)"
ir = Int[]; gold = Float64[]; bsi = Float64[]; dsf = Float64[]
for ln in eachline(LOG)
    m = match(rowre, ln); m === nothing && continue
    push!(ir, parse(Int, m[1])); push!(gold, parse(Float64, m[2]))
    push!(bsi, parse(Float64, m[3])); push!(dsf, parse(Float64, m[4]))
end
isempty(ir) && error("no data rows parsed from $LOG")
nmis = count(>(0.02), dsf)
@printf("parsed %d radii from %s  (%d mismatch >2%%)\n", length(ir), basename(LOG), nmis)

mis = dsf .> 0.02
p = plot(ir, gold; seriestype=:line, marker=:circle, ms=5, lw=2, color=:black,
         label="dense grid (golden)", xlabel="radial index IR_EXP", ylabel="marginal sfmin",
         yscale=:log10, legend=:topleft, title="sfmin vs radius  (nb16, full grid 8×8×4×4)",
         size=(950, 560), left_margin=6Plots.mm, bottom_margin=5Plots.mm)
plot!(p, ir, bsi; seriestype=:line, marker=:diamond, ms=5, lw=2, color=:royalblue,
      label="batched_si (fixed-shift hybrid)")
# highlight mismatched radii
any(mis) && scatter!(p, ir[mis], bsi[mis]; marker=:xcross, ms=9, color=:red, lw=3,
                     label="mismatch (>2%)")
for i in eachindex(ir)
    mis[i] && annotate!(p, ir[i], bsi[i]*1.25, text(@sprintf("%.0f%%", 100*dsf[i]), 7, :red, :left))
end
out = normpath(@__DIR__, "sfmin_vs_radius.png")
savefig(p, out)
println("wrote ", out)
