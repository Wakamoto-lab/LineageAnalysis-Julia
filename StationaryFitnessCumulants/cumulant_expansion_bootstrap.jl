### Codes for cumulant expansions of stationary data ###
#=
Julia Version 1.6.2
Platform Info:
  OS: macOS (x86_64-apple-darwin19.6.0)
  CPU: Intel(R) Core(TM) i9-10910 CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
=#

using StatsBase, Statistics, OffsetArrays, TaylorSeries, DelimitedFiles

struct FundFunc
    ff
    ff′
end

## Define functions ##
# dc for division counts in retrospective view
ff(ξ,dc;rs=true) = rs ? log(mean(2.0.^((ξ-1)*dc))) - log(mean(2.0.^-dc)) : log(mean(2.0.^(ξ*dc)))
ff′(ξ,dc;rs=true) = rs ? mean(2.0.^((ξ-1)*dc).*dc*log(2))/mean(2.0.^((ξ-1)*dc)) : mean(2.0.^(ξ*dc).*dc*log(2))/mean(2.0.^(ξ*dc))
ff′_norm(ξ,ff′) = (ff′(ξ) - ff′(0))/(ff′(1) - ff′(0))
ssrel₁(ff,ff′) = 1 - ff′(0)/ff(1)
ssrel₂(ff,ff′) = ff′(1)/ff(1) - 1


## Import data ##
# ATTENTION: Rewrite the file paths properly before running this script
# dc_early for division counts from early stationary
dc_early = readdlm("early/ndiv_ealry.csv", ',', Int, '\n')[:]
# dc_late for division counts from late stationary
dc_late = readdlm("late/ndiv_late.csv", ',', Int, '\n')[:]


## 20000 times random resampling with retrospective weights ##
num_bs = 20000 #Number of random resamplimng
ff_bs_early = Vector{FundFunc}(undef,num_bs)
ff_bs_late = similar(ff_bs_early)
size_resample_early = length(dc_early)
size_resample_late = length(dc_late)
for i=1:num_bs
    dc_resample_early_rs = sample(dc_early,size_resample_early)
    dc_resample_late_rs = sample(dc_late,size_resample_late)

    ff_bs_early[i] = FundFunc(ξ->ff(ξ,dc_resample_early_rs),ξ->ff′(ξ,dc_resample_early_rs))
    ff_bs_late[i] = FundFunc(ξ->ff(ξ,dc_resample_late_rs),ξ->ff′(ξ,dc_resample_late_rs))
end


## Bootstrap estimates for cumulant expansions ##
n = 20
t = Taylor1(Float64,n)

# Value at i=-1 is log growth ratio ff(1) = τΛ
# Value at i∈0:n is i-th order cumulant/i! (i-th coefficient of Taylor series)
taylorcoef_ff_early = OffsetArray(zeros(n+2,num_bs),-1:n,1:num_bs)
taylorcoef_ff_late = OffsetArray(zeros(n+2,num_bs),-1:n,1:num_bs)
for i=1:num_bs
    taylorcoef_ff_early[-1,i] = ff_bs_early[i].ff(1)
    taylorcoef_ff_late[-1,i] = ff_bs_late[i].ff(1)
    taylorcoef_ff_early[0:n,i] .= ff_bs_early[i].ff(t)[:]
    taylorcoef_ff_late[0:n,i] .= ff_bs_late[i].ff(t)[:]
end

wn_bs_early = zeros(n,num_bs)#OffsetArray(zeros(n,num_bs),0:n,1:num_bs)
wn_bs_late = similar(wn_bs_early)
for i=1:num_bs
    wn_bs_early[:,i] .= cumsum(taylorcoef_ff_early[1:n,i],dims=1)/taylorcoef_ff_early[-1,i]
    wn_bs_late[:,i] .= cumsum(taylorcoef_ff_late[1:n,i],dims=1)/taylorcoef_ff_late[-1,i]
end


## Save results ##
open("cumulant_stationary_bs_retrospective_20000_NEW.txt", "w") do io
    writedlm(io, ["order" "Early(mean)" "Early(2SD)" "Late(mean)" "Late(2SD)"])
    writedlm(io,
        [1:n mean(wn_bs_early,dims=2) 2*std(wn_bs_early,dims=2) mean(wn_bs_late,dims=2) 2*std(wn_bs_late,dims=2)])
end
