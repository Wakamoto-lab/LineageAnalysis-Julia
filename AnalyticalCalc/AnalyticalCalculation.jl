### Scripts for Figure 3 ###
#=
Julia Version 1.6.2
Platform Info:
  OS: macOS (x86_64-apple-darwin19.6.0)
  CPU: Intel(R) Core(TM) i9-10910 CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
=#

using TaylorSeries, PyPlot, Printf, Roots, Distributions

## Define functions ##
# Binomial/Negative-Binomila scheme
ff_bin(ξ,ϕ) = ϕ==1 ? 2^ξ-1 : log((1-ϕ)*2^ξ+ϕ)/(1-ϕ)
ffprime_bin(ξ,ϕ) = log(2)/(1-(1-2.0^(-ξ))*ϕ)
cumcontrib_bin(ϕ,n) = cumsum(ff_bin(Taylor1(Float64,n),ϕ)[1:end])/ff_bin(1,ϕ)

ssrel₁_bin(ϕ) = 1 - ffprime_bin(0,ϕ)/ff_bin(1,ϕ)
ssrel₂_bin(ϕ) = ffprime_bin(1,ϕ)/ff_bin(1,ϕ) - 1

ssequal(ϕ) = (2-ϕ)*log(2-ϕ)/(4-ϕ)/(1-ϕ) - 0.5log(2)


## Settings for plots ##
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams.o

rcParams["font.family"] = ["sans-serif"]#["serif"]#
rcParams["mathtext.fontset"] = "dejavuserif"#dejavusans

ls_arr = ["-","--",":","-."]
markers=["o","x","s","^"]


## Generate Figure 3 ##
ϕ₀ = find_zero(ssequal,(0.1,0.9))
#@show ϕ₀

# Figure 3A
subplot(221); text(-6,0.26,"A", fontsize=15)
Dmax = 20
D̄ = Dmax*(1-ϕ₀)
#@show D̄
xrange=0:Dmax
ϕ_arr = [0.3,ϕ₀,1,1.6]
for i in eachindex(ϕ_arr)
    ϕ = ϕ_arr[i]
    if ϕ < 1
        scatter(xrange,pdf.(Binomial(round(Int64,D̄/(1-ϕ)),1-ϕ),xrange),label=@sprintf("\$\\phi=%.2f\$",ϕ),
            marker=markers[i],s=15,lw=1)
    elseif ϕ==1
        scatter(xrange,pdf.(Poisson(D̄),xrange),label=@sprintf("\$\\phi=%.2f\$",ϕ),
            marker=markers[i],s=15,lw=1)
    else
        scatter(xrange,pdf.(NegativeBinomial(D̄/(ϕ-1),1/ϕ),xrange),label=@sprintf("\$\\phi=%.2f\$",ϕ),
            marker=markers[i],s=15,lw=1)
    end
end
xlabel(L"D",fontsize=10)
ylabel(L"Q_{\mathrm{cl}}\left(D\right)",fontsize=10)
xticks(0:5:20,fontsize=10)
yticks(fontsize=10)
#legend(fontsize=8,framealpha=1)
legend(("\$\\phi=0.3\$","\$\\phi=\\phi_0\$","\$\\phi=1\$","\$\\phi=1.6\$"),fontsize=8,framealpha=1)

# Figure 3B
subplot(222); text(-1.9,1.07,"B", fontsize=15)
n = 10
for i in eachindex(ϕ_arr)
    ϕ = ϕ_arr[i]
    scatter(1:n, cumcontrib_bin(ϕ,n),label=@sprintf("\$\\phi\$=%.2f",ϕ),
        marker=markers[i],s=15,lw=1)
end
plot(1:n,ones(n),ls=":",color="grey")
xticks(1:10,fontsize=10)
yticks(fontsize=10)
xlabel(L"n",fontsize=10)
ylabel(L"W_n^{\left(D\right)}",fontsize=10)
#legend(fontsize=10)

# Figure 3C
subplot(223); text(-0.17,1.36,"C", fontsize=15)
ϕrange_bin = 0:1e-2:0.99
ϕrange_pois = [1.0]
ϕrange_negbin = 1.01:1e-2:1.6
ϕrange_ϕ₀ = [ϕ₀]
l1, = plot(ssrel₁_bin.(ϕrange_bin),ssrel₂_bin.(ϕrange_bin),lw=1,zorder=1)
label1 = "Binomial (\$\\phi<1\$)"
l2 = scatter(ssrel₁_bin.(ϕrange_pois),ssrel₂_bin.(ϕrange_pois),color="black",lw=1,s=15,zorder=3)
label2 = "Poisson (\$\\phi=1\$)"
l3, = plot(ssrel₁_bin.(ϕrange_negbin),ssrel₂_bin.(ϕrange_negbin),lw=1,zorder=2)
label3 = "Negative binomial (\$\\phi>1\$)"
l4 = scatter(ssrel₁_bin.(ϕrange_ϕ₀),ssrel₂_bin.(ϕrange_ϕ₀),
    color="black",marker="o",lw=1,s=15,facecolor="none",zorder=4)
label4 = L"S^{\left(1\right)}_\mathrm{KL}\left[D\right]=S^{\left(2\right)}_\mathrm{KL}\left[D\right]\ (\phi=\phi_0)"

plot(0:1e-2:0.55,0:1e-2:0.55,ls=":",color="grey",lw=1)
xticks(fontsize=10)
yticks(fontsize=10)
xlabel(L"S^{\left(1\right)}_\mathrm{KL}\left[D\right]/\tau\Lambda",fontsize=10)
ylabel(L"S^{\left(2\right)}_\mathrm{KL}\left[D\right]/\tau\Lambda",fontsize=10)
#legend(fontsize=7,framealpha=1)
legend((l1,l2,l3,l4),(label1,label2,label3,label4),fontsize=6,framealpha=1)

# Figure 3D
subplot(224); text(-0.32,1.08,"D", fontsize=15)
ξrange = 0:1e-2:1
for i in eachindex(ϕ_arr)
    ϕ = ϕ_arr[i]
    plot(ξrange, (ffprime_bin.(ξrange,ϕ).-ffprime_bin(0,ϕ))/(ffprime_bin(1,ϕ)-ffprime_bin(0,ϕ)),
        label=@sprintf("\$\\phi\$=%.2f",ϕ),lw=1)
end
plot(ξrange,ξrange,ls=":",color="grey",lw=1)
xticks(fontsize=10)
yticks(0:0.5:1,fontsize=10)
xlabel(L"\xi",fontsize=10)
#ylabel(L"(K^\prime(\xi)-K^\prime(0))/(K^\prime(1)-K^\prime(0))",fontsize=10)
ylabel("Rescaled \$K^\\prime\\left(\\xi\\right)\$",fontsize=10)
#legend(fontsize=10)

subplots_adjust(wspace=0.4,hspace=0.4)
savefig("AnalyticalCalculation_demo.pdf")
;
