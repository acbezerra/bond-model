
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
using PyPlot
# using PyCall
using Seaborn
using LaTeXStrings

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
modls = ["Batch", "ModelObj", 
         "AnalyticFunctions", "BondPrInterp",
         "EqFinDiff", "FullInfoEq",
         "ModelPlots", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


# Set Safe Type
s_iota = 0.0 # 2.5 * 1e-4

# ###############################################################################
# Get OTC Results ###############################################################
# ###############################################################################
# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [x for x in Batch.cvm_param_values_dict[:iota] 
                                                  if .&(x >= s_iota, x <= 20. * 1e-4)])
svmdict = deepcopy(cvmdict)
svmdict[:kappa] = [25 * 1e-4]
svmdict[:lambda] = [.2]
svmdict[:iota] = [.0]
svmdict[:sigmah] = Batch.svm_param_values_dict[:sigmah]
# ##############################################################################


# Get Safe and Risky Firms' Full Info Optimal Results ##########################
firm_obj_fun = :firm_value 
cvmdf, _, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)
# ##############################################################################

# Safe Type's Firm Value as a function of Transaction Costs Kappa
cond = abs.(cvmdf[:, :iota] .- s_iota) .< 1e-5
fi_fv_fun = Dierckx.Spline1D(cvmdf[cond, :kappa] .* 1e4, cvmdf[cond, :firm_value]; 
                             k=3, bc="extrapolate")
# ###############################################################################


# ###############################################################################
# Load Data #####################################################################
# ###############################################################################
# Set Path
jks_fpath = string(main_path, "/Julia/Results/JEQ/kappa_25.00_bp__sigmal_0.150/m_1.00_pcr_12.00")

# Load Results
fidf = JointEq.process_results(jks_fpath, "full_info")
misrepdf = JointEq.process_results(jks_fpath, "misrep")
pooldf = JointEq.process_results(jks_fpath, "pooling")
sepdf = JointEq.process_results(jks_fpath, "separating")

# Keep only the results for which mu_s = .2
mu_s = .2
pooldf = pooldf[pooldf[:, :mu_s] .== mu_s, :]
sepdf = sepdf[sepdf[:, :mu_s] .== mu_s, :]

# Iota and kappa in basis points
fidf[!, :iota] .= fidf[:, :iota] .* 1e4
fidf[!, :kappa] .= fidf[:, :kappa] .* 1e4
for df in [misrepdf, sepdf, pooldf]
    for ft in [:s_, :r_]
        df[!, Symbol(ft, :iota)] .= df[:, Symbol(ft, :iota)] .* 1e4
    end
    
    df[!, :kappa] .= df[:, :kappa] .* 1e4
end
# ###############################################################################


# ###############################################################################
# Pick Safe Type and Interpolate Firm Value, MBR and Leverage on iota and sigmah 
# ###############################################################################
cond = [abs.(iota .- s_iota * 1e4) .> 1e-5 for iota in fidf[:, :iota]]
fi_fd = ModelPlots.interp_z_values(fidf[cond, :])
mp_fd = ModelPlots.interp_z_values(misrepdf)
pool_fd = ModelPlots.interp_z_values(pooldf)
sep_fd = ModelPlots.interp_z_values(sepdf)
# ###############################################################################


# ###############################################################################
# Generate Payoff Functions ##################################################### 
# ###############################################################################
fi_funs = ModelPlots.fi_payoff_functions(fi_fd)
mp_funs = ModelPlots.misrep_payoff_functions(deepcopy(fi_funs), mp_fd)
pool_funs = ModelPlots.jeq_payoff_functions(deepcopy(fi_funs), pool_fd; eq_type="pooling")
sep_funs = ModelPlots.jeq_payoff_functions(deepcopy(fi_funs), sep_fd; eq_type="separating")
# ###############################################################################


# #########################################################################
# Equilibrium Type Indicators and Equilibrium-Specific Payoffs ############
# #########################################################################
fi_fv=fidf[1, :firm_value]
println(string("s_fi_fv 1: ", fi_fv))
println(string("s_fi_fv 2: ", fidf[abs.(fidf[:iota] .- s_iota * 1e4) .< 1e-5, :firm_value][1]))
k_ep = 25.
k_otc=32.5


# Safe Type's Firm Value in the EP Full Information Equilibrium:
s_fi_fv = fi_fv_fun(k_ep) # fidf[abs.(fidf[:iota] .- s_iota * 1e4) .< 1e-5, :firm_value][1]
println(string("s_fi_fv 3: ", s_fi_fv))

# Safe Type's Firm Value in the EP Full Information Equilibrium:
s_otc_fv = fi_fv_fun(k_otc)
println(string("s_otc_fv: ", s_otc_fv))


fi_fv=fidf[1, :firm_value]
k_otc=32.5
fun_dict = ModelPlots.get_contour_equilibria_funs(fi_funs,
                                                  mp_funs,
                                                  pool_funs,
                                                  sep_funs,
                                                  fi_fv,
                                                  fi_fv_fun,
                                                  k_otc)

yvals = [x for x in pooldf[:, :r_sigmah] if .&(!isnan(x), x>.0)]
xvals = [x for x in pooldf[:, :r_iota] if .&(!isnan(x), x>.0)]
pltd = ModelPlots.get_eq_contour_mesh_grid(xvals, yvals, fun_dict, N=10^3)
