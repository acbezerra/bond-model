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
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "ModelPlots", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


# * Load Full Information Results
# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [.15, .2, .25, .30],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :kappa  => [25 * 1e-4],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [x for x in Batch.cvm_param_values_dict[:iota] 
                                                  if x <= 20. * 1e-4])

# anything goes here. For the CVM, we don't
# use the SVM df.
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.2]

# #########################################################
# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 

# #########################################################
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)

# * Adjustments
# ** Form Safe and Risky Firms' DataFrames               
sf_cond = .&(abs.(cvmdf[:, :iota] .- 10 * 1e-4) .< 1e-5,
             abs.(cvmdf[:, :sigmal] .- .15) .< 1e-5)      
sfdf = deepcopy(cvmdf[sf_cond, :])

rf_cond = abs.(cvmdf[:, :iota] .- .0) .< 1e-5 #),
             #abs.(cvmdf[:, :sigmal] .- .25) .< 1e-5)
rfdf = deepcopy(cvmdf[rf_cond, :])

# ** Make CVM results work with EP constructor
rfdf[!, :sigmah] = rfdf[!, :sigmal]

# ** Set Iota and Kappa variables to basis points
# Iota and Kappa in Basis Points
for x in [:iota, :kappa]
    sfdf[!, x] = sfdf[:, x] .* 1e4
    rfdf[!, x] = rfdf[:, x] .* 1e4
end

# * Find Cut-off values               
# Cut-off Values ##########################################
xgrid = range(minimum(rfdf[:, :sigmah]), stop=maximum(rfdf[:, :sigmah]), length=10^5)
                
# iota : RM Firm Value = NRM Firm Value
fv_sigmah = ModelPlots.get_cutoff_value(rfdf, :sigmah,
                                      :firm_value, sfdf[1, :firm_value];
                                      xgrid=xgrid)

# iota : RM Firm Value = NRM Firm Value
mbr_sigmah = ModelPlots.get_cutoff_value(rfdf, :sigmah,
                                       :MBR, sfdf[1, :MBR];
                                       xgrid=xgrid)
# #########################################################

# * Plots
# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:sigmah, [:firm_value, :MBR], 
                                       deepcopy(sfdf), deepcopy(rfdf);
                                       xgrid=xgrid,
                                       fv_xvar=fv_sigmah, mbr_xvar=mbr_sigmah,
                                       color_conflict_region=true,
                                       save_fig=true)
