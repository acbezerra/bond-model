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
module_path = string(main_path, "/", "bond-model/modules/")
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "ModelPlots", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# * Load Full Information Results
# Set Parameters ##########################################
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                               :m  => [1.],
                               :gross_delta => [0.02],
                               :kappa  => [25 * 1e-4],
                               :mu_b => [1.0],
                               :xi => [1.0],
                               :iota => [.0])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = Batch.svm_param_values_dict[:lambda]
# svmdict[:kappa] = [25 * 1e-4]
svmdict[:sigmah] = [.225]
# #########################################################


# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)

# * Adjustments
# ** Form Sage and Risky Firm's DataFrames
# Set Targeted Safe Firm
sf_model = "cvm"
sf_comb_num = cvmdf[1, :comb_num]
rf_comb_nums = svmdf[:, :comb_num]

# sf_lambda = unique(svmdf[:, :lambda])[2]
# loc = abs.(svmdf[:, :lambda] .- sf_lambda) .< 1e-4
# if svmdf[loc, :firm_value][1] > cvmdf[1, :firm_value]
#     sf_model = "svm"
#     sf_comb_num = svmdf[loc, :comb_num][1]
#     rf_comb_nums = [x for x in rf_comb_nums if x != sf_comb_num]
# end

tmp = deepcopy(DataFrame(cvmdf[1, :]))
tmp[!, :eq_deriv_min_val] .= NaN
tmp[!, :eq_negative] .= false
svmdf = vcat(tmp, svmdf)
svmdf[1, :sigmah] = svmdf[2, :sigmah]
svmdf[1, :lambda] = 0.0
cvmdf = deepcopy(DataFrame(cvmdf[1, :]))
# #########################################################

# ** Set Iota and Kappa variables to basis points
# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    cvmdf[!, x] = cvmdf[:, x] .* 1e4
    svmdf[!, x] = svmdf[:, x] .* 1e4
end
# #########################################################

# * Find Cut-off values
# Cut-off Values ##########################################
xgrid = range(minimum(svmdf[:, :lambda]), stop=maximum(svmdf[:, :lambda]), length=10^5)

# lambda : RM Firm Value = NRM Firm Value
fv_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                        :firm_value, cvmdf[1, :firm_value];
                                        xgrid=xgrid)

# lambda : RM Firm Value = NRM Firm Value
mbr_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                         :MBR, cvmdf[1, :MBR];
                                         xgrid=xgrid)
# #########################################################
# * Restrict x-axis values
max_xval = .5 # maximum(svmdf[:, :lambda])
svmdf_input = deepcopy(svmdf[abs.(svmdf[:, :lambda]) .<= max_xval, :])

# * Plots
# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value, :MBR], 
                                       deepcopy(cvmdf), svmdf_input;
                                       nrmp=true,
                                       fv_xvar=fv_lambda,
                                       mbr_xvar=mbr_lambda,
                                       color_rm_region=false,
                                       color_nrm_region=false,
                                       color_conflict_region=false,
                                       color_misrep_region=false, 
                                       save_fig=true)
