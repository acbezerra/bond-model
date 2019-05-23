
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
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end

# Set Parameters ##########################################
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                               :m  => [1.],
                               :gross_delta => [0.02],
                               :kappa  => [25 * 1e-4],
                               :mu_b => [1.0],
                               :xi => [1.0],
                               :iota => [1 * 1e-3])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = Batch.svm_param_values_dict[:lambda]
# svmdict[:kappa] = [25 * 1e-4]
svmdict[:iota] = [.0]
svmdict[:sigmah] = [.25]
# #########################################################


# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)

# Set Targeted Safe Firm
sf_model = "cvm"
sf_comb_num = cvmdf[1, :comb_num]
rf_comb_nums = svmdf[:comb_num]

sf_lambda = unique(svmdf[:lambda])[2]
loc = abs.(svmdf[:lambda] .- sf_lambda) .< 1e-4
if svmdf[loc, :firm_value][1] > cvmdf[1, :firm_value] 
    sf_model = "svm"
    sf_comb_num = svmdf[loc, :comb_num][1]
    rf_comb_nums = [x for x in rf_comb_nums if x != sf_comb_num]
end
# #########################################################

