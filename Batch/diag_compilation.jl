
using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
modls = ["ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "Batch"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# #########################################################
# ######################## INPUTS #########################
# #########################################################
# Choose model
model="cvm"

# Debt Maturity
m = 1.

# Firms' Objective Functions
obj_funs = Dict{Symbol, Symbol}(:fv => :firm_value)

# Run Diagnostics
diagnose=true

# Compile Results to find Optimal Capital Structure
compilation=true

# Dictionary to store results
res = Dict{Symbol, Any}(:model => model,
                        :m => m)

# Wether to recompute results
recompute=true

# Whether to return the results to user interface
return_results=true
# #########################################################


bt = Batch.get_bt(; model=model, comb_num=1, display_msgs=false)
# #########################################################
# ###################### DIAGNOSTICS ######################
# #########################################################
if diagnose
    diagdf = Batch.diagnosis(bt)
    res[:diag] = diagdf
    # push!(DFs, diagdf)
end
# #########################################################


# #########################################################
# ###################### COMPILATION ######################
# #########################################################
if compilation
    for of in keys(obj_funs)
        if model == "svm"
            optDF = @time Batch.compile_svm_opt_results(bt; m=m,
                                                        firm_obj_fun=obj_funs[of],
                                                        recompute_comb_opt_res=recompute)
        elseif model == "cvm"
            if recompute
                optDF = @time Batch.compile_cvm_opt_results(; m=m,
                                                            firm_obj_fun=obj_funs[of])
            else
                optDF = Batch.load_cvm_opt_results_df(; m=m)
            end
        end
        res[Symbol(:opt_, of)] = optDF
        # push!(DFs, optDF)
    end
end
# #########################################################

if return_results
    return res
end
