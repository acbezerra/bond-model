
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
# Diagnostics Script Inputs
# 1. (OPTIONAL) choose model: cvm = 0, svm = 1 (default)
# 2. (OPTIONAL) recompute: false = 0 (default), true = 1
# 3. (OPTIONAL) diganose: false = 0, true = 1 (default)
# 4. (OPTIONAL) maturity: m = 1. (default)
# 5. (OPTIONAL) compilation: false = 0, true = 1 (default)
# 6. (OPTIONAL) return results: false = 0, true = 1 (default)

# ################ SYS ARGUMENTS ################
println(string("ARGUMENTS: ", ARGS))
# Choose model
model = "svm"
if size(ARGS, 1) > 0
    model = parse(Int, ARGS[1]) == 1 ? "svm" : "cvm"
end

# Wether to recompute results
recompute = false
if size(ARGS, 1) > 1
    recompute = parse(Bool, ARGS[2])
end

# Run Diagnostics
diagnose = true
if size(ARGS, 1) > 2
    diagnose = parse(Bool, ARGS[3])
end

# Debt Maturity
m = 1.
if size(ARGS, 1) > 3
    m = parse(Float64, ARGS[4])
end

# Compile Results to find Optimal Capital Structure
compilation = true
if size(ARGS, 1) > 4
    compilation = parse(Bool, ARGS[5])
end

# Whether to return the results to user interface
return_results = true
if size(ARGS, 1) > 5
    return_results = parse(Bool, ARGS[6])
end

# Firms' Objective Functions
obj_funs = Dict{Symbol, Symbol}(:fv => :firm_value)


# Dictionary to store results
res = Dict{Symbol, Any}(:model => model,
                        :m => m)
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
