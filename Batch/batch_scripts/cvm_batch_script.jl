using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using Dierckx
using JLD
using CSV

start_tic = time_ns()

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
# include(string(module_path, "/", "TestFunctions.jl"))
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


println(string("ARGUMENTS: ", ARGS))
# ################ SYS ARGUMENTS ################
# Capture Combination Number:
m_comb_num = parse(Int, ARGS[1])

# ###############################################

firm_obj_fun=:firm_value

# Start Timer
tic = time_ns()

for m in Batch.cvm_param_values_dict[:m]
    # Create Batch & Model Objects
    bt, cvm = Batch.get_bt_cvm(; m=m, m_comb_num=m_comb_num)

    # Create Directories
    bt = Batch.mk_comb_res_dirs(bt)

    # Solve Debt at Par for all Coupon Values
    soldf = @time Batch.get_cvm_debt_at_par(bt, cvm;
                                            mu_b=cvm.mu_b, m=cvm.m,
                                            save_soldf=true,
                                            soldf_name=Batch.soldf_name)

    # Solve for Optimal Capital Structure
    optdf = Batch.optimal_cvm_capital_struct(bt, cvm;
                                             firm_obj_fun=firm_obj_fun,
                                             df=soldf,
                                             save_optdf=true,
                                             optdf_name=Batch.optdf_name)

end

println(string("Total Script Run Time: ", (time_ns() - start_tic)/1e9/60., " minute(s)."))
