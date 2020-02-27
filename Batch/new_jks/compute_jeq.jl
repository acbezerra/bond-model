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
         "ModelPlots", "JointEqStructs", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# * Create Joint Firm Objects
k_otc = 50 * 1e-4
k_ep = 25 * 1e-4

# Measure of Safe Firms
mu_s = .2

rt_iota = 2 * 1e-4
rt_lambda = .2
rt_sigmah = .225 #Batch.svm_param_values_dict[:sigmah][2]

jfotc, jfep = JointEq.otc_ep_jfs(k_otc, k_ep;
                                 mu_s=mu_s,
                                 rt_iota=rt_iota,
                                 rt_lambda=rt_lambda,
                                 rt_sigmah=rt_sigmah)


# * Make Directories and File Names
# Make results directories
jks_fpath = JointEq.make_jeq_jks_fpath(jfep)


# * Full Information - Bond Contract and Equilibrium Results
fidf_fpath_name = string(jks_fpath, "/",
                         JointEqStructs.eq_type_dict["full_info"][:dfn], ".csv")


fidf, fieqdf = JointEq.get_fi_results(jfep, fidf_fpath_name;
                                      load_df=true,
                                      recompute_df=false,
                                      save_df=false)

# Get unique results
uq_fidf = JointEqStructs.get_unique_df(jks_fpath, :fi; del_old_files=true)


# * Misrepresentaion
# Main results file
misrepdf_fpath_name = string(jks_fpath, "/", JointEqStructs.eq_type_dict["misrep"][:dfn], ".csv")
misrepdf = JointEq.get_misrep_results(jfep, fidf, misrepdf_fpath_name;
                                      recompute_df=false, save_df=false)


# Get unique Results
uq_misrepdf = JointEqStructs.get_unique_df(jks_fpath, :misrep)#; del_old_files=true)


# * Joint Equilibrium
# ** Identify the Joint Equilibrium Case
jeqid = JointEq.get_joint_eq_inputs(jfep, jks, jks_fpath)


# ** Compute Pooling Equilibrium Payoffs
pooldf = JointEq.get_jeq_results(jfep, jks, 
                                :pool, :exp_firm_value, 
                                jeqid; load_df=true,
                                recompute_df=false,
                                 save_df=true)


# ** Compute Separating Equilibrium Payoffs
sepdf = JointEq.get_jeq_results(jfep, jks, 
                                :sep, :st_firm_value, 
                                jeqid; load_df=false,
                                recompute_df=true,
                                save_df=true)

