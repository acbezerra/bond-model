# %% codecell
using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
# using PyPlot
using Seaborn

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
# include(string(module_path, "/", "TestFunctions.jl"))
modls = ["Batch", "ModelObj", "AnalyticFunctions",
         "BondPrInterp", "EqFinDiff",
         "JointEqStructs", "FullInfoEq", "JointEq", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

ENV["LINES"] = 750
ENV["COLUMNS"] = 10000

# %%
fpath = string(main_path, "/Julia/Batch/new_jks/")
include(string(fpath, "/jeq_calculator.jl"))



# %%
tcdf

# %%
misrepdf

# %%
df

# %% codecell
# * SYS & Other Args
# println(string("ARGUMENTS: ", ARGS))
# ################ SYS ARGUMENTS ################
# Position in the Safe Firm Measure array:
# pair_num = parse(Int, ARGS[1])
# rerun_fi = parse(Bool, ARGS[2])
# rerun_misrep = parse(Bool, ARGS[3])
# rerun_pool = parse(Bool, ARGS[4])
# rerun_sep = parse(Bool, ARGS[5])

pair_num=15
rerun_fi=false
rerun_misrep=false
rerun_pool=false
rerun_sep=false
# ##################################################

load_df=true
save_df=true

pool_obj_fun=:st_firm_value
sep_obj_fun=:st_firm_value


println(string("load_df: ", load_df))
println(string("rerun_fi: ", rerun_fi))
println(string("rerun_misrep: ", rerun_misrep))
println(string("rerun_pool: ", rerun_pool))
println(string("rerun_sep: ", rerun_sep))
println(string("save_df: ", save_df))
println(string("pool_obj_fun: ", pool_obj_fun))

st_rt_comb = Dict{Symbol, Float64}()

# Transaction Costs and Measure of safe firms
st_rt_comb[:k_ep] = 25 * 1e-4
st_rt_comb[:k_otc] = 50 * 1e-4

# ATTENTION!!!!!! ############################################
# Value must be in the svm_param_values_dict, otherwise there are
# no precompiled results for the SVM objects!
st_rt_comb[:lambda] = Batch.svm_param_values_dict[:lambda][4]

# * Form Joint Firms Objects
# ** Extract list with all possible combinations
# Set Objective Functions
firm_obj_fun = :firm_value

# Transaction Costs and Volatility Risk Parameters
pardict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :mu_b => [1.0],
                                        :xi => [1.0])

# Risk Type Parameters
st_iota_vec = [.0]
rt_iota_vec = [.0]
iota_vec = unique(vcat(st_iota_vec, rt_iota_vec))
lambda_vec = [st_rt_comb[:lambda]]
sigmah_vec = Batch.svm_param_values_dict[:sigmah]
mu_s_vec = [.1, .2, .5, .75, .9]

# %%

# ############################################################################
# From HERE ON, CODE SHOULD BE THE SAME! #####################################
# ############################################################################

# %%
# * Functions
function get_types_mu_s_comb_df(types_dict::Dict{Symbol, Array{Float64, N} where N})
    types_order = [:mu_s, :st_iota, :rt_iota,
                   :st_lambda, :rt_lambda, :st_sigmah, :rt_sigmah]

    types_combs = [Array([x1, x2, x3, x4, x5, x6, x7]) for x1 in types_dict[types_order[1]],
                                                           x2 in types_dict[types_order[2]],
                                                           x3 in types_dict[types_order[3]],
                                                           x4 in types_dict[types_order[4]],
                                                           x5 in types_dict[types_order[5]],
                                                           x6 in types_dict[types_order[6]],
                                                           x7 in types_dict[types_order[7]]]

    df = DataFrame(hcat(types_combs...)')
    rename!(df, types_order)
    df[!, :pair_num] .= 1:size(df, 1)

    return df
end

function combination_row_finder(df::DataFrame, pardict::Dict{Symbol,Array{Float64,N} where N})
    cond  = [true for x in 1:size(tcdf, 1)]
    tester(df, x, y) = isnan(y) ? isnan.(df[:, x]) : Array(abs.(df[:, x] .- y ) .< 1e-5)

    for x in keys(pardict)
        bool_sum = sum([tester(df, x, y) for y in pardict[x]]) .> 0
        cond = .&(cond, bool_sum)
    end

    return DataFrame(df[cond, :])
end

# *
# Load Unconstrained Results
cvmdf, svmdf = JointEqStructs.get_joint_types_dfs(pardict, iota_vec,
                                                  lambda_vec, sigmah_vec;
                                                  firm_obj_fun=firm_obj_fun)

# Extract Results
all_typesdict = Dict{Symbol, Array{Float64, N} where N}(:st_iota => st_iota_vec,
                                                        :st_lambda => [NaN],
                                                        :st_sigmah => [NaN],
                                                        :rt_iota => rt_iota_vec,
                                                        :rt_lambda => lambda_vec,
                                                        :rt_sigmah => sigmah_vec,
                                                        :mu_s => mu_s_vec)
tcdf = get_types_mu_s_comb_df(all_typesdict)

# ** Optional
# # Set Parameters for the case:
# st_iota = [.0]
# rt_iota = [.0]
# rt_lambda = [.2]
# rt_sigmah = [.225]
# combdict = Dict(:st_iota => st_iota,
#                 :st_lambda => [NaN],
#                 :st_sigmah => [NaN],
#                 :rt_iota => rt_iota,
#                 :rt_lambda => rt_lambda,
#                 :rt_sigmah => rt_sigmah)

# tmp = combination_row_finder(tcdf, combdict)
# row = 1
# df = tcdf[tcdf[:, :pair_num] .== tmp[row, :pair_num], :]

# ** Choose Result
df = tcdf[tcdf[:, :pair_num] .== pair_num, :]

# Set up Parameters
for x in names(df)
    st_rt_comb[x] = df[1, x]
end


# ** Create Joint Firms Objects
jfotc, jfep = JointEq.otc_ep_jfs(st_rt_comb[:k_otc], st_rt_comb[:k_ep];
                                 mu_s=st_rt_comb[:mu_s],
                                 st_iota=st_rt_comb[:st_iota],
                                 st_lambda=st_rt_comb[:st_lambda],
                                 st_sigmah=st_rt_comb[:st_sigmah],
                                 rt_iota=st_rt_comb[:rt_iota],
                                 rt_lambda=st_rt_comb[:rt_lambda],
                                 rt_sigmah=st_rt_comb[:rt_sigmah])

# * Make Directories and File Names
# Make results directories
jks_fpath = JointEq.make_jeq_jks_fpath(jfep)

# Full Information Equilibrium
fidf_fpath_name = string(jks_fpath, "/",
                         JointEqStructs.eq_type_dict["full_info"][:dfn], ".csv")

# Misrepresentation
misrepdf_fpath_name = string(jks_fpath, "/", JointEqStructs.eq_type_dict["misrep"][:dfn], ".csv")


# * Full Information - Bond Contract and Equilibrium Results
fidf, fieqdf = JointEq.get_fi_results(jfep, fidf_fpath_name;
                                      load_df=load_df,
                                      recompute_df=rerun_fi,
                                      save_df=save_df)


sf, rf, jks = JointEqStructs.form_firms_jks(jfep, fieqdf)

# * Misrepresentation
misrepdf = JointEq.get_misrep_results(jfep, fidf, misrepdf_fpath_name;
                                     load_df=load_df,
                                     recompute_df=rerun_misrep,
                                     save_df=save_df)

# * Prepare for Joint Equilibrium Computation
jeqid = JointEq.get_joint_eq_inputs(jfep, jks, jks_fpath;
                                fidf=fidf, fieqdf=fieqdf,
                                misrepdf=misrepdf)

# %%
