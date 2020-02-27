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

# * SYS & Other Args
println(string("ARGUMENTS: ", ARGS))
# ################ SYS ARGUMENTS ################
# Position in the Safe Firm Measure array:
pair_num = parse(Int, ARGS[1])
rerun_fi = parse(Bool, ARGS[2])
rerun_misrep = parse(Bool, ARGS[3])
rerun_pool = parse(Bool, ARGS[4])
rerun_sep = parse(Bool, ARGS[5])

# pair_num=7
# rerun_fi=false
# rerun_misrep=false
# rerun_pool=false
# rerun_sep=false
# ##################################################

load_df=true
save_df=true

pool_obj_fun=:exp_firm_value
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
st_rt_comb[:mu_s] = .2


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
lambda_vec = Batch.svm_param_values_dict[:lambda]
sigmah_vec = Batch.svm_param_values_dict[:sigmah]

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
                                                        :rt_sigmah => sigmah_vec)
tcdf = FullInfoEq.get_types_comb_df(all_typesdict)

# ** Optional
function combination_row_finder(df::DataFrame, pardict::Dict{Symbol,Array{Float64,N} where N})
    cond  = [true for x in 1:size(tcdf, 1)]
    tester(df, x, y) = isnan(y) ? isnan.(df[:, x]) : Array(abs.(df[:, x] .- y ) .< 1e-5)
    
    for x in keys(pardict)
        bool_sum = sum([tester(df, x, y) for y in pardict[x]]) .> 0
        cond = .&(cond, bool_sum)
    end

    return DataFrame(df[cond, :])
end

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


# # * Misrepresentation
# misrepdf = JointEq.get_misrep_results(jfep, fidf, misrepdf_fpath_name;
#                                       load_df=load_df,
#                                       recompute_df=rerun_misrep,
#                                       save_df=save_df)


# # * Prepare for Joint Equilibrium Computation
# jeqid = JointEq.get_joint_eq_inputs(jfep, jks, jks_fpath;
#                                     fidf=fidf, fieqdf=fieqdf,
#                                     misrepdf=misrepdf)


# # * Pooling Equilibrium
# pooldf = JointEq.get_jeq_results(jfep, jks,
#                                  :pool, pool_obj_fun,
#                                  jeqid; load_df=load_df,
#                                  recompute_df=rerun_pool,
#                                  save_df=save_df)


# # * Separating Equilibrium
# sepdf = JointEq.get_jeq_results(jfep, jks,
#                                 :sep, sep_obj_fun,
#                                 jeqid; load_df=load_df,
#                                 recompute_df=rerun_sep,
#                                 save_df=save_df)


# * Old Code ===========================================
# ** Create Joint Firm Objects
# # INPUTS ###################################################
# k_otc = 50 * 1e-4
# k_ep = 25 * 1e-4

# # Measure of Safe Firms
# mu_s = .2

# # Safe Firm's Risk Management Costs
# st_iota = 0.0 # 2.5 * 1e-4

# firm_obj_fun = :firm_value


# # Transaction Costs and Volatility Risk Parameters
# pardict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
#                                         :m  => [1.],
#                                         :gross_delta => [0.02],
#                                         :kappa  => [k_ep, k_otc],
#                                         :mu_b => [1.0],
#                                         :xi => [1.0])

# # Risk Type Parameters
# iota_vec = [x for x in Batch.vm_param_values_dict[:iota]
#             if .&(x >= st_iota, x <= 20. * 1e-4)]
# lambda_vec = [.2]
# sigmah_vec = Batch.svm_param_values_dict[:sigmah]

# # Load Unconstrained Results
# cvmdf, svmdf = JointEqStructs.get_joint_types_dfs(pardict, iota_vec,
#                                                   lambda_vec, sigmah_vec;
#                                                   firm_obj_fun=firm_obj_fun)

# all_typesdict = Dict{Symbol, Array{Float64, N} where N}(:st_iota => [st_iota],
#                                                         :st_lambda => [NaN],
#                                                         :st_sigmah => [NaN],
#                                                         :rt_iota => iota_vec,
#                                                         :rt_lambda => lambda_vec,
#                                                         :rt_sigmah => sigmah_vec)
# tcdf = FullInfoEq.get_types_comb_df(all_typesdict)

# # My choice
# tc_pair = tcdf[.&(abs.(tcdf[:, :rt_iota] .- iota_vec[5]) .<1e-5,
#                   abs.(tcdf[:, :rt_sigmah] .- sigmah_vec[2]) .<1e-5), :pair_num][1]

# # Form types dictionary
# tmp = tcdf[abs.(tcdf[:, :pair_num] .- tc_pair) .< 1e-5, :]
# typesdict = Dict{Symbol, Float64}([x=> tmp[:, x][1] for x in names(tmp)])
# # #########################################################

# jfotc = JointEqStructs.form_joint_types(k_otc, mu_s, pardict, typesdict;
#                                         cvmdf=cvmdf, svmdf=svmdf,
#                                         firm_obj_fun=firm_obj_fun,
#                                         bc=nothing)

# jfep = JointEqStructs.form_joint_types(k_ep, mu_s, pardict, typesdict;
#                                        cvmdf=cvmdf, svmdf=svmdf,
#                                        firm_obj_fun=firm_obj_fun,
#                                        bc=nothing)
# *** Old Approach
# # INPUTS ###################################################
# # Measure of Safe Firms
# mu_s = .2

# # Safe Firm's Risk Management Costs
# st_rm_iota = 0.0 # 2.5 * 1e-4

# # Transaction Costs and Volatility Risk Parameters
# cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
#                                         :m  => [1.],
#                                         :gross_delta => [0.02],
#                                         :kappa  => [25 * 1e-4],
#                                         :mu_b => [1.0],
#                                         :xi => [1.0],
#                                         :iota => [x for x in Batch.cvm_param_values_dict[:iota]
#                                                   if .&(x >= st_rm_iota, x <= 20. * 1e-4)])
# svmdict = deepcopy(cvmdict)
# svmdict[:lambda] = [.2]
# svmdict[:iota] = [.0]
# svmdict[:sigmah] = Batch.svm_param_values_dict[:sigmah]
# # #########################################################

# # Get Safe and Risky Firms' Full Info Optimal Results #####
# firm_obj_fun = :firm_value
# cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
#                                              firm_obj_fun=firm_obj_fun)
# # #########################################################


# # Firm Types ##############################################
# # Safe Type - RMP-Contingent Combinations
# sf_rm_comb_num = cvmdf[abs.(cvmdf[:, :iota] .- st_rm_iota) .< 1e-5, :comb_num][1]
# sf_nrm_comb_num = 0

# # Risky Type - RMP-Contingent Combinations
# rf_iota = cvmdict[:iota][5]
# rf_lambda = svmdict[:lambda][1]
# rf_sigmah = svmdict[:sigmah][2]
# rf_rm_comb_num = cvmdf[abs.(cvmdf[:, :iota] .- rf_iota) .< 1e-5, :comb_num][1]
# rf_nrm_comb_num = svmdf[.&(abs.(svmdf[:, :lambda] .- rf_lambda) .< 1e-5,
#                            abs.(svmdf[:, :sigmah] .- rf_sigmah) .< 1e-5), :comb_num][1]

# # Form Types
# st = JointEq.firm_type_constructor(; cvm_comb_num=sf_rm_comb_num,
#                                    svm_comb_num=sf_nrm_comb_num,
#                                    set_rmpc_opt_k_struct=true,
#                                    cvmdf=cvmdf, svmdf=svmdf)

# rt = JointEq.firm_type_constructor(; cvm_comb_num=rf_rm_comb_num,
#                                    svm_comb_num=rf_nrm_comb_num,
#                                    set_rmpc_opt_k_struct=true,
#                                    cvmdf=cvmdf, svmdf=svmdf)
# # #########################################################

# # #########################################################
# # Optimal p/c Ratio: sf.optKS.p/sf.optKS.c ~~ 12
# st_opt_rmp = JointEq.get_otc_fi_opt_rmp(st)
# st_optKS = st_opt_rmp == :rm ? st.rm.fr.optKS : st.nrm.fr.optKS
# ep_c = JointEq.round_value(st_optKS.c)
# ep_ratio = JointEq.round_value(st_optKS.p/st_optKS.c)
# ep_p = ep_ratio * ep_c

# # # Set Capital Structure
# # jks = JointEq.JointKStruct(mu_s,
# #                            st_optKS.mu_b,
# #                            st_optKS.m, ep_c, ep_p)


# # tcp = JointEq.get_types_common_params(st)

# # # Joint Firm Object
# # jf = JointEq.JointFirms(jks, st, rt, tcp,
# #                         cvmdf, svmdf)
# # # #########################################################

# # # #########################################################


# # # Make Directories and File Names ##########################
# # # Make results directories

# ** Directories
# Make results directories
# jks_fpath = JointEq.make_jeq_jks_fpath(jfep)

# # Full Information Equilibrium
# fi_fpath_name = JointEq.get_jeq_contour_fname(jks_fpath, 1; eq_type="full_info")

# tmp = CSV.read(fi_fpath_name)


# ** Misrepresentation
# fpathname = string(pwd(), "/tmp.csv")
# # CSV.write(fpathname, ep_misrep_eqdf)
# ep_misrep_eqdf = CSV.read(fpathname)
# JointEq.reshape_joint_df(ep_misrep_eqdf)
#
# *** Old
# Safe Type Copies Risky Type #######

# # Opt Capital Structure
# rt_opt_mub = fieqdf[fieqdf[:, :type] .== :rt, :mu_b][1]
# opt_bc = jfep.bc

# # Find Optimal VB
# st_misrep_vb = FullInfoEq.find_full_info_vb(jfep.st.rm.fr, opt_bc, rt_opt_mub)

# # Compute Misrepresentation Payoffs
# st_misrep_df = EqFinDiff.eq_fd(jfep.st.rm.fr; vbl=st_misrep_vb, mu_b=rt_opt_mub,
#                 m=opt_bc.m, c=opt_bc.c, p=opt_bc.p)

# ** Joint Equilibria - Pooling
# sf, rf, jks = JointEqStructs.form_firms_jks(jfep, fieqdf)

# mu_bN = 20
# min_mu_b = .75 * minimum(fidf[:, :mu_b])
# max_mu_b = 1.25 * maximum(fidf[:, :mu_b])
# mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN)

# # Compute Optimal VB for each mu_b candidate
# jks2 = deepcopy(jks)

# spline_k = 3
# spline_bc = "extrapolate"
# N = 10^5

#resd = JointEq.find_joint_payoffs(sf, rf, jks2, [x for x in mu_b_grid])

# dfl = @time fetch(@spawn [JointEq.find_joint_optimal_vb(sf, rf, jks2;
#                                                     mu_b=mu_b,
#                                                     rerun_fi_vb=true)
#                               for mu_b in mu_b_grid])

# # Store results in a DataFrame
# dfl2 = vcat(dfl...)

# CSV.write(string(jks_fpath, "/dfl"), dfl2)

# df = CSV.read(string(jks_fpath, "/dfl"))

# # Form Refined mu_b grid
# ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=N)

# # Form Dictionary to Store Results
# resd = Dict{Symbol, Any}(:st => Dict{Symbol, Any}(),
#                          :rt => Dict{Symbol, Any}(),
#                          :mu_b_grid => ref_mu_b_grid)

# # Interpolate Objective Functions in mu_b ################################
# for ft in [:st, :rt]
#     aft = ft == :st ? :rt : :st
#     ftdf = df[isnan.(df[:, Symbol(aft, :_vb)]), :]
#     for fun in [:debt, :equity, :MBR, :firm_value]
#         resd[ft][fun] =  Dierckx.Spline1D(ftdf[:, :mu_b],
#                                           ftdf[:, fun];
#                                           k=spline_k, bc=spline_bc)
#     end
# end
# resd

# ** Joint Equilibria - Separating
# rf_objf = :MBR
# Notice here that the RMP might differ from the optimal
# RMP for the Risky Type under full information
# fi_rf_mu_b = fieqdf[fieqdf[:, :type] .== :rt, :mu_b][1]
# rf_ic_objf_val = fieqdf[fieqdf[:, :type] .== :rt, rf_objf][1]
# jeq_objf = :sf_firm_value
# N1 = 20

# # Need to add Bond Contract!
# bc = jfep.bc

# sepdf = JointEq.get_sep_eqdf(sf, rf, bc, jks,
#                              Array(mu_b_grid),
#                              fi_rf_mu_b,
#                              rf_ic_objf_val;
#                              jeq_objf=jeq_objf)


# ** Functions
# *** Get Misrep Payoffs
# function get_misrep_payoffs(jf, misrep_ft, fidf::DataFrame, fieqdf::DataFrame)
#     copied_ft = misrep_ft == :st ? :rt : :st

#     # Opt Capital Structure
#     opt_mub = fieqdf[fieqdf[:, :type] .== :copied_ft, :mu_b][1]
#     opt_bc = jf.bc

#     for rmp in [:rm, :nrm]
#         fr = getfield(getfield(jfep, misrep_ft), rmp).fr

#         misrep_df = FullInfoEq.get_empty_fidf()
#         if !isnothing(fr)
#             # Find Optimal VB
#             misrep_vb = FullInfoEq.find_full_info_vb(jfep.st.rm.fr, opt_bc, opt_mub)

#             # Compute Misrepresentation Payoffs
#             tmp = EqFinDiff.eq_fd(fr; vbl=misrep_vb, mu_b=opt_mub,
#                                         m=opt_bc.m, c=opt_bc.c, p=opt_bc.p)

#             tmp[!, :eq_type] .= :misrep
#             tmp[!, :datetime] .= Dates.now()
#             tmp[!, :type] .= misrep_ft
#             tmp[!, :rmp] .= rmp
#             tmp = vcat([misrep_df, tmp]...)
#         end

#         return misrep_df
#     end
# end

# *** Find Joint Payoffs
# function find_joint_payoffs(sf, rf, jks,
#                             mu_b_grid::Array{Float64,1};
#                             spline_k::Int64=3,
#                             spline_bc::String="extrapolate",
#                             N::Int64=10^5)


#     dfl = @time fetch(@spawn [JointEq.find_joint_optimal_vb(sf, rf, jks;
#                                                     mu_b=mu_b,
#                                                     rerun_fi_vb=true)
#                               for mu_b in mu_b_grid])

#     # Store results in a DataFrame
#     df = vcat(dfl...)

#     # Form Refined mu_b grid
#     ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=N)

#     # Form Dictionary to Store Results
#     res = Dict{Symbol, Any}(:st => Dict{Symbol, Any}(),
#                             :rt => Dict{Symbol, Any}(),
#                             :mu_b_grid => ref_mu_b_grid)

#     # Interpolate Objective Functions in mu_b ################################
#     for ft in [:st, :rt]
#         aft = ft == :st ? :rt : :st
#         ftdf = df[isnan.(df[:, Symbol(aft, :_vb)]), :]
#         for fun in [:debt, :equity]
#             res[ft][fun] =  Dierckx.Spline1D(ftdf[:, :mu_b],
#                                              ftdf[:, fun];
#                                              k=spline_k, bc=spline_bc)
#         end
#     end

#     return res
# end
#
# *** Form Objective Function
# function form_obj_fun(resd::Dict{Symbol, Any}, obj_fun::Symbol,
#                       ftype::Symbol, mu_s::Float64, V0::Float64)
#    ft_fv(ft, x) = resd[ft][:debt](x) + resd[ft][:equity](x)
#    exp_fv(x) = mu_s * ft_fv(:st, x) + (1 - mu_s) * ft_fv(:rt, x)
#    ft_mbr(ft, x) = resd[:ft][:equity](x) / (V0 - resd[:ft][:debt](x))

#     if ftype != :joint
#         if obj_fun == :firm_value
#             objf(x) = ft_fv(ftype, x)
#             return objf
#         elseif obj_fun == :mbr
#             return ft_mbr(ftype, x)
#         end
#     else
#         return exp_fv
#     end
# end
# st_fv = form_obj_fun(resd, :firm_value, :st, .2, 100.)
# st_fv(.5)

# ** Other
# Form Grid of mu_b candidates
# min_mu_b = .75 * minimum([fi_sf_mu_b, fi_rf_mu_b])
# max_mu_b = 1.25 * maximum([fi_sf_mu_b, fi_rf_mu_b])
# mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN)

# opt_st_vb, opt_rt_vb = JointEq.interp_optimal_vbs(misrep_jks, sfdf)
# st_vb, rt_vb = JointEq.get_type_contingent_vbs(opt_st_vb,
#                                            misrep_jks.fi_st_vb,
#                                            misrep_jks.fi_rt_vb;
#                                            sf_defaults_first=true)
# s1 = JointEq.joint_eq_fd(sf, rf, misrep_jks, st_vb=st_vb, rt_vb=rt_vb)
# st_vb, rt_vb = JointEq.get_type_contingent_vbs(opt_rt_vb,
#                                            misrep_jks.fi_st_vb,
#                                            misrep_jks.fi_rt_vb;
#                                            sf_defaults_first=true)
# s2 = JointEq.joint_eq_fd(sf, rf, misrep_jks, st_vb=st_vb, rt_vb=rt_vb)

# opt_st_vb, opt_rt_vb = JointEq.interp_optimal_vbs(misrep_jks, rfdf)

# # vbl that sets E'_s(vbl) to zero
# st_vb, rt_vb = JointEq.get_type_contingent_vbs(opt_st_vb,
#                                        misrep_jks.fi_st_vb,
#                                        misrep_jks.fi_rt_vb;
#                                        sf_defaults_first=false)
# r1 = JointEq.joint_eq_fd(sf, rf, misrep_jks; st_vb=st_vb, rt_vb=rt_vb)

# # vbl that sets E'_r(vbl) to zero
# st_vb, rt_vb = JointEq.get_type_contingent_vbs(opt_rt_vb,
#                                        misrep_jks.fi_st_vb,
#                                        misrep_jks.fi_rt_vb;
#                                        sf_defaults_first=false)
# r2 = JointEq.joint_eq_fd(sf, rf, misrep_jks; st_vb=st_vb, rt_vb=rt_vb)
# if abs.(r2[isnan.(r2[:, :st_vb]), :eq_deriv][1]) > 1e-2
#     r2 = JointEq.refine_contingent_vbs(sf, rf, misrep_jks, st_vb, rt_vb)
# end

# ** Optimal Mu_b
# Form Grid of mu_b candidates
# min_mu_b = .75 * minimum([fi_sf_mu_b, fi_rf_mu_b])
# max_mu_b = 1.25 * maximum([fi_sf_mu_b, fi_rf_mu_b])
# mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN)

# rf_objf = :MBR
# rf_ic_objf_val = fieqdf[fieqdf[:, :type] .== :rt, rf_objf][1]
# jeq_objf = :exp_firm_value
# N1 = 20

# # The Risky-Type's Incentive Compatibility Condition
# # Range of mu_b values in pooling equilibrium
# rt_ic_cond = resd[:rt][rf_objf](resd[:mu_b_grid]) .>= rf_ic_objf_val
# filtered_mu_b_grid_1, filtered_mu_b_grid_2 = JointEq.find_mu_b_intervals(resd[:mu_b_grid],
#                                                                  rt_ic_cond; N=N1)

# # Form Objective Function
# objf = JointEq.get_joint_objf(resd, jks, jeq_objf)

# # Maximize the Objective Function s.t.
# # the Incentive Compatibility Constraint
# mu_b_opt = JointEq.find_opt_mu_b(objf,
#                                  Array(filtered_mu_b_grid_1),
#                                  Array(filtered_mu_b_grid_2))

# #     # Compute Results
# tmp = JointEq.find_joint_optimal_vb(sf, rf, jks;
#                                  mu_b=mu_b_opt, rerun_fi_vb=true)
# ** Joint Payoffs
# resd = JointEq.find_joint_payoffs(sf, rf, jks2, [x for x in mu_b_grid])
#
# ** Optimal VB
# Get Candidates
# df = JointEq.compile_opt_vb_results(sf, rf, misrep_jks, sfdf, rfdf)
