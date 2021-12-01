

using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
# using PyPlot
# using PyCall
# using Seaborn
# using LaTeXStrings

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "bond-model/modules/")
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "ModelPlots",
         "FullInfoEq", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


println(string("ARGUMENTS: ", ARGS))
# ################ SYS ARGUMENTS ################
# Capture Combination Number:
# cn = parse(Int, ARGS[1])
# rerun_fi = parse(Bool, ARGS[2])
# run_misrep = parse(Bool, ARGS[3])
# run_pool = parse(Bool, ARGS[4])
# run_sep = parse(Bool, ARGS[5])

cn = 1
rerun_fi = true
run_misrep = true
run_pool = false
run_sep = false


# INPUTS ###################################################
# Measure of Safe Firms
mu_s = .2

# Safe Firm's Risk Management Costs
sf_iota = 0.0 # 2.5 * 1e-4

# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :kappa  => [25 * 1e-4],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [x for x in Batch.cvm_param_values_dict[:iota] 
                                                  if .&(x >= sf_iota, x <= 20. * 1e-4)])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.2]
svmdict[:iota] = [.0]
svmdict[:sigmah] = Batch.svm_param_values_dict[:sigmah]
# #########################################################


# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)
# #########################################################


# Form Safe Firm ##########################################
sf_model = "cvm"
sf_df = (sf_model == "cvm") ? cvmdf : svmdf
sf_comb_num = sf_df[abs.(sf_df[:, :iota] .- sf_iota) .< 1e-5, :comb_num][1]
sf_bt, sf = Batch.get_bt_mobj(; model=sf_model, comb_num=sf_comb_num)
sf = ModelObj.set_opt_k_struct(sf, sf_df)
# #########################################################


# Form Risky Firm #########################################
rf_model = "cvm"
if cn + 1 <= size(cvmdf, 1)
    rf_comb_num = cvmdf[cn + 1, :comb_num]
    rf_bt, rf = Batch.get_bt_cvm(; comb_num=rf_comb_num)
    rf = ModelObj.set_opt_k_struct(rf, cvmdf)

    println(string("rf iota: ", rf.pm.iota))  
else #cn > size(cvmdf, 1)
    svm_row = cn + 1 - size(cvmdf, 1)
    rf_comb_num = svmdf[svm_row, :comb_num]
    rf_bt, rf = Batch.get_bt_svm(; comb_num=rf_comb_num)
    rf = ModelObj.set_opt_k_struct(rf, svmdf)
end
# #########################################################


# #########################################################
# Optimal p/c Ratio: sf.optKS.p/sf.optKS.c ~~ 12
ep_c = JointEq.round_value(sf.optKS.c)
ep_ratio = JointEq.round_value(sf.optKS.p/sf.optKS.c)
ep_p = ep_ratio * ep_c

# Pick Capital Structure
jks = JointEq.JointKStruct(mu_s, 
                           sf.optKS.mu_b,
                           sf.optKS.m, ep_c, ep_p, 
                           NaN, NaN, NaN, NaN, NaN)

# Joint Firm Object
jf = JointEq.JointFirms(jks, sf, rf,
                        sf_bt,
                        rf_bt,  
                        cvmdf, svmdf)
# #########################################################


# Make Directories and File Names ##########################
# Make results directories
jks_fpath = JointEq.make_jeq_jks_fpath(jf)

# Full Information Equilibrium
fi_fpath_name = JointEq.get_jeq_contour_fname(jks_fpath, cn; eq_type="full_info")

# Misrepresentation
misrep_fpath_name = JointEq.get_jeq_contour_fname(jks_fpath, cn; eq_type="misrep")

# Pooling Equilibrium
pool_fpath_name = JointEq.get_jeq_contour_fname(jks_fpath, cn; eq_type="pooling")

# Separating Equilibrium
sep_fpath_name = JointEq.get_jeq_contour_fname(jks_fpath, cn; eq_type="separating")
# #########################################################


# Joint Equilibrium Parameters ############################
jep = JointEq.store_joint_eq_parameters(mu_s, sf.pm.kappa, sf.pm.kappa;
                                        s_iota=sf.pm.iota,
                                        s_lambda=sf.pm.lambda,
                                        s_sigmah=sf.pm.sigmah,
                                        r_iota=rf.pm.iota,
                                        r_lambda=rf.pm.lambda,
                                        r_sigmah=rf.pm.sigmah)
# #########################################################


if !isfile(fi_fpath_name) 
    rerun_fi = true
end
if !isfile(misrep_fpath_name) 
    run_misrep = true
end


jeq = JointEq.ep_constructor(jep, jf.cvm_bt, jf.svm_bt;
                             ep_jks=jks,
                             run_misrep=run_misrep,
                             run_pool_eq=false,
                             run_sep_eq=false,                       
                             sf_obj_fun=:firm_value,
                             rf_obj_fun=:firm_value,
                             rerun_full_info=rerun_fi,
                             rerun_pool=false,
                             rerun_sep=false)


# Save Full Information Eq. Results
if rerun_fi
    fi_tmp = vcat(jeq.sfdf, jeq.rfdf)
    fi_tmp[!, :eq_type] .= "full_info"
    fi_tmp[!, :datetime] .= Dates.now()
    CSV.write(fi_fpath_name, fi_tmp)
    
    
    # Update
    fidf = CSV.read(string(jks_fpath, "/", JointEq.eq_type_dict["full_info"][:dfn], ".csv"))
    fidf = JointEq.remove_dup_save_df(vcat(fidf, fi_tmp), "full_info", jks_fpath)
end

# Save Misrepresentation Results
if run_misrep
    jeq.misrep[!, :datetime] .= Dates.now()
    CSV.write(misrep_fpath_name, jeq.misrep)
   
    # Update
    misrepdf_fpath_name = string(jks_fpath, "/", JointEq.eq_type_dict["misrep"][:dfn], ".csv")
    if isfile(misrepdf_fpath_name)
        # Load file
        misrepdf = CSV.read(misrepdf_fpath_name)

        # New Update
        # jeq.misrep[!, :sigmal] .= jeq.misrep[:, :s_sigmal]
        # cols = [x for x in names(misrepdf) if !(x in [:s_sigmal, :r_sigmal])]

        # Update
        misrepdf = JointEq.remove_dup_save_df(vcat(misrepdf[:, cols], jeq.misrep[:, cols]), "misrep", jks_fpath)
    else
        misrepdf = CSV.write(misrep_fpath_name, jeq.misrep)
    end
end


# ##########################
if !isfile(pool_fpath_name)
    run_pool = true
end

if !isfile(sep_fpath_name)
    run_sep = true
end

# #########################################################
# Compute Pooling Equilibrium #############################
jeq = JointEq.ep_constructor(jep, jf.cvm_bt, jf.svm_bt;
                             ep_jks=jks,
                             run_misrep=false,
                             run_pool_eq=run_pool,
                             run_sep_eq=run_sep,                       
                             sf_obj_fun=:firm_value,
                             rf_obj_fun=:MBR,
                             fi_fpath_name=fi_fpath_name,
                             rerun_full_info=false,
                             rerun_pool=false,
                             rerun_sep=false)
# #########################################################


# Save Results ############################################
# Save Pooling Eq. Results
if run_pool
    CSV.write(pool_fpath_name, jeq.pool)
end

# Save Separating Eq. Results
if run_sep
    CSV.write(sep_fpath_name, jeq.sep)
end
# ###########################################################

