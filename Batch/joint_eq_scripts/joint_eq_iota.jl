
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
module_path = string(main_path, "/", "Julia/modules/")
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "ModelPlots", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end


println(string("ARGUMENTS: ", ARGS))
# # ################ SYS ARGUMENTS ################
# # Position in the Safe Firm Measure array:
mu_s_pos = parse(Int, ARGS[1])
rerun_fi = parse(Bool, ARGS[2])
run_misrep = parse(Bool, ARGS[3])
run_pool = parse(Bool, ARGS[4])
run_sep = parse(Bool, ARGS[5])


# INPUTS ###################################################
# Measure of Safe Firms
mu_s_grid = range(.05, stop=.95, length=10)
mu_s = mu_s_grid[mu_s_pos]

# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                               :m  => [1.],
                               :gross_delta => [0.02],
                               :kappa  => [25 * 1e-4],
                               :mu_b => [1.0],
                               :xi => [1.0],
                               :iota => [x for x in Batch.cvm_param_values_dict[:iota] 
                                         if x <= 20. * 1e-4])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.2]
svmdict[:iota] = [.0]
svmdict[:sigmah] = [.225]
# #########################################################


# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)
# #########################################################


# Form Safe Firm ##########################################
sf_model = "cvm"
sf_iota = 2.5 * 1e-4
sf_df = (sf_model == "cvm") ? cvmdf : svmdf
sf_comb_num = sf_df[abs.(sf_df[:iota] .- sf_iota) .< 1e-5, :comb_num][1]
sf_bt, sf = Batch.get_bt_mobj(; model=sf_model, comb_num=sf_comb_num)
sf = ModelObj.set_opt_k_struct(sf, sf_df)
# #########################################################


# Form Risky Firm #########################################
rf_model = "svm"
rf_comb_num = svmdf[1, :comb_num]
rf_bt, rf = Batch.get_bt_svm(; comb_num=rf_comb_num)
rf_df = (rf_model == "cvm") ? cvmdf : svmdf
rf = ModelObj.set_opt_k_struct(rf, rf_df)
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
fi_fpath_name = JointEq.get_jeq_mus_fname(jks_fpath; eq_type="full_info")

# Misrepresentation
misrep_fpath_name = JointEq.get_jeq_mus_fname(jks_fpath; eq_type="misrep")

# Pooling Equilibrium
pool_fpath_name = JointEq.get_jeq_mus_fname(jks_fpath; eq_type="pooling", mu_s=jf.jks.mu_s)

# Separating Equilibrium
sep_fpath_name = JointEq.get_jeq_mus_fname(jks_fpath; eq_type="separating")
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

    
# Full Information and Misrepresentation Results ##########
# Wait to load Full Information Results
if mu_s_pos > 1
    if isfile(fi_fpath_name)
        mod_time = unix2datetime(stat(fidf_file).mtime)
        minute_diff = (Dates.now() - mod_time).value/10^3/60
        
        if minute_diff > 10
            rerun_fi = true
        end
    else
        tic = time()
        while(.&(!isfile(fi_fpath_name), (time() -tic)/60 < 1))
            sleep(5)
        end
                
        if !isfile(fi_fpath_name)
            rerun_fi = true
        end
    end
end

        

if any([mu_s_pos == 1, rerun_fi])
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
        fi_tmp[:eq_type] = "full_info"
        fi_tmp[:datetime] = Dates.now()
        CSV.write(fi_fpath_name, fi_tmp)


        # Update
        fidf = CSV.read(string(jks_fpath, "/", JointEq.eq_type_dict["full_info"][:dfn], ".csv"))
        fidf = JointEq.remove_dup_save_df(vcat(fidf, fi_tmp), "full_info", jks_fpath)
    end
    
    # Save Misrepresentation Results
    if run_misrep
        jeq.misrep[:datetime] = Dates.now()
        CSV.write(misrep_fpath_name, jeq.misrep)

        # Update
        misrepdf = CSV.read(string(jks_fpath, "/", JointEq.eq_type_dict["misrep"][:dfn], ".csv"))
        misrepdf = JointEq.remove_dup_save_df(vcat(misrepdf, jeq.misrep), "misrep", jks_fpath)
    end
end    
# ##########################################################


# Pooling and Separating Equilibria ########################
skip_cond = .&(!run_pool, !run_sep)

if !skip_cond
    if !isfile(pool_fpath_name)
        run_pool = true
    end
    
    if .&(!isfile(sep_fpath_name), mu_s_pos == size(mu_s_grid, 1))
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


    # Update Results DataFrames ################################# 
    if mu_s_pos == size(mu_s_grid, 1)
        pool_list = JointEq.collect_joint_eq_files(jks_fpath; eq_type="pooling")
        sep_list = JointEq.collect_joint_eq_files(jks_fpath; eq_type="separating")
        # files_num = minimum([size(pool_list,1), size(sep_list,1)])

        tic = time()
        while .&(size(pool_list,1) < mu_s_pos, size(sep_list, 1) < 1, (time() -tic)/60 < 10)
            sleep(5)
            pool_list = JointEq.collect_joint_eq_files(jks_fpath; eq_type="pooling")
            sep_list = JointEq.collect_joint_eq_files(jks_fpath; eq_type="separating")
            files_num = minimum([size(pool_list,1), size(sep_list,1)])
        end
        
        # Collect files
        pooldf, sepdf = JointEq.joint_eq_form_dataframes(; pool_list=pool_list, 
                                                           sep_list=sep_list)
        
        pooldf_all = JointEq.remove_dup_save_df(pooldf, "pooling", jks_fpath)
        sepdf_all = JointEq.remove_dup_save_df(sepdf, "separating", jks_fpath)
    end
    # ###########################################################

    
    return jeq
end

# return jeq
