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
    include(string(joinpath(module_path, modl), ".jl"))
end

plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_lambda"
rerun_misrep = false
save_misrepdf = true


# * Load Full Information Results
# Set Parameters ##########################################
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :kappa  => [25 * 1e-4],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [0.0, 1 * 1e-3])

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
cvmdf = deepcopy(DataFrame(cvmdf[2, :]))
# #########################################################

# * Compute Misrepresentation DF
# Misrepresentation Payoffs ###############################
script_dir = string(plot_script_path, "/", plots_xvar_dir)
misrepdf_fn = "misrepdf.csv"

if rerun_misrep | !(misrepdf_fn in readdir(script_dir))
    # Form Misrepresentation DF 
    LL = [ ]
    # Preliminary Objects #####################################
    sf_bt, sf = Batch.get_bt_mobj(; model=sf_model, comb_num=sf_comb_num)
    # sf_df = (sf_model == "cvm") ? cvmdf : svmdf
    sf_df = tmp
    sf = ModelObj.set_opt_k_struct(sf, sf_df)
    
    # Capital Structure -> Fixed
    jks = JointEq.JointKStruct(1., 
                               sf.optKS.mu_b,
                               sf.optKS.m, sf.optKS.c, sf.optKS.p, 
                               NaN, NaN, NaN, NaN, NaN)
    #  #########################################################
    
    for rf_comb_num in rf_comb_nums
        rf_bt, rf = Batch.get_bt_svm(; comb_num=rf_comb_num)
        
        # Joint Equilibrium Parameters
        jep = JointEq.store_joint_eq_parameters(jks.mu_s, sf.pm.kappa, sf.pm.kappa;
                                                s_iota=sf.pm.iota,
                                                s_lambda=sf.pm.lambda,
                                                s_sigmah=sf.pm.sigmah,
                                                r_iota=rf.pm.iota,
                                                r_lambda=rf.pm.lambda,
                                                r_sigmah=rf.pm.sigmah)
        
        # Compute Misrepresentation
        jeq = JointEq.ep_constructor(jep, sf_bt, rf_bt;
                                     ep_jks=jks,
                                     run_pool_eq=false,
                                     run_sep_eq=false,                       
                                     sf_obj_fun=firm_obj_fun,
                                     rf_obj_fun=firm_obj_fun,
                                     rerun_full_info=false,
                                     run_misrep=true,
                                     rerun_pool=false,
                                     rerun_sep=false)
        
        
        push!(LL, getfield(jeq, :misrep))
    end
    
    # Form Misrepresentation DataFrame
    misrepdf = vcat(LL...)


    cols = [:eq_deriv, :eq_min_val, :mu_b, :eq_deriv_min_val, 
            :eq_negative, :eq_vb, :MBR, :debt, :equity, :firm_value, 
            :leverage, :iota, :lambda, :sigmah, :delta, :sigmal, :obj_fun]
    tmp2 = DataFrame(misrepdf[1, :])
    for col in cols
        tmp2[!, Symbol(:r_, col)] .= tmp2[:, Symbol(:s_, col)]
    end
    tmp2[!, :rf_vb] .= tmp2[:, :sf_vb]
    tmp2[1, :r_lambda] = 0.0
    tmp2[!, :r_obj_fun] .= "misrep"
    misrepdf = vcat(tmp2, misrepdf)
    misrepdf[!, :s_sigmah] .= misrepdf[:, :s_sigmal]
    misrepdf[!, :s_lambda] .= 0.0
    # Save the DataFrame
    if save_misrepdf
        CSV.write(string(script_dir, "/", misrepdf_fn), misrepdf)
    end
else
    misrepdf = CSV.read(string(script_dir, "/", misrepdf_fn))#, types=JointEq.mps_col_types)
end
# #########################################################


# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    cvmdf[!, x] = cvmdf[:, x] .* 1e4
    svmdf[!, x] = svmdf[:, x] .* 1e4
end
misrepdf[!, :kappa] = misrepdf[:, :kappa] .* 1e4
for pf in [:s_, :r_]
    misrepdf[!, Symbol(pf, :iota)] = misrepdf[:, Symbol(pf, :iota)] .* 1e4
end
# #########################################################


# Cut-off Values ##########################################
xgrid = range(minimum(svmdf[:, :lambda]), stop=maximum(svmdf[:, :lambda]), length=10^5)

# sigmah : RM Firm Value = NRM Firm Value
fv_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                        :firm_value, cvmdf[1, :firm_value];
                                        xgrid=xgrid)

# sigmah : RM Firm Value = NRM Firm Value
mbr_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                         :MBR, cvmdf[1, :MBR];
                                         xgrid=xgrid)

# sigmah : FI MBR = Misrep MBR
cvm_misrep_lambda = 0.0
svm_misrep_lambda = ModelPlots.get_misrep_cutoff_value(:lambda, :MBR, 
                                                       deepcopy(cvmdf),
                                                       deepcopy(svmdf),
                                                       deepcopy(misrepdf),
                                                       xgrid=xgrid)
# #########################################################

# * Plots
# # Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value], 
                                   deepcopy(cvmdf), deepcopy(svmdf),
                                   interp_yvar=true,
                                   misrepdf=deepcopy(misrepdf),
                                   fv_xvar=fv_lambda,
                                   mbr_xvar=mbr_lambda,
                                   cvm_misrep_xvar=cvm_misrep_lambda,
                                   svm_misrep_xvar=svm_misrep_lambda,
                                   color_rm_region=false,
                                   color_nrm_region=false,
                                   color_conflict_region=false,
                                   color_misrep_region=true, 
                                   save_fig=true)


# # # Market-to-Book Ratio
mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:MBR], 
                                    deepcopy(cvmdf), deepcopy(svmdf),
                                    interp_yvar=true,
                                    misrepdf=deepcopy(misrepdf),
                                    fv_xvar=fv_lambda,
                                    mbr_xvar=mbr_lambda,
                                   cvm_misrep_xvar=cvm_misrep_lambda,
                                   svm_misrep_xvar=svm_misrep_lambda,
                                    color_rm_region=false,
                                    color_nrm_region=false,
                                    color_conflict_region=false,
                                    color_misrep_region=true, 
                                    save_fig=true)

# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value, :MBR], 
                                       deepcopy(cvmdf), deepcopy(svmdf),
                                       interp_yvar=true,
                                       misrepdf=deepcopy(misrepdf),
                                       fv_xvar=fv_lambda,
                                       mbr_xvar=mbr_lambda,
                                       cvm_misrep_xvar=cvm_misrep_lambda,
                                       svm_misrep_xvar=svm_misrep_lambda,
                                       color_rm_region=false,
                                       color_nrm_region=false,
                                       color_conflict_region=false,
                                       color_misrep_region=true, 
                                       save_fig=true)
