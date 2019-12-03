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
plots_xvar_dir = "rmp_iota"
rerun_misrep =  true
save_misrepdf = true

# * Load Full Information Results
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
# svmdict[:kappa] = [25 * 1e-4]
svmdict[:iota] = [.0]
svmdict[:sigmah] = [.225]

# #########################################################
# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)



# Set Targeted Safe Firm
sf_model = "cvm"
sf_iota = 0.0 # 2.5 * 1e-4
sf_comb_num = cvmdf[abs.(cvmdf[:, :iota] .- sf_iota) .< 1e-5, :comb_num][1]
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
    sf_df = (sf_model == "cvm") ? cvmdf : svmdf
    sf = ModelObj.set_opt_k_struct(sf, sf_df)
    
    # Capital Structure -> Fixed
    jks = JointEq.JointKStruct(1., 
                               sf.optKS.mu_b,
                               sf.optKS.m, sf.optKS.c, sf.optKS.p, 
                               NaN, NaN, NaN, NaN, NaN)
    
    rf_bt, rf = Batch.get_bt_svm(; comb_num=svmdf[1, :comb_num])
    
    # Joint Firm Object
    jf = JointEq.JointFirms(jks, sf, rf,
                            sf_bt,
                            rf_bt,  
                            cvmdf, svmdf)    
    # #########################################################
    
    for rf_comb_num in svmdf[:, :comb_num] #svm_combs
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
        jeq = JointEq.ep_constructor(jep, jf.cvm_bt, jf.svm_bt;
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

    # Save the DataFrame
    if save_misrepdf
        CSV.write(string(script_dir, "/", misrepdf_fn), misrepdf)
    end
else
    misrepdf = CSV.read(string(script_dir, "/", misrepdf_fn))#, types=JointEq.mps_col_types)
end
# #########################################################

# * Adjustments
# ** Form Safe and Risky Firms' DataFrames
# Safe and Risky Firms' and Misrepresentation DFs #########
sfdf = deepcopy(cvmdf)
rfdf = deepcopy(svmdf)
# #########################################################

# ** Set Iota and Kappa variables to basis points
# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    sfdf[!, x] = sfdf[:, x] .* 1e4
    rfdf[!, x] = rfdf[:, x] .* 1e4
end
misrepdf[!, :kappa] = misrepdf[!, :kappa] .* 1e4
for pf in [:s_, :r_]
    misrepdf[!, Symbol(pf, :iota)] = misrepdf[:, Symbol(pf, :iota)] .* 1e4
end
# #########################################################

# * Find Cut-off values
# Cut-off Values ##########################################
xgrid = range(minimum(sfdf[:, :iota]), stop=maximum(sfdf[:, :iota]), length=10^5)
# iota : RM Firm Value = NRM Firm Value
fv_iota = ModelPlots.get_cutoff_value(sfdf, :iota,
                                      :firm_value, rfdf[1, :firm_value];
                                      xgrid=xgrid)

# iota : RM Firm Value = NRM Firm Value
mbr_iota = ModelPlots.get_cutoff_value(sfdf, :iota,
                                       :MBR, rfdf[1, :MBR];
                                       xgrid=xgrid)

# iota : FI MBR = Misrep MBR
# iota : min(iota) s.t. Misrep-MBR(iota) >= FI-MBR(iota)
cvm_misrep_iota = fv_iota # ModelPlots.get_cutoff_value(sfdf, :iota, :MBR, misrepdf[1, :r_MBR];
                          #                xgrid=xgrid)
# #########################################################

# * Plots
# Firm Value
# fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:firm_value], 
#                                    deepcopy(sfdf), deepcopy(rfdf),
#                                    interp_yvar=true,
#                                    misrepdf=deepcopy(misrepdf),
#                                    fv_xvar=fv_iota,
#                                    cvm_misrep_xvar=cvm_misrep_iota,
#                                    svm_misrep_xvar=NaN,
#                                    color_rm_region=false,
#                                    color_nrm_region=false,
#                                    color_conflict_region=false,
#                                    color_misrep_region=true)

# # Market-to-Book Ratio
# fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:MBR], 
#                                    deepcopy(sfdf), deepcopy(rfdf),
#                                    interp_yvar=true,
#                                    misrepdf=deepcopy(misrepdf),
#                                    fv_xvar=fv_iota,
#                                    cvm_misrep_xvar=cvm_misrep_iota,
#                                    svm_misrep_xvar=NaN,
#                                    color_rm_region=false,
#                                    color_nrm_region=false,
#                                    color_conflict_region=false,
#                                    color_misrep_region=true)

# Firm Value and MBR Multiplot
fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:firm_value, :MBR], 
                                   deepcopy(sfdf), deepcopy(rfdf),
                                   interp_yvar=true,
                                   misrepdf=deepcopy(misrepdf),
                                   fv_xvar=fv_iota,
                                   cvm_misrep_xvar=cvm_misrep_iota,
                                   svm_misrep_xvar=NaN,
                                   color_rm_region=false,
                                   color_nrm_region=false,
                                   color_conflict_region=false,
                                   color_misrep_region=true)



