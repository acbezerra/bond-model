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
plots_xvar_dir = "rmp_sigmah"
rerun_misrep = false
save_misrepdf = true

# * Load Full Information Results
# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [.15, .2, .25, .30],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :kappa  => [25 * 1e-4],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [x for x in Batch.cvm_param_values_dict[:iota] 
                                                  if x <= 20. * 1e-4])

# anything goes here. For the CVM, we don't
# use the SVM df. 
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.2]

# #########################################################
# Get Safe and Risky Firms' Full Info Optimal Results #####
firm_obj_fun = :firm_value 

# #########################################################
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)
                
cvmdf[!, :sigmah] = cvmdf[!, :sigmal]


# Set Targeted Safe Firm
sf_model = "cvm"
sf_cond = .&(abs.(cvmdf[:, :iota] .- 10 * 1e-4) .< 1e-5,
             abs.(cvmdf[:, :sigmal] .- .15) .< 1e-5)      
sf_comb_num = cvmdf[sf_cond, :comb_num][1]

rf_model = "cvm"
rf_cond = abs.(cvmdf[:, :iota] .- .0) .< 1e-5
rf_comb_nums = cvmdf[rf_cond, :comb_num]
# #########################################################


# * Compute Misrepresentation DF
# Misrepresentation Payoffs ###############################
script_dir = string(plot_script_path, "/", plots_xvar_dir)
misrepdf_fn = "cvm_misrepdf.csv"

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
    #  #########################################################
    
    for rf_comb_num in rf_comb_nums
        rf_bt, rf = Batch.get_bt_cvm(; comb_num=rf_comb_num)
        
        # Joint Equilibrium Parameters
        jep = JointEq.store_joint_eq_parameters(jks.mu_s, sf.pm.kappa, sf.pm.kappa;
                                                s_iota=sf.pm.iota,
                                                s_lambda=sf.pm.lambda,
                                                s_sigmal=sf.pm.sigmal,
                                                s_sigmah=sf.pm.sigmah,
                                                r_iota=rf.pm.iota,
                                                r_sigmal=rf.pm.sigmal,
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

    # Save the DataFrame
    if save_misrepdf
        CSV.write(string(script_dir, "/", misrepdf_fn), misrepdf)
    end
else
    misrepdf = CSV.read(string(script_dir, "/", misrepdf_fn)) #, types=JointEq.mps_col_types)
end
# #########################################################

# * Adjustments
# ** Form Safe and Risky Firms' DataFrames
sfdf = cvmdf[cvmdf[:, :comb_num] .== sf_comb_num, :]
rfdf_cond = [in(x, rf_comb_nums) for x in cvmdf[:, :comb_num]]
rfdf = cvmdf[rfdf_cond, :]

# ** Make CVM results work with EP constructor
misrepdf[!, :s_sigmah] = misrepdf[:, :s_sigmal]
misrepdf[!, :r_sigmah] = misrepdf[:, :r_sigmal]

# ** Set Iota and Kappa variables to basis points
# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    sfdf[!, x] = sfdf[:, x] .* 1e4
    rfdf[!, x] = rfdf[:, x] .* 1e4
end
misrepdf[!, :kappa] = misrepdf[:, :kappa] .* 1e4
for pf in [:s_, :r_]
    misrepdf[!, Symbol(pf, :iota)] = misrepdf[:, Symbol(pf, :iota)] .* 1e4
end
#########################################################

# * Find Cut-off values
# Cut-off Values ##########################################
xgrid = range(minimum(rfdf[:, :sigmah]), stop=maximum(rfdf[:, :sigmah]), length=10^5)

# sigmah : RM Firm Value = NRM Firm Value
fv_sigmah = ModelPlots.get_cutoff_value(rfdf, :sigmah,
                                        :firm_value, sfdf[1, :firm_value];
                                        xgrid=xgrid)

# sigmah : RM Firm Value = NRM Firm Value
mbr_sigmah = ModelPlots.get_cutoff_value(rfdf, :sigmah,
                                         :MBR, sfdf[1, :MBR];
                                         xgrid=xgrid)

# sigmah : FI MBR = Misrep MBR
cvm_misrep_sigmah, svm_misrep_sigmah = ModelPlots.get_misrep_cutoff_value(:sigmah, :MBR, 
                                                                          deepcopy(sfdf),
                                                                          deepcopy(rfdf),
                                                                          deepcopy(misrepdf),
                                                                          xgrid=xgrid)
xvar = :sigmah
yvar = :MBR

xvals = rfdf[:, xvar]
if size(xgrid, 1) == 0
    xgrid = range(minimum(xvals), stop=maximum(xvals), length=10^5)
end

misrep_yval_interp = Dierckx.Spline1D(misrepdf[:, Symbol(:r_, xvar)], 
                                      misrepdf[:, Symbol(:r_, yvar)]; 
                                      k=3, bc="extrapolate")
fi_fv_interp = Dierckx.Spline1D(xvals, rfdf[:, :firm_value]; k=3, bc="extrapolate")
fi_svm_mbr_interp = Dierckx.Spline1D(xvals, rfdf[:, :MBR]; k=3, bc="extrapolate")
fi_yval_interp = [(sfdf[1, :firm_value] > fi_fv_interp(x)) ? sfdf[1, :MBR] : fi_svm_mbr_interp(x)
                  for x in xgrid]


svm_cv = NaN
svm_diff = abs.(fi_yval_interp .- misrep_yval_interp(xgrid))
if !isempty(svm_diff .< 1e-5)
    cv1 = xgrid[argmin(svm_diff)]
    cv2 = xgrid[argmax(svm_diff)]
end
# #########################################################


# * Plots
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:sigmah, [:firm_value, :MBR], 
                                       deepcopy(sfdf), deepcopy(rfdf),
                                       interp_yvar=true,
                                       misrepdf=deepcopy(misrepdf),
                                       fv_xvar=fv_sigmah,
                                       mbr_xvar=mbr_sigmah,
                                       cvm_misrep_xvar=cv1,
                                       svm_misrep_xvar=fv_sigmah,
                                       color_rm_region=false,
                                       color_nrm_region=false,
                                       color_conflict_region=false,
                                       color_misrep_region=true, 
                                       save_fig=true)

