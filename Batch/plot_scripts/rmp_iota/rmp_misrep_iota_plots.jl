main_path = "/home/artur/BondPricing"
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_iota"
plot_script_name = "rmp_iota_preamble.jl"
include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))


# Set Targeted Safe Firm
sf_model = "cvm"
sf_iota = 2.5 * 1e-4
sf_comb_num = cvmdf[abs.(cvmdf[:iota] .- sf_iota) .< 1e-5, :comb_num][1]
# #########################################################


# Misrepresentation Payoffs ###############################
rerun_misrep = false
save_misrepdf = true
script_dir = string(plot_script_path, "/", plots_xvar_dir)
misrepdf_fn = "misrepdf.csv"
plot_script_name = "rmp_misrep_iota_script.jl"

if rerun_misrep | !(misrepdf_fn in readdir(script_dir))
    # Form Misrepresentation DF 
    misrepdf = include(string(script_dir, "/", plot_script_name))
else
    misrepdf = CSV.read(string(script_dir, "/", misrepdf_fn))#, types=JointEq.mps_col_types)
end
# misrepdf = include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))
# #########################################################


# Safe and Risky Firms' and Misrepresentation DFs #########
sfdf = deepcopy(cvmdf)
rfdf = deepcopy(svmdf)
# #########################################################


# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    sfdf[x] = sfdf[x] .* 1e4
    rfdf[x] = rfdf[x] .* 1e4
end
misrepdf[:kappa] = misrepdf[:kappa] .* 1e4
for pf in [:s_, :r_]
    misrepdf[Symbol(pf, :iota)] = misrepdf[Symbol(pf, :iota)] .* 1e4
end
# #########################################################


# Cut-off Values ##########################################
xgrid = range(minimum(sfdf[:iota]), stop=maximum(sfdf[:iota]), length=10^5)
# iota : RM Firm Value = NRM Firm Value
fv_iota = ModelPlots.get_cutoff_value(sfdf, :iota,
                                      :firm_value, rfdf[1, :firm_value];
                                      xgrid=xgrid)

# iota : RM Firm Value = NRM Firm Value
mbr_iota = ModelPlots.get_cutoff_value(sfdf, :iota,
                                       :MBR, rfdf[1, :MBR];
                                       xgrid=xgrid)

# iota : FI MBR = Misrep MBR
cvm_misrep_iota = ModelPlots.get_cutoff_value(sfdf, :iota, :MBR, misrepdf[1, :r_MBR];
                                          xgrid=xgrid)
# #########################################################


# Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:firm_value], 
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

# Market-to-Book Ratio
fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:MBR], 
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



