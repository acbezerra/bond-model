plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_lambda"
plot_script_name = "rmp_lambda_preamble.jl"
include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))


# Misrepresentation Payoffs ###############################
rerun_misrep = false
save_misrepdf = true
script_dir = string(plot_script_path, "/", plots_xvar_dir)
misrepdf_fn = "misrepdf.csv"
plot_script_name = "rmp_misrep_lambda_script.jl"

if rerun_misrep | !(misrepdf_fn in readdir(script_dir))
    # Form Misrepresentation DF 
    misrepdf = include(string(script_dir, "/", plot_script_name))
else
    misrepdf = CSV.read(string(script_dir, "/", misrepdf_fn), types=JointEq.mps_col_types)
end
# #########################################################


# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    cvmdf[x] = cvmdf[x] .* 1e4
    svmdf[x] = svmdf[x] .* 1e4
end
misrepdf[:kappa] = misrepdf[:kappa] .* 1e4
for pf in [:s_, :r_]
    misrepdf[Symbol(pf, :iota)] = misrepdf[Symbol(pf, :iota)] .* 1e4
end
# #########################################################


# Cut-off Values ##########################################
xgrid = range(minimum(svmdf[:lambda]), stop=maximum(svmdf[:lambda]), length=10^5)

# sigmah : RM Firm Value = NRM Firm Value
fv_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                        :firm_value, cvmdf[1, :firm_value];
                                        xgrid=xgrid)

# sigmah : RM Firm Value = NRM Firm Value
mbr_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                         :MBR, cvmdf[1, :MBR];
                                         xgrid=xgrid)

# sigmah : FI MBR = Misrep MBR
misrep_lambda = ModelPlots.get_misrep_cutoff_value(:lambda, :MBR, 
                                                   deepcopy(cvmdf),
                                                   deepcopy(svmdf),
                                                   deepcopy(misrepdf),
                                                   xgrid=xgrid)
# #########################################################

# return cvmdf, svmdf, misrepdf, fv_lambda, mbr_lambda, misrep_lambda

# Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value], 
                                   deepcopy(cvmdf), deepcopy(svmdf),
                                   interp_yvar=true,
                                   misrepdf=deepcopy(misrepdf),
                                   fv_xvar=fv_lambda,
                                   mbr_xvar=mbr_lambda,
                                   misrep_xvar=misrep_lambda,
                                   color_rm_region=false,
                                   color_nrm_region=false,
                                   color_conflict_region=false,
                                   color_misrep_region=true, 
                                   save_fig=true)


# # Market-to-Book Ratio
mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:MBR], 
                                    deepcopy(cvmdf), deepcopy(svmdf),
                                    interp_yvar=true,
                                    misrepdf=deepcopy(misrepdf),
                                    fv_xvar=fv_lambda,
                                    mbr_xvar=mbr_lambda,
                                    misrep_xvar=misrep_lambda,
                                    color_rm_region=false,
                                    color_nrm_region=false,
                                    color_conflict_region=false,
                                    color_misrep_region=true, 
                                    save_fig=true)

# # Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value, :MBR], 
                                       deepcopy(cvmdf), deepcopy(svmdf),
                                       interp_yvar=true,
                                       misrepdf=deepcopy(misrepdf),
                                       fv_xvar=fv_lambda,
                                       mbr_xvar=mbr_lambda,
                                       misrep_xvar=misrep_lambda,
                                       color_rm_region=false,
                                       color_nrm_region=false,
                                       color_conflict_region=false,
                                       color_misrep_region=true, 
                                       save_fig=true)
