plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_lambda"
plot_script_name = "rmp_lambda_preamble.jl"
include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))

println(svmdf)

# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    cvmdf[x] = cvmdf[x] .* 1e4
    svmdf[x] = svmdf[x] .* 1e4
end
# #########################################################


# Cut-off Values ##########################################
xgrid = range(minimum(svmdf[:lambda]), stop=maximum(svmdf[:lambda]), length=10^5)

# lambda : RM Firm Value = NRM Firm Value
fv_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                        :firm_value, cvmdf[1, :firm_value];
                                        xgrid=xgrid)

# lambda : RM Firm Value = NRM Firm Value
mbr_lambda = ModelPlots.get_cutoff_value(svmdf, :lambda,
                                         :MBR, cvmdf[1, :MBR];
                                         xgrid=xgrid)
# #########################################################

# return cvmdf, svmdf, fv_lambda, mbr_lambda

# Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value], 
                          deepcopy(cvmdf), deepcopy(svmdf),
                          fv_xvar=fv_lambda,
                          color_rm_region=true,
                          color_nrm_region=true,
                          color_conflict_region=false,
                          color_misrep_region=false, 
                          save_fig=true)

# Market-to-Book Ratio
mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:MBR], 
                                    deepcopy(cvmdf), deepcopy(svmdf),
                                    fv_xvar=fv_lambda,
                                    mbr_xvar=mbr_lambda,
                                    color_rm_region=true,
                                    color_nrm_region=true,
                                    color_conflict_region=true,
                                    color_misrep_region=false, 
                                    save_fig=true)



# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:lambda, [:firm_value, :MBR], 
                                       deepcopy(cvmdf), deepcopy(svmdf),
                                       fv_xvar=fv_lambda,
                                       mbr_xvar=mbr_lambda,
                                       color_rm_region=true,
                                       color_nrm_region=true,
                                       color_conflict_region=true,
                                       color_misrep_region=false, 
                                       save_fig=true)
