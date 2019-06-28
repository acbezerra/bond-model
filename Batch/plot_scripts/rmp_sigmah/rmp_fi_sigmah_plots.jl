main_path = "/home/artur/BondPricing"
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_sigmah"
plot_script_name = "rmp_sigmah_preamble.jl"
include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))

# Iota and Kappa in Basis Points ##########################
for x in [:iota, :kappa]
    cvmdf[x] = cvmdf[x] .* 1e4
    svmdf[x] = svmdf[x] .* 1e4
end
# #########################################################


# Cut-off Values ##########################################
xgrid = range(minimum(svmdf[:sigmah]), stop=maximum(svmdf[:sigmah]), length=10^5)

# sigmah : RM Firm Value = NRM Firm Value
fv_sigmah = ModelPlots.get_cutoff_value(svmdf, :sigmah,
                                        :firm_value, cvmdf[1, :firm_value];
                                        xgrid=xgrid)

# sigmah : RM Firm Value = NRM Firm Value
mbr_sigmah = ModelPlots.get_cutoff_value(svmdf, :sigmah,
                                         :MBR, cvmdf[1, :MBR];
                                         xgrid=xgrid)
# #########################################################

# Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:sigmah, [:firm_value], 
                          deepcopy(cvmdf), deepcopy(svmdf),
                          fv_xvar=fv_sigmah,
                          color_rm_region=true,
                          color_nrm_region=true,
                          color_conflict_region=false,
                          color_misrep_region=false, 
                          save_fig=true)

# Market-to-Book Ratio
mbr_fig = ModelPlots.rmp_fi_plotfun(:sigmah, [:MBR], 
                                    deepcopy(cvmdf), deepcopy(svmdf),
                                    fv_xvar=fv_sigmah,
                                    mbr_xvar=mbr_sigmah,
                                    color_rm_region=true,
                                    color_nrm_region=true,
                                    color_conflict_region=true,
                                    color_misrep_region=false, 
                                    save_fig=true)



# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:sigmah, [:firm_value, :MBR], 
                                       deepcopy(cvmdf), deepcopy(svmdf),
                                       fv_xvar=fv_sigmah,
                                       mbr_xvar=mbr_sigmah,
                                       color_rm_region=true,
                                       color_nrm_region=true,
                                       color_conflict_region=true,
                                       color_misrep_region=false, 
                                       save_fig=true)
