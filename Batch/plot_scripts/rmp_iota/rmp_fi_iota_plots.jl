main_path = "/home/artur/BondPricing"
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts")
plots_xvar_dir = "rmp_iota"
plot_script_name = "rmp_iota_preamble.jl"
include(string(plot_script_path, "/", plots_xvar_dir, "/", plot_script_name))


# Safe and Risky Firms' and Misrepresentation DFs #########
sfdf = deepcopy(cvmdf)
rfdf = deepcopy(svmdf)
# #########################################################

# Iota and Kappa in Basis Points
for x in [:iota, :kappa]
    sfdf[x] = sfdf[x] .* 1e4
    rfdf[x] = rfdf[x] .* 1e4
end

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
# #########################################################


# Firm Value
fv_fig = ModelPlots.rmp_fi_plotfun(:iota, [:firm_value], 
                                   deepcopy(sfdf), deepcopy(rfdf);
                                   xgrid=xgrid,
                                   fv_xvar=fv_iota)

# Market-to-Book Ratio
mbr_fig = ModelPlots.rmp_fi_plotfun(:iota, [:MBR], 
                                    deepcopy(sfdf), deepcopy(rfdf);
                                    xgrid=xgrid,
                                    fv_xvar=fv_iota, mbr_xvar=mbr_iota,
                                    color_conflict_region=true)

# Firm Value and MBR Multiplot
fv_mbr_fig = ModelPlots.rmp_fi_plotfun(:iota, [:firm_value, :MBR], 
                                       deepcopy(sfdf), deepcopy(rfdf);
                                       xgrid=xgrid,
                                       fv_xvar=fv_iota, mbr_xvar=mbr_iota,
                                       color_conflict_region=true)
