
# * Directories & File Names ###############################################
main_dir_path = form_main_dir_path(main_dir)
plots_dir="Plots"
contour_plots_path = string(main_dir_path, "/Plots/Contour")

cvm_vs_svm_plots_dir = "CVMvsSVM"
rmp_plots_dir = "RMP"
heat_surf_graph_dir="HeatSurf"
obj_fun_dict = Dict{Symbol, String}(:firm_value => "fv",
                                    :MBR => "mbr")
# ########################################################################


# * Auxiliary Functions ####################################################
function par_val_printer(x::Symbol)
    return !(x in [:iota, :kappa]) ? string(x, "_", xylabels[x][2]) : string(x, "_bp_", xylabels[x][2])
end

function par_val_adj(x::Symbol, val::Float64)
    return !(x in [:iota, :kappa]) ? val : val * 1e4
end
# ########################################################################


vartitles = Dict{Symbol, String}(:vb => "\$V^B\$",
                                 :c => "C",
                                 :p => "\$P = Debt\$",
                                 :equity => "Equity",
                                 :firm_value => "Firm Value",
                                 :leverage => "Leverage (\$\\%\$)",
                                 :MBR => "MBR (\$\\%\$)")


# * CVM Plots ##############################################################
cvm_plots_title_params_order = [:mu_b, :m, :iota, :xi, :kappa, :sigmal]


# * SVM HeatMap and Surface Plots ##########################################
xylabels = Dict{Symbol, Array{String,1}}(:mu_b => ["\\mu_b", comb_folder_dict[:mu_b][2]],
                                         :m => ["m", comb_folder_dict[:m][2]],
                                         :iota => ["\\iota \\, (b.p.)", comb_folder_dict[:iota][2]],
                                         :xi => ["\\xi", comb_folder_dict[:xi][2]],
                                         :kappa => ["\\kappa \\, (b.p.)", comb_folder_dict[:kappa][2]],
                                         :lambda => ["\\lambda", comb_folder_dict[:lambda][2]],
                                         :sigmal => ["\\sigma_l", comb_folder_dict[:sigmal][2]],
                                         :sigmah => ["\\sigma_h", comb_folder_dict[:sigmah][2]])


zlabels = Dict{Symbol, Array{String,1}}(:c => ["Coupon", "%.2f"],
                                        :p => ["Principal", "%.2f"],
                                        :vb => ["VB", "%.1f"],
                                        :debt => ["Debt", "%.1f"],
                                        :equity => ["Equity", "%.1f"],
                                        :firm_value => ["Debt + Equity", "%1d"],
                                        :leverage => ["Leverage", "%1d"],
                                        :MBR => ["Market-to-Book Ratio", "%1d"])

svm_plots_title_params_order = [:mu_b, :m, :iota, :xi, :kappa, :lambda, :sigmal, :sigmah]

heat_surf_graph_format = "eps"
# ########################################################################


# * CVM v.s. SVM & Misrepresentation Plots ###############################
rmp_fname_ext = "png"
fixed_vars = [:mu_b, :m, :xi, :sigmal]
cvs_xvars = [:kappa, :lambda, :sigmah]


# Subplots #######################################
ax_subplots = Dict{Int64, Array{Int64,1}}(1 => [111],
                                          2 => [211, 212])


# Axes ###########################################
cvs_xlabels = Dict{Symbol, Array{String,1}}(:mu_s => ["\\mu_{s}", "%.2f"],
                                            :kappa_ep => ["\\kappa^{{EP}}", comb_folder_dict[:kappa][2]],
                                            :kappa_otc => ["\\kappa^{OTC}", comb_folder_dict[:kappa][2]],
                                            :iota => ["\\iota", comb_folder_dict[:iota][2]],
                                            :lambda => ["\\lambda", comb_folder_dict[:lambda][2]],
                                            :sigmah => ["\\sigma_h", comb_folder_dict[:sigmah][2]])

cvs_ylabels = Dict(zip([:firm_value, :equity, :debt, 
                        :c, :p, :vb, :leverage, :MBR],
                       ["Firm Value", "Equity", "Debt", "Coupon", 
                        "Principal", "\$ V^B\$", "Leverage", "Market-to-Book Ratio (\$\\%\$)"]))


# Markers and Line Styles ########################
cvmlinestyles = ["-", "-.", "--"]
cvmmarkers = ["d", "1", "2"]
svmlinestyles = ["-", "-.", "--"]
svmmarkers = ["", "d", "o"]


# Title ###########################################
cvs_plots_title_params_order = [:mu_b, :m, :xi, :sigmal]
tlabels = deepcopy(xylabels)


# Curves ##########################################
cvm_curve_label = "\\overline{\\iota} \\geqslant 0"
cvm_curve_color = "green"
svm_curve_color = "blue"
misrep_curve_color = "red"
# svm_curve_fv_label_ypos = 


# Vertical Lines ##################################
fv_color = "black"
mbr_color = "blueviolet"
misrep_color = "red"
                                              
function return_xsym_dict(xvar::Symbol)
    sub_sup = (xvar == :iota) ? "_" : "^"
    return Dict{Symbol, String}(:fv => string("\n \$", cvs_xlabels[xvar][1], sub_sup, "{fv}\$"),
                                :mbr => string("\n \$", cvs_xlabels[xvar][1], sub_sup, "{mbr}\$"),
                                :misrep => string("\n \$", cvs_xlabels[xvar][1],
                                                  sub_sup, "{mp}\$"),
                                :misrep1 => string("\n \$", cvs_xlabels[xvar][1],
                                                   sub_sup, "{mp, 1}\$"),
                                :misrep2 => string("\n \$", cvs_xlabels[xvar][1],
                                                   sub_sup, "{mp, 2}\$"))
end

function vlines_labels_dict(xvar; fv_xvar::Float64=NaN,
                            fv_color::String=fv_color,
                            mbr_xvar::Float64=NaN,
                            mbr_color::String=mbr_color,
                            cvm_misrep_xvar::Float64=NaN,
                            svm_misrep_xvar::Float64=NaN,
                            misrep_color::String=misrep_color)
    xsym_dict = return_xsym_dict(xvar)

    if any([isnan(cvm_misrep_xvar), isnan(svm_misrep_xvar)])
        misrep_xvar = isnan(cvm_misrep_xvar) ? svm_misrep_xvar : cvm_misrep_xvar
        
        return Dict(:firm_value => Dict(zip([:value, :xsym, :color],
                                          [fv_xvar, xsym_dict[:fv], fv_color])),
                    :MBR => Dict(zip([:value, :xsym, :color],
                                   [mbr_xvar, xsym_dict[:mbr], mbr_color])),
                    :misrep => Dict(zip([:value, :xsym, :color],
                                        [misrep_xvar, xsym_dict[:misrep], misrep_color])))
    else
        misrep1_xvar = minimum([cvm_misrep_xvar, svm_misrep_xvar])
        misrep2_xvar = maximum([cvm_misrep_xvar, svm_misrep_xvar])
        
        return Dict(:firm_value => Dict(zip([:value, :xsym, :color],
                                          [fv_xvar, xsym_dict[:fv], fv_color])),
                    :MBR => Dict(zip([:value, :xsym, :color],
                                   [mbr_xvar, xsym_dict[:mbr], mbr_color])),
                    :misrep1 => Dict(zip([:value, :xsym, :color],
                                         [misrep1_xvar, xsym_dict[:misrep1], misrep_color])),
                    :misrep2 => Dict(zip([:value, :xsym, :color],
                                        [misrep2_xvar, "", misrep_color])))
    end
end


# Region Colors ####################################
rm_region_color = "blue"
nrm_region_color = "#76D7C4"
conflict_region_color = "#EB984E"
misrep_region_color = "#F1948A"
box_color = "#EE0839"


# ##################################################
# tlabels = Dict(zip(vcat(fixed_vars, yvars), 
#                         ["\\mu_b", "m", "\\xi", "\\sigma_l", 
#                          "\\kappa^{EP}", "\\lambda", "\\sigma_h"]))
# title_params = join([string("\$", tlabels[x], "= \$ ", svmcombs[1][x])
#                      for x in vcat(fixed_vars, yvars) if 
#                      !(x in [:sigmah, :lambda, :kappa, :iota])], ", ")
# ########################################################################


# * Misrepresentation Plots ################################################
rmp_plots_title_params_order =  [:mu_b, :m, :xi, :kappa, :sigmal]
rmp_fn_prefix = "rmp"
rmp_full_info_prefix = "fi"
rmp_misrep_prefix = "misrep"
rmp_fig_aspect = .5
rmp_multi_plot_fig_size = (10., 8.)
rmp_multi_plot_figpad = .9

# ########################################################################

jeq_xlabels = Dict{Symbol, Array{String,1}}(:mu_s => ["\\mu_{s}", "%.2f"])

# Markers and Line Styles ########################
fi_linestyles = ["-", "-.", "--"]
fi_markers = ["d", "1", "2"]


# Line Styles
jeq_linestyles = ["-", "-.", "--"]


# Curve Colors
fi_ep_curve_color = "blue"
fi_otc_curve_color = "#DC7633"
# fi_otc_curve_color = "#A569BD"
pool_curve_color = "red"
sep_curve_color = "#17A589"


jeq_markers = ["", "d", "o"]
fi_curve_color = "blue"

jeq_plots_title_params_order =  [:m, :xi, :sigmal]
otc_region_color = "#F0B27A"
# ########################################################################


# * Contour Plots ##########################################################
contour_fname_ext = "eps"

resdict = Dict{Symbol,Any}(:xvar => :Symbol,
                           :yvar => :Symbol,
                           :zvar => :Symbol,
                           :x_fun => Spline1D,
                           :y_fun => Spline1D,
                           :xy_fun => Spline1D,
                           :df => DataFrame())


contour_tlabels = Dict{Symbol, Array{String,1}}(:mu_s => ["\\mu_s", "%.2f"],
                                                :lambda => ["\\lambda", "%.3f"],
                                                :m => ["m", "%.2f"],
                                                :kappa => ["\\kappa^{EP} \\, (b.p.)", "%.2f"],
                                                :kappa_otc => ["\\kappa^{OTC} \\, (b.p.)", "%.2f"],
                                                :xi => ["\\xi", "%.2f"],
                                                :sigmal => ["\\sigma_l", "%.3f"],
                                                :c => ["Coupon", "%.2f"],
                                                :p => ["Principal", "%.2f"],
                                                :vb => ["VB", "%.1f"],
                                                :debt => ["Debt", "%.1f"],
                                                :equity => ["Equity", "%.1f"],
                                                :firm_value => ["Firm Value", "%1d"],
                                                :leverage => ["Leverage", "%1d"],
                                                :MBR => ["Market-to-Book Ratio", "%1d"],
                                                :pcr => ["p/c", "%.1f"])


contour_diff = ["v.s. Full Information Eq. ", "differential"]


contour_xvar = :iota
contour_yvar = :sigmah
contour_zvars_sym=Dict{Symbol, Symbol}(:MBR => :mbr, :firm_value => :fv, :leverage => :lev)
contour_zvars=Array([x for x in keys(contour_zvars_sym)])

contour_firm_types = Dict{Symbol, Symbol}(:safe => :s_, :risky => :r_)

# contour_plots_title_params_order = [:m, :xi, :kappa, :lambda, :sigmal]
contour_plots_title_params_order = [:m, :pcr, :xi, :kappa, :lambda, :sigmal]

eq_type_title = Dict{String, Array{Any,1}}("full_info" => [:fi, "Full Information"],
                                           "misrep" => [:mp, "Misrepresentation"],
                                           "pooling" => [:pool, "Pooling"],
                                           "separating" => [:sep, "Separating"])


iso_cmaps = Dict{String, Any}("full_info" => Seaborn.get_cmap("YlGnBu_r"),
                              "misrep" => Seaborn.palplot("Reds"),
                              "pooling" => "BuPu",
                              "separating" => "RdPu")


iso_plt_inputs = Dict{Symbol,Any}(:seaborn_style => "white", 
                                  :iso_levels => 20, 
                                  :heat_levels => 25, 
                                  :fig_aspect => .4, 
                                  :iso_fontsize => 9.,
                                  :use_subgrid => true,
                                  :subgrid_rows => 3, 
                                  :iso_cols => 6, 
                                  :heat_cols => 4,
                                  :title_font_size => 14.5,
                                  :fig_dpi => 300,
                                  :tight_pad => 3., 
                                  :h_pad => .75,
                                  :w_pad => .75)


eq_cat_dict = Dict{Symbol, Array{Any,1}}(:fi => [4, "FI"],
                                         :sep => [3, "SEP"],
                                         :pool => [2, "POOL"],
                                         :otc => [1, "OTC"])

# obj_fun_symbol = Dict{Symbol, Symbol}(:MBR => :mbr, :firm_value => :fv)
# ########################################################################
