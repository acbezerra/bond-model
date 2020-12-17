# vim: set fdm=marker :

dir_vec = rsplit(pwd(), "/")
pos = findall(x -> x .== "bond-model", dir_vec)[1]
main_path = join(dir_vec[1:pos], "/")
module_path=string(main_path, "/modules")
push!(LOAD_PATH, module_path)
# modnames = ["Batch"]
# for modl in modnames
#     if !(joinpath(module_path, modl) in LOAD_PATH)
#         push!(LOAD_PATH, joinpath(module_path, modl))
#     end
# end


module ModelPlots

using Distributed
using Dierckx

using Parameters
using Printf
using DataFrames
using CSV

# Plots
using PyPlot
using Seaborn
using LaTeXStrings

# User Defined Modules
using Batch: BatchStruct,
             main_dir,
             mat_dir_prefix,
             comb_folder_dict,
             form_main_dir_path,
             load_cvm_opt_results_df,
             load_svm_opt_results_df,
             str_format_fun, get_bt,
             form_main_dir_path, main_dir,
             BatchObj, get_bt,
             get_batch_comb_numbers,
             svm_param_values_dict

# * Plot Inputs {{{1
# include("_plots_inputs.jl")

# ** Directories & File Names {{{2
main_dir_path = form_main_dir_path(main_dir)
plots_dir="Plots"
contour_plots_path = string(main_dir_path, "/Plots/Contour")

cvm_vs_svm_plots_dir = "CVMvsSVM"
rmp_plots_dir = "RMP"
nrmp_plots_dir = "NRMP"
heat_surf_graph_dir="HeatSurf"
obj_fun_dict = Dict{Symbol, String}(:firm_value => "fv",
                                    :MBR => "mbr")
# ########################################################################


# ** Auxiliary Functions {{{2
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
                                 :MBR => "MBR (\$\\%\$)",
                                 :yield_spd => "Bond Spreads (b.p.)")



cvm_plots_title_params_order = [:mu_b, :m, :iota, :xi, :kappa, :sigma]




# ** SVM HeatMap and Surface Plots {{{2
xylabels = Dict{Symbol, Array{String,1}}(:mu_b => ["\\mu_b", comb_folder_dict[:mu_b][2]],
                                         :m => ["m", comb_folder_dict[:m][2]],
                                         :iota => ["\\iota \\, (b.p.)", comb_folder_dict[:iota][2]],
                                         :xi => ["\\xi", comb_folder_dict[:xi][2]],
                                         :kappa => ["\\kappa \\, (b.p.)", comb_folder_dict[:kappa][2]],
                                         :lambda => ["\\lambda", comb_folder_dict[:lambda][2]],
                                         :sigmal => ["\\sigma_l", comb_folder_dict[:sigmal][2]],
                                         :sigmah => ["\\sigma_h", comb_folder_dict[:sigmah][2]],
                                         :mu_s => ["\\mu_s", "%.2f"])


zlabels = Dict{Symbol, Array{String,1}}(:c => ["Coupon", "%.2f"],
                                        :p => ["Principal", "%.2f"],
                                        :vb => ["VB", "%.1f"],
                                        :debt => ["Debt", "%.1f"],
                                        :equity => ["Equity", "%.1f"],
                                        :firm_value => ["Debt + Equity", "%1d"],
                                        :fv => ["Debt + Equity", "%1d"],
                                        :leverage => ["Leverage", "%1d"],
                                        :MBR => ["Market-to-Book Ratio", "%1d"],
                                        :mbr => ["Market-to-Book Ratio", "%1d"],
                                        :yield_spd => ["Bond Spreads (b.p.)", "%1d"])

svm_plots_title_params_order = [:mu_b, :m, :iota, :xi, :kappa, :lambda, :sigmal, :sigmah]

heat_surf_graph_format = "eps"
# ########################################################################


# ** CVM v.s. SVM & Misrepresentation Plots {{{2
rmp_fname_ext = "png"
fixed_vars = [:mu_b, :m, :xi, :sigmal]
cvs_xvars = [:kappa, :lambda, :sigmah]

# Subplots
ax_subplots = Dict{Int64, Array{Int64,1}}(1 => [111],
                                          2 => [211, 212])


# Axes
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

        if .&(isnan(misrep_xvar), abs(fv_xvar - mbr_xvar) < 1e-5)

            xvar_label = xvar == :sigmah ? latexstring("\$ \\sigma_l \$") : ""
            return Dict(:misrep => Dict(zip([:value, :xsym, :color],
                                            [misrep_xvar, xvar_label, misrep_color])))
        elseif .&(abs(fv_xvar - misrep_xvar) < 1e-5, abs(mbr_xvar - misrep_xvar) < 1e-5)
            # if isnan(cvm_misrep_xvar)
            #     xvar_label = xvar == :sigmah ? latexstring("\$ \\sigma_l \$") : ""
            # else
            #     xvar_label = latexstring(xsym_dict[:fv],  " = ",
            #                              xsym_dict[:mbr], " = ", xsym_dict[:misrep])
            # end
            xvar_label = xvar == :sigmah ? latexstring("\$ \\sigma_l \$") : ""

            return Dict(:misrep => Dict(zip([:value, :xsym, :color],
                                             [misrep_xvar, xvar_label, misrep_color])))
       else
            return Dict(:firm_value => Dict(zip([:value, :xsym, :color],
                                                [fv_xvar, xsym_dict[:fv], fv_color])),
                        :MBR => Dict(zip([:value, :xsym, :color],
                                         [mbr_xvar, xsym_dict[:mbr], mbr_color])),
                        :misrep => Dict(zip([:value, :xsym, :color],
                                            [misrep_xvar, xsym_dict[:misrep], misrep_color])))
        end
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
                                        [misrep2_xvar, xsym_dict[:misrep2], misrep_color])))
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


# ** Misrepresentation Plots {{{2
rmp_plots_title_params_order =  [:mu_b, :m, :xi, :kappa, :sigmal]
rmp_fn_prefix = "rmp"
nrmp_fn_prefix = "nrmp"
full_info_prefix = "fi"
misrep_prefix = "misrep"
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


# ** Contour Plots {{{2
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
                                                :fv => ["Firm Value", "%1d"],
                                                :leverage => ["Leverage", "%1d"],
                                                :mbr => ["Market-to-Book Ratio", "%1d"],
                                                :yield => ["Bond Yield", "%.2f"],
                                                :yield_spd => ["Bond Spread (b.p.)", "%.f"],
                                                :pcr => ["p/c", "%.1f"])


contour_diff = ["v.s. Full Information Eq. ", "differential"]


contour_xvar = :iota
contour_yvar = :sigmah
contour_zvars_sym=Dict{Symbol, Symbol}(:MBR => :mbr, :firm_value => :fv, :leverage => :lev,
                                      :yield => :yield, :yield_spd => :spd)
contour_zvars=Array([x for x in keys(contour_zvars_sym)])

contour_firm_types = Dict{Symbol, Symbol}(:safe => :s_, :risky => :r_)

# contour_plots_title_params_order = [:m, :xi, :kappa, :lambda, :sigmal]
contour_plots_title_params_order = [:m, :pcr, :xi, :kappa, :lambda, :sigmal]

eq_type_title = Dict{String, Array{Any,1}}("full_info" => [:fi, "Full Information"],
                                           "misrep" => [:mp, "Misrepresentation"],
                                           "pooling" => [:pool, "Pooling"],
                                           "separating" => [:sep, "Separating"],
                                           "ep_market" => [:ep_market, "Prevailing EP Market Equilibria"],
                                           "dual_market" => [:dual_market, "Prevailing Dual Market Equilibria"])


iso_cmaps = Dict{String, Any}("full_info" => Seaborn.get_cmap("YlGnBu_r"),
                              "misrep" => Seaborn.palplot("Reds"),
                              "pooling" => "BuPu",
                              "separating" => "RdPu")

rmp_cmaps = Dict{String, Any}("misrep" => "GnBu",
                              "pooling" => "BuGn")

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


# * Plot Constructors {{{1
@with_kw mutable struct PlotStruct
    # Batch Obj
    bt

    # Results DataFrame
    cvm_data::DataFrame
    svm_data::DataFrame

    # Firm Objective Function
    obj_fun::Symbol

    # Surface Slice
    _svm_surf::DataFrame

    # Labels
    xylabels::Dict{Symbol, Array{String,1}}
    zlabels::Dict{Symbol,Array{String,1}}
end


function PlotsObj(bt;
                  load_cvm_data::Bool=true,
                  load_svm_data::Bool=true,
                  firm_obj_fun::Symbol=:firm_value, #:MBR,
                  cvm_m::Float64=NaN,
                  svm_m::Float64=NaN,
                  xylabels::Dict{Symbol, Array{String, 1}}=xylabels,
                  zlabels::Dict{Symbol,Array{String,1}}=zlabels)

    cvm_data = DataFrame()
    if load_cvm_data
        try
            cvm_data = load_cvm_opt_results_df(; m=cvm_m,
                                               firm_obj_fun=firm_obj_fun)
        catch
            println("Unable to load CVM data")
        end
    end

    svm_data = DataFrame()
    if load_svm_data
        try
            svm_data = load_svm_opt_results_df(bt;
                                               firm_obj_fun=firm_obj_fun,
                                               m=svm_m)
        catch
            println("Unable to load SVM data")
        end
    end

    return PlotStruct(bt,
                      cvm_data,
                      svm_data,
                      firm_obj_fun,
                      DataFrame(),
                      xylabels, zlabels)
end


# * Auxiliary Methods {{{1
# include("_plots_auxiliary_methods.jl")
function set_cvm_data(pt::PlotStruct, df::DataFrame)
    if isempty(pt.cvm_data)
        println("Setting the Data")
    else
        println("Updating the Data")
    end

    pt.cvm_data=df
end


function plots_form_combinations(bt, fig_name_vars::Array{Symbol,1})
    value_lists = [bt.bp._param_values_dict[x] for x in fig_name_vars]

    # Remaining Params
    # rparams = [x for x in bt._params_order if !(x in fig_name_vars)]

    combs  = [[x1, x2, x3, x4, x5, x6] for
              x1 in bt.bp._param_values_dict[fig_name_vars[1]],
              x2 in bt.bp._param_values_dict[fig_name_vars[2]],
              x3 in bt.bp._param_values_dict[fig_name_vars[3]],
              x4 in bt.bp._param_values_dict[fig_name_vars[4]],
              x5 in bt.bp._param_values_dict[fig_name_vars[5]],
              x6 in bt.bp._param_values_dict[fig_name_vars[6]]]

    # combinations  = hcat(combs...)'

    return combs
end


function plot_svmdf_slicer(pt, fixed_params::Dict{Symbol, Float64})
    locs = [abs.(pt.svm_data[:, x] .- fixed_params[x]) .< 1e-6 for x in keys(fixed_params)]
    srows = sum(hcat(locs...), dims=2) .== size(hcat(locs...), 2)
    return pt.svm_data[vcat(srows...), :]
end


function get_cvm_svm_dfs(cvmdict::Dict{Symbol,Array{Float64,1}},
                         svmdict::Dict{Symbol,Array{Float64,1}};
                         firm_obj_fun::Symbol=:firm_value)
    m = (:m in keys(cvmdict)) ? cvmdict[:m][1] : NaN

    bt = BatchObj()
    pt = PlotsObj(bt; firm_obj_fun=firm_obj_fun, cvm_m=m, svm_m=m)
    sbt = get_bt(; model="svm", m=m, m_comb_num=1)
    cbt = get_bt(; model="cvm", m=m, m_comb_num=1)

    # Set Batch Objects
    cbt = get_bt(;model="cvm", m=m, m_comb_num=1)
    sbt = get_bt(;model="svm", m=m, m_comb_num=1)

    # Get Combination Numbers
    cvm_combs = get_batch_comb_numbers(cbt, cvmdict)[:, :comb_num]
    svm_combs = get_batch_comb_numbers(sbt, svmdict)[:, :comb_num]
    # #########################################################


    # Safe and Risky Firms' DFs ###############################
    cvmdf = pt.cvm_data[[(x in cvm_combs) for x in pt.cvm_data[:, :comb_num]], :]
    svmdf = pt.svm_data[[(x in svm_combs) for x in pt.svm_data[:, :comb_num]], :]
    # #########################################################

    # x-axis variable:
    xvar =  [x for x in cvs_xvars if size(unique(svmdf[:, x]), 1) == size(unique(pt.svm_data[:, x]), 1)]
    if !isempty(xvar)
        if size(xvar, 1) > 1
            println("Multiple xvars!")
        end
        xvar = xvar[1]
    end

    return cvmdf, svmdf, xvar
end


# * CVM Plot Methods {{{1
# include("_cvm_plot_methods.jl")
function cvm_data_handler(pt)
    kappa_vals = unique(pt.cvm_data[:, :kappa])
    sigmal_vals = unique(pt.cvm_data[:, :sigmal])
    iota_vals = [ i*10^4 for i in unique(pt.cvm_data[:, :iota]) if i < maximum(pt.cvm_data[:, :iota])]
    # scalarMap = Seaborn.color_palette("cool", n_colors=size(sigma_vals, 1))
    scalarMap = Seaborn.get_cmap("cool", size(sigmal_vals, 1))
    l_styles = ["-" , ":", "--", "-."]

    return kappa_vals, sigmal_vals, iota_vals, scalarMap, l_styles
end


function cvm_single_axis_plot(pt, fig, plot_sub_num::Int64, varname::Symbol;
                              kappa::Float64=NaN,
                              sigma_low::Float64=NaN,
                              legend::Bool=true)
    (kappa_vals, sigmal_vals,
     iota_vals, scalarMap, l_styles) = cvm_data_handler(pt)

    if isnan(kappa)
        kappa = kappa_vals[1]
    end
    if isnan(sigma_low)
        sigma_low = minimum(sigmal_vals)
    end

    ax = fig.add_subplot(plot_sub_num)
    for j in 1:size(sigmal_vals, 1)
        colorVal = scalarMap(j)
        pos = .&(abs.(pt.cvm_data[:, :kappa] .- kappa) .< 1e-6,
                 abs.(pt.cvm_data[:, :sigmal] .- sigmal_vals[j]) .< 1e-6)

        if (abs.(sigmal_vals[j] .- sigma_low) .< 1e-6)
            pos = .&(pos, [(x .* 10^4 in iota_vals) for x in pt.cvm_data[:, :iota]])

            ax.plot(pt.cvm_data[pos, :iota] .* 10^4,
                    pt.cvm_data[pos, varname],
                    linewidth=1.1,
                    linestyle=l_styles[j],
                    color=colorVal,
                    label=sigmal_vals[j])
        else
            firm_val = pt.cvm_data[pos, varname][1]
            ax.axhline(y=firm_val,
                       linewidth=1.1,
                       linestyle=l_styles[j],
                       color=colorVal,
                       label=sigmal_vals[j],
                       xmin=.01, xmax=.99)
         end
    end

    if legend
        ax.legend(loc=0,
                   title="Volatility \$\\sigma\$",
                   ncol=1,
                   frameon=true,
                   shadow=true,
                   framealpha=.85,
                   edgecolor="white",
                   facecolor="white",
                   fancybox=true)
        # bbox_to_anchor=(.95, 0.85),

        ax.get_legend().get_title().set_color("#34495E")
    end
    # Axes' Labels:
    # if ylabel:
    #     ax.set_ylabel('Values for $\sigma = \overline{\sigma}$ and $\iota\geqslant 0$',
    #                    fontsize=12, labelpad=10)
    ax.set_xlabel("Risk Management Cost, \$\\iota\$ (b.p.)",
                   fontsize=12, labelpad=10)
    return ax
end


function cvm_double_axes_plot(pt, fig,
                              plot_sub_num::Int64, varname::Symbol;
                              kappa::Float64=NaN,
                              sigma_low::Float64=NaN,
                              y1label::Bool=true, y2label::Bool=true,
                              legend::Bool=true)
    (kappa_vals, sigmal_vals,
     iota_vals, scalarMap, l_styles) = cvm_data_handler(pt)

    if isnan(kappa)
        kappa = kappa_vals[1]
    end
    if isnan(sigma_low)
        sigma_low = minimum(sigmal_vals)
    end

    ax1 = fig.add_subplot(plot_sub_num)
    ax2 = ax1.twinx()
    for j in 1:size(sigmal_vals, 1)
        colorVal = scalarMap(j)
        pos = .&(abs.(pt.cvm_data[:, :kappa] .- kappa) .< 1e-6,
                 abs.(pt.cvm_data[:, :sigmal] .- sigmal_vals[j]) .< 1e-6)

        if (abs.(sigmal_vals[j] .- sigma_low) .< 1e-6)
            pos = .&(pos, [(x .* 10^4 in iota_vals) for x in pt.cvm_data[:, :iota]])


            ax1.plot(pt.cvm_data[pos, :iota] .* 10^4,
                     pt.cvm_data[pos, varname],
                     linewidth=1.1,
                     linestyle=l_styles[j],
                     color=colorVal,
                     label=sigmal_vals[j])
        else
            firm_val = pt.cvm_data[pos, varname][1]
            ax2.axhline(y=firm_val,
                        linewidth=1.1,
                        linestyle=l_styles[j],
                        color=colorVal,
                        label=sigmal_vals[j],
                        xmin=.01, xmax=.99)
        end
    end


    if legend
        # handles,labels = ax1.get_legend_handles_labels()
        handles, labels = [ ], [ ]
        # for ax in fig.axes:
        for ax in [ax1, ax2]
            for (h ,l) in zip(ax.get_legend_handles_labels()...)
                push!(handles, h)
                push!(labels, l)
            end
        end

        ax1.legend(handles, labels,
                   loc=0,
                   title="Volatility \$\\sigma\$",
                   ncol=1,
                   frameon=true,
                   shadow=true,
                   framealpha=.85,
                   edgecolor="white",
                   facecolor="white",
                   fancybox=true)
                   # bbox_to_anchor=(.95, 0.85),

        ax1.get_legend().get_title().set_color("#34495E")
    end

    # Axes' Labels:
    if y1label
        ax1.set_ylabel(string("Values for \$\\sigma = \\underline{\\sigma}\$",
                              " and \$\\iota \\geq 0\$"), usetex=true,
                       fontsize=12, labelpad=10)
    end

    if y2label
        ax2.set_ylabel(string("Values for \$\\sigma > \\underline{\\sigma}\$",
                              " and \$\\iota = 0\$"), usetex=true,
                       fontsize=12, labelpad=10)
    end

    ax1.set_xlabel("Risk Management Cost, \$\\iota\$ (b.p.)",
                   fontsize=12, labelpad=10)

    ax2.grid(nothing)

    return ax1

end


function cvm_plot_path_fname(pt, var::Symbol, fixed_params::Dict{Symbol, Any};
                             title_params_order::Array{Symbol,1}=cvm_plots_title_params_order,
                             main_dir::String=main_dir,
                             plots_dir::String=plots_dir,)

    # Paths to Plot ############################################
    # Main Directory
    main_dir_path = form_main_dir_path(main_dir)

    comb_dir_params = [x for x in title_params_order if x in keys(fixed_params)]
    comb_dir = join([str_format_fun(par_val_printer(x),
                                    par_val_adj(x, fixed_params[x]))
                     for x in comb_dir_params], "__")

    # Form Paths ###############################################
    dirs = [plots_dir, "CVM", comb_dir]
    graphs_path = main_dir_path
    for dir in dirs
        graphs_path = string(graphs_path, "/", dir)
        if !isdir(graphs_path)
            mkdir(graphs_path)
        end
    end
    # ##########################################################

    # Filename #################################################

    file_name = string("cvm_", obj_fun_dict[fixed_params[:obj_fun]],
                       "_opt_", var)

    return [graphs_path, file_name]
end


function plot_cvm_optimal_solutions(pt, var::Symbol;
                                    kappa::Float64=NaN,
                                    sigma_low::Float64=NaN,
                                    title_params_order::Array{Symbol,1}=cvm_plots_title_params_order,
                                    figaspect::Float64=.55,
                                    facecolor::String="w",
                                    save_fig::Bool=true,
                                    return_fig::Bool=true,
                                    fig_dpi::Int64=400,
                                    graph_format::String="png")

    if isnan(kappa)
        kappa = minimum(unique(pt.cvm_data[:, :kappa]))
    end
    if isnan(sigma_low)
        sigma_low = minimum(unique(pt.cvm_data[:, :sigmal]))
    end

    Seaborn.set_style("darkgrid")
    fig = PyPlot.figure(figsize=PyPlot.figaspect(figaspect), facecolor=facecolor)
    if !(var in [:firm_value, :MBR])
        ax1 = ModelPlots.cvm_double_axes_plot(pt, fig, 111, var;
                                              kappa=kappa, sigma_low=sigma_low)
    else
        ax1 = ModelPlots.cvm_single_axis_plot(pt, fig, 111, var;
                                              kappa=kappa, sigma_low=sigma_low)
    end

    ax1.set_title(string("Constant Volatility Model Optimal ", vartitles[var]), fontsize=14)
    PyPlot.tight_layout()


    if save_fig
        fixed_params = Dict{Symbol, Any}(:obj_fun => Symbol(pt.cvm_data[1, :obj_fun]),
                                         :mu_b => unique(pt.cvm_data[:, :mu_b])[1],
                                         :m => unique(pt.cvm_data[:, :m])[1],
                                         :xi => unique(pt.cvm_data[:, :xi])[1],
                                         :kappa => kappa,
                                         :sigmal => sigma_low)

        folder_path, file_name = cvm_plot_path_fname(pt, var, fixed_params;
                                                     title_params_order=title_params_order)
        PyPlot.savefig(string(folder_path, "/", file_name, ".", graph_format), dpi=fig_dpi, format=graph_format)
    end

    if return_fig
        return fig
    end
end


# * SVM Plot Methods {{{1
# include("_svm_plot_methods.jl")

# ** Slice and Interpolate Surface {{{2
function set_svm_surf(pt, fixed_params::Dict{Symbol, Float64})
    if isempty(pt._svm_surf)
        println("Setting SVM Slice Data")
    else
        println("Updating SVM Slice Data")
    end

    pt._svm_surf = plot_svmdf_slicer(pt, fixed_params)

    return pt
end


function svm_interp_smooth_surface(pt, fixed_params::Dict{Symbol,Float64},
                                   z_var::Symbol;
                                   interp::Bool=true,
                                   smooth::Bool=true)

    # Variables
    xy_vars = [:kappa, :lambda, :sigmah]

    kappa_vec = unique(pt._svm_surf[:, :kappa])
    sigmah_vec = unique(pt._svm_surf[:, :sigmah])
    lambda_vec = unique(pt._svm_surf[:, :lambda])

    # xy axes - pick 2 variables from (kappa, lambda, sigmah):
    xy_list = [var for var in xy_vars if
                !(var in keys(fixed_params))]

    # Sigmah on Y-axis, Kappa on X-axis:
    # if Sigmah & Lambda -> lambda on X-Axis
    if .&((:kappa in xy_list), (:sigmah in xy_list))
        x_vec = kappa_vec
        y_vec = sigmah_vec
        var_axes = [:kappa, :sigmah]
    elseif .&((:kappa in xy_list), (:lambda in xy_list))
        x_vec = kappa_vec
        y_vec = lambda_vec
        var_axes = [:kappa, :lambda]
    else
        x_vec = lambda_vec
        y_vec = sigmah_vec
        var_axes = [:lambda, :sigmah]
    end

    # Form Combinations:
    comb = hcat([Array([x, y]) for x in x_vec, y in  y_vec]...)'

    # #####################################################################
    # ####################### Extract z-var values: #######################
    # #####################################################################
    z_val = Array{Float64}(undef, size(x_vec, 1) .* size(y_vec, 1), 1)

    # Boolean Indexer for Fixed Parameter Values:
    pos_loc = [(abs.(pt._svm_surf[:, x] .- fixed_params[x]) .< 1e-6)
                                for x in keys(fixed_params)]
    bool_loc =  .&(pos_loc...) # sum(pos_loc) .== size(pos_loc, 1)

    for i in 1:size(comb, 1)
        if .&((:kappa in xy_list), (:sigmah in xy_list))
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:, :kappa] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:, :sigmah] .- comb[i, 2]) .< 1e-6),
                                    bool_loc), z_var][1]
        elseif .&((:kappa in xy_list), (:lambda in xy_list))
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:, :kappa] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:, :lambda] .- comb[i, 2]) .< 1e-6),
                                     bool_loc), z_var][1]
        else
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:, :lambda] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:, :sigmah] .- comb[i, 2]) .< 1e-6),
                                    bool_loc), z_var][1]
        end

        # Keep in mind that not all combinations of x_vec and
        # y_vec values may be available.
        if isempty(value)
            z_val[i] = NaN
        else
            z_val[i] = value
        end
    end

    # #####################################################################
    # ########################### Form Surface ############################
    # #####################################################################
    # Form Surface
    xyz = hcat(comb, z_val)

    # Remove combinations for which output is NaN:
    xyz[vec(sum(isnan.(xyz), dims = 2) .== 0), :]

    # Interpolate Surface
    if interp
        # Spline: Setting s to length xyz gives a SmoothBivariateSpline
        xyz_interp = Dierckx.Spline2D(xyz[:, 1], xyz[:, 2], xyz[:, 3];
                                      kx=3, ky=3,s=size(xyz, 1))

        # Increase Data Points

        # To spot inconsistencies and errors in the
        # manipulation of grid points, I use
        # rectangular xy grid:
        x_step_num = 100
        y_step_num = 105
        x_grid_ref = range(minimum(xyz[:, 1]), stop=maximum(xyz[:, 1]), length=x_step_num)
        y_grid_ref = range(minimum(xyz[:, 2]), stop=maximum(xyz[:, 2]), length=y_step_num)
        xy_grid_ref = hcat([Array([x, y]) for x in x_grid_ref, y in y_grid_ref]...)'
        grid_x = xy_grid_ref[:, 1]
        grid_y = xy_grid_ref[:, 2]

        # Back-out Interpolated Values:
        grid_z = evaluate(xyz_interp, grid_x, grid_y)
    else
        grid_x = xyz[:, 1]
        grid_y = xyz[:, 2]
        grid_z = xyz[:, 3]
    end

    # Kernel Surface Smoother
#     f = None
#     if smooth:
#         f = Rbf(grid_x, grid_y, grid_z, smooth=0.15, epsilon=1.)

    return [xyz, grid_x, grid_y, grid_z, var_axes]
end

# ** Heatmap Plot {{{2
function svm_plot_heatmap(pt, xy_list, z_var::Symbol,
                          fixed_params::Dict{Symbol,Float64};
                          plt_cmap::String="viridis",
                          order_ascend::Bool=true,
                          reverse_x_axis::Bool=false,
                          reverse_y_axis::Bool=false,
                          make_title::Bool=true,
                          add_title::Bool=false)
    # ############################################
    # ############### Axis Format: ###############
    # ############################################
    title_list = [string("\$", pt.xylabels[x][1], "=\$ ", pt.xylabels[x][2])
                        for x in keys(pt.xylabels) if !(x in vcat(xy_list, :mu_s))]
    title_values = [fixed_params[x] for x in keys(pt.xylabels) if !(x in vcat(xy_list, :mu_s))]


    svmdf = deepcopy(pt._svm_surf)

    # Localize entries:
    # Localize before changing kappa to basis points:
    cond = true
    for key in keys(fixed_params)
        cond = .&(cond, (abs.(svmdf[:, key] .- fixed_params[key]) .< 1e-6))
    end

    # Sigmah on Y-axis, Kappa on X-axis:
    # if Sigmah & Lambda -> lambda on X-Axis
    # Adjust value of Kappa
    if .&((:sigmah in xy_list), (:kappa in xy_list))
        cols = vcat([:sigmah, :kappa], [z_var])
        # Kappa in Basis Points:
        svmdf[!, :kappa] .= convert.(Int64, svmdf[:, :kappa] .* 1e4)
    elseif .&((:sigmah in xy_list), (:lambda in xy_list))
        cols = vcat([:sigmah, :lambda], [z_var])
        # Kappa in Basis Points
        # fixed_params[:kappa] = #convert.(Int64, fixed_params[:kappa] .* 1e4)
    else
        cols = vcat([:lambda, :kappa], [z_var])
        # Kappa in Basis Points:
        svmdf[!, :kappa] .=  convert.(Int64, svmdf[:, :kappa] .* 1e4)
    end

    # Pivot Table for Plot: ##############################
    xtickvals = unique(svmdf[:, cols[2]])
    ytickvals = unique(svmdf[:, cols[1]])
    pivotdf = unstack(svmdf[cond, cols], cols[1], cols[2], cols[3])
    Z = convert(Matrix, pivotdf[:, [Symbol(x) for x in names(pivotdf) if !(Symbol(x) in cols)]])

    if reverse_x_axis
        Z = reverse(Z, dims=2)
        xtickvals = reverse(xtickvals)
    end

    if reverse_y_axis
        Z = reverse(Z, dims=1)
        ytickvals = reverse(ytickvals)
    end

    ax = Seaborn.heatmap(Z, cmap=plt_cmap,
                         xticklabels=xtickvals,
                         yticklabels=ytickvals)
    ax.set_xlabel(string("\$", pt.xylabels[cols[2]][1], "\$"))
    ax.set_ylabel(string("\$", pt.xylabels[cols[1]][1], "\$"))
    # ##################################################


    # ############### TITLE ###############
    ax_title = " "
    if make_title | add_title
        # str_format_fun(a, b) = @eval @sprintf($a, $b)
        formatted_var_string = join([str_format_fun(title_list[i],
                                                    title_values[i])
                                     for i in 1:size(title_list, 1)], ", ")
        ax_title = string(pt.zlabels[z_var][1], " values for ",
                          formatted_var_string)
    end

    if add_title
        ax.set_title(ax_title)
    end

    return [ax, ax_title]

end

# ** Surface Plot {{{2
function svm_plot_surface(pt, ax, xv, yv, zv,
                          xy_list::Array{Symbol,1},
                          z_var::Symbol,
                          fixed_params::Dict{Symbol,Float64};
                          title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                          plt_cmap::String="viridis",
                          seaborn_style::String="darkgrid",
                          make_title=true, add_title=false,
                          zpad=10, view_elev=25., view_azim=210)
                          #xylabels::Array{Symbol,1}=[:mu_b, :m, :xi, :kappa, :lambda, :sigmal, :sigmah])
    # Documentation available at: https://matplotlib.org/mpl_toolkits/mplot3d/api.html

    # ATTENTION! -> only adjust kappa after computing zv:
    # Set Kappa to Basis Points:
    if xy_list[2] == :kappa
        xv = convert.(Int64, xv .* 10^4)
    elseif xy_list[1] == :kappa
        yv = convert.(Int64, yv .* 10^4)
    end

    cp_fixed_params = deepcopy(fixed_params)
    if :kappa in keys(cp_fixed_params)
        cp_fixed_params[:kappa] = convert.(Int64, cp_fixed_params[:kappa] .* 10^4)
    end
    # ############################################
    # ############### Axis Format: ###############
    # ############################################
    major_xticks = range(minimum(xv), stop=maximum(xv), length=6)
    major_yticks = range(minimum(yv), stop=maximum(yv), length=6)
    major_zticks = range(ceil(minimum(zv)), stop=ceil(maximum(zv)), length=10)
    if z_var == :c
        major_zticks = range(minimum(zv), stop=maximum(zv), length=7)
    end

    xy_format = [[string("\$",  pt.xylabels[x][1], "\$"), pt.xylabels[x][2]]
                    for x in keys(pt.xylabels) if (x in xy_list)]

    title_list = [string("\$", pt.xylabels[x][1], "=\$ ", pt.xylabels[x][2])
                        for x in title_params_order if !(x in xy_list)]
    title_values = [cp_fixed_params[x] for x in title_params_order if !(x in xy_list)]

    if isempty(seaborn_style)
        Seaborn.reset_orig()
    else
        Seaborn.set(style="darkgrid")
    end
    surf = ax.plot_trisurf(xv, yv, zv,
                           linewidth=0.2,
                           cmap=plt_cmap,
                           antialiased=true)

    # Customize the x axis
    # ax.xaxis_inverted()
    ax.invert_xaxis()
    ax.xaxis.set_rotate_label(false)  # disable automatic rotation
    ax.set_xlabel(xy_format[1][1], labelpad=10, rotation=0)
    ax.xaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(xy_format[1][2]))
    ax.set_xticks(major_xticks)

    # Customize the y axis
    ax.yaxis.set_rotate_label(false)  # disable automatic rotation
    ax.set_ylabel(xy_format[2][1], labelpad=10, rotation=0)
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(xy_format[2][2]))
    ax.set_yticks(major_yticks)

    # Customize the z axis
    ax.zaxis.set_rotate_label(false)  # disable automatic rotation
    ax.zaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter("%.2f"))
    ax.set_zlabel(pt.zlabels[z_var][1], rotation=90, labelpad=zpad)
    ax.set_zticks(major_zticks)
    ax.tick_params(axis="z", pad=zpad/2)

    ax.minorticks_on()
    ax.grid(which="minor", alpha=0.3)
    # ############### TITLE ###############
    ax_title = " "
    if make_title | add_title
        # str_format_fun(a, b) = @eval @sprintf($a, $b)
        formatted_var_string = join([str_format_fun(title_list[i],
                                                    title_values[i])
                                     for i in 1:size(title_list, 1)], ", ")
        ax_title = latexstring(pt.zlabels[z_var][1], " values for ",
                          formatted_var_string)
    end

    if add_title
        ax.set_title(ax_title)
    end

#     # Add a color bar which maps values to colors.
#     # fig.colorbar(surf, shrink=0.6, aspect=6)

    ax.view_init(view_elev, view_azim)

    return surf, ax, ax_title
end

# ** Heat Surf Auxiliary Functions {{{2
function svm_heat_surf_plot_path_fname(pt, xy_list::Array{Symbol,1}, z_var::Symbol,
                                       fixed_params::Dict{Symbol,Float64};
                                       title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                                       main_dir::String=main_dir,
                                       plots_dir::String=plots_dir,
                                       graph_dir::String=heat_surf_graph_dir,
                                       mat_dir_prefix::String=mat_dir_prefix)

    # Functions to Format Numeric-to-String
    # str_format_fun(a::String,b::Float64) = @eval @sprintf($a, $b)
    # par_val_printer(x::Symbol) = !(x in [:iota, :kappa]) ? string(x, "_", pt.xylabels[x][2]) : string(x, "_(bp)_", pt.xylabels[x][2])
    # par_val_adj(x::Symbol, val::Float64) = !(x in [:iota, :kappa]) ? val : val * 1e4

    # Paths to Plot ##########################################

    # Main Directory
    main_dir_path = form_main_dir_path(main_dir)

    # Maturity Directory
    mat_dir = string(mat_dir_prefix, str_format_fun(pt.xylabels[:m][2], fixed_params[:m]))

    # Get Graph Type
    graph_type = ""
    graph_sub_folder = ""
    if .&((:kappa in xy_list), (:sigmah in xy_list))
        graph_type = "kappa_sigmah"
        graph_sub_folder = str_format_fun(string("lambda_", pt.xylabels[:lambda][2]),
                                          fixed_params[:lambda])
    elseif .&((:kappa in xy_list), (:lambda in xy_list))
        graph_type = "kappa_lambda"
        graph_sub_folder = str_format_fun(string("sigmah_", pt.xylabels[:sigmah][2]),
                                          fixed_params[:sigmah])
    else
        graph_type = "lambda_sigmah"
        graph_sub_folder = str_format_fun(string("kappa_bp_", pt.xylabels[:kappa][2]),
                                          fixed_params[:kappa] * 10^4)
    end

    # Form Paths ##############################################
    dirs = [plots_dir, "SVM", graph_dir, mat_dir,
            graph_type, graph_sub_folder]
    graphs_path = main_dir_path
    for dir in dirs
        graphs_path = string(graphs_path, "/", dir)
        if !isdir(graphs_path)
            mkdir(graphs_path)
        end
    end
    # ##########################################################

    # Filename #################################################
    if isempty(title_params_order)
        title_params_order = keys(fixed_params)
    end

    title_params = [x for x in title_params_order if x in keys(fixed_params)]
    par_values_str = join([str_format_fun(par_val_printer(k),
                           par_val_adj(k, fixed_params[k]))
                      for k in title_params], "__")
    file_name = string("svm_", obj_fun_dict[Symbol(pt.svm_data[1, :obj_fun])],
                       "_", z_var, "_", par_values_str)

    return [graphs_path, file_name]
end

# ** HeatMap + Surface Plot {{{2
function svm_plot_heatmap_surf(pt, xy_list::Array{Symbol, 1},
                               z_var::Symbol,
                               fixed_params::Dict{Symbol,Float64};
                               title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                               save_fig::Bool=true,
                               elev::Float64=25.,
                               azim::Float64=210.,
                               zpad::Float64=10.,
                               plt_cmap::String="viridis", # 'Spectral', 'cool',
                               seaborn_style::String="darkgrid",
                               heat_reverse_x::Bool=false,
                               heat_reverse_y::Bool=true,
                               interp_bool::Bool=false,
                               smooth_bool::Bool=false,
                               ax1_dist::Float64=8.5,
                               axs_wspace::Float64=.1,
                               cbaxes::Array{Float64,1}=[.975, 0.15, 0.015, 0.675],
                               return_fig::Bool=true,
                               fig_dpi::Int64=400,
                               graph_format::String=heat_surf_graph_format)
# sup_title_x::Float64=.575,

    # ###################################################
    # ################### MULTI PLOT ####################
    # ###################################################
    using3D()
    fig = PyPlot.figure("pyplot_surfaceplot", figsize=PyPlot.figaspect(.4), facecolor="w")

    # ###################################################
    # ##################### SURFACE #####################
    # ###################################################
    ax1 = fig.add_subplot(1, 2, 2, projection="3d")

    xyz, grid_x, grid_y, grid_z, var_axes = svm_interp_smooth_surface(pt,
                                                                      fixed_params,
                                                                      z_var,
                                                                      interp=interp_bool,
                                                                      smooth=smooth_bool)

    # Notice that variables X and Y are listed in var_axes and follow the
    # convention of sigmah on the Y-Axis and Kappa on the X-Axis. When passing
    # the list of variables to the plot_surface function, use var_axes list:
    surf, ax1, fig_title = svm_plot_surface(pt, ax1, xyz[:,1], xyz[:,2], xyz[:,3],
                                            xy_list, z_var, fixed_params;
                                            title_params_order=title_params_order,
                                            plt_cmap=plt_cmap,
                                            seaborn_style=seaborn_style,
                                            make_title=true,
                                            add_title=false,
                                            zpad=zpad, view_elev=elev, view_azim=azim)

    # Set the background color of the pane YZ
    ax1.w_xaxis.set_pane_color(PyPlot.matplotlib.colors.hex2color("#d5d8dc"))

    # Add a color bar which maps values to colors.
    #    cb = fig.colorbar(surf, aspect=20, ax=ax1)
    cbax = fig.add_axes(cbaxes)
    cb = fig.colorbar(surf, aspect=20, cax=cbax)
    cb.outline.set_visible(false)

    ax1.patch.set_facecolor("white")
    ax1.dist = ax1_dist
    # ###################################################

    # ###################################################
    # ##################### HEATMAP #####################
    # ###################################################
    ax2 = fig.add_subplot(1, 2, 1)

    ax2, _ = svm_plot_heatmap(pt, xy_list, z_var, fixed_params;
                              plt_cmap=plt_cmap,
                              reverse_x_axis=heat_reverse_x,
                              reverse_y_axis=heat_reverse_y)
    # ###################################################

    fig.suptitle(fig_title, fontsize=14) #, x=sup_title_x)
    PyPlot.subplots_adjust(wspace=axs_wspace)
    #PyPlot.show()

    if save_fig
        m = unique(pt._svm_surf[:, :m])[1]
        folder_path, file_name = svm_heat_surf_plot_path_fname(pt, xy_list,
                                                               z_var, fixed_params;
                                                               title_params_order=title_params_order)

        PyPlot.savefig(string(folder_path, "/", file_name, ".", graph_format), dpi=fig_dpi, format=graph_format)
    end

    if return_fig
        return fig
    end
end


# * CVM v.s. SVM Plot Methods {{{1
# include("_cvm_vs_svm_plot_methods.jl")
function form_cvm_svm_combinations(cvm_m_comb_nums::Array{Int64,1},
                           svm_m_comb_nums::Array{Int64,1};
                           m::Float64=1.)
    id_vars = [:comb_num, :m_comb_num]

    sbt = get_bt(; model="svm", m=m, m_comb_num=1)
    cbt = get_bt(; model="cvm", m=m, m_comb_num=1)

    cvmls = [ ]
    for i in cvm_m_comb_nums
        tmp = cbt.bp.df[.&(cbt.bp.df[:m].==m, cbt.bp.df[:m_comb_num].== i), :]
        cols = [x for x in Symbol.(names(tmp)) if !(x in vcat(id_vars, [:lambda, :sigmah]))]
        push!(cvmls, Dict{Symbol, Float64}(cols .=> [tmp[1, var] for var in cols] ))
    end

    svmls = [ ]
    for i in svm_m_comb_nums
        tmp = sbt.bp.df[.&(sbt.bp.df[:m].==m, sbt.bp.df[:m_comb_num].== i), :]
        cols = [x for x in Symbol.(names(tmp)) if !(x in vcat(id_vars, :sigmah))]
        push!(svmls, Dict{Symbol, Float64}(cols .=> [tmp[1, var] for var in cols] ))
    end

    return Dict("svm" => svmls,
                "cvm" => cvmls)
end


function curve_legend(df_slice, var2::Symbol; model::String="svm",
                      xlabels::Dict{Symbol, Array{String,1}}=cvs_xlabels)

    var1_symbol = (model == "svm") ? Symbol(:kappa, :_ep) : Symbol(:kappa, :_otc)
    var2_val = (var2 == :iota) ? df_slice[1, var2] * 1e4 : df_slice[1, var2]

    return latexstring("\$\\left(", xlabels[var1_symbol][1], ", ",
                       xlabels[var2][1], "\\right) = \$ (",
                       df_slice[1, :kappa] * 1e4,
                       ", ", var2_val, ")")
end


function cvm_vs_svm_plot_dirs(yvar::Symbol, dvar::Symbol, xvar::Symbol;
                              m::Float64=NaN,
                              main_dir_path::String=main_dir_path,
                              plots_dir::String=plots_dir,
                              cvm_vs_svm_plots_dir::String=cvm_vs_svm_plots_dir)
    plots_path = string(main_dir_path, "/", plots_dir)
    fig_path = string(plots_path, "/", cvm_vs_svm_plots_dir)


    m_val = "all"
    if !isnan(m)
        m_val = str_format_fun(comb_folder_dict[:m][2], m)
    end
    fig_folder = string(fig_path, "/m_", m_val)

    # Create Directories
    for fdir in [plots_path, fig_path, fig_folder]
        if !isdir(fdir)
            mkdir(fdir)
        end
    end

    # FileName
    fig_name = string("cvm_vs_svm_", yvar, "_", dvar, "_", xvar,".png")

    return fig_folder, fig_name
end


function cvm_vs_svm_plotfun(cvmdf, svmdf,
                            xvar::Symbol, yvar::Symbol, dvar::Symbol;
                            figaspect::Float64=.5,
                            figpad::AbstractFloat=1.8,
                            plot_vlines::Bool=false,
                            xvars::Array{Symbol,1}=cvs_xvars,
                            cvmlinestyles::Array{String,1}=cvmlinestyles,
                            cvmmarkers::Array{String,1}=cvmmarkers,
                            svmlinestyles::Array{String,1}=svmlinestyles,
                            svmmarkers::Array{String,1}=svmmarkers,
                            save_fig::Bool=true,
                            fig_dpi::Int64=300,
                            main_dir_path::String=main_dir_path,
                            plots_dir::String=plots_dir,
                            cvm_vs_svm_plots_dir::String=cvm_vs_svm_plots_dir)


    # cvmdf = pt.cvm_data[[(x in cvm_combs) for x in pt.cvm_data[:comb_num]], :]
    # svmdf = pt.svm_data[[(x in svm_combs) for x in pt.svm_data[:comb_num]], :]

    # x-axis variable:
    # xvar =  [x for x in cvs_xvars if size(unique(svmdf[x]), 1)  == size(unique(pt.svm_data[x]), 1)][1]


    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=Tuple(PyPlot.figaspect(figaspect)))
    ax = fig.add_subplot(111)

    # Plot SVM Curves #############################################
    pchipflist = []
    xpos=0.
    i = 1
    for dval in unique(svmdf[:, dvar])
        # Slice SVM DataFrame
        svm_slice = svmdf[abs.(svmdf[:, dvar] .- dval) .< 1e-6, :]

        # Plot Curve
        ax.plot(svm_slice[:, xvar], svm_slice[:, yvar],
                color="blue",
                linewidth=1,
                linestyle=svmlinestyles[i],
                marker=svmmarkers[i],
                markersize=3)

        # Add Legend to the Curves
        svm_label = curve_legend(svm_slice, :lambda; model="svm")
        ax.text(svm_slice[end, xvar], svm_slice[end, yvar],
                  svm_label, fontsize=10, va="bottom")


        # Interpolate SVM Curves and Store Them
        pchipf = Dierckx.Spline1D(svm_slice[:, xvar], svm_slice[:, yvar];
                                  k=3, bc="extrapolate")
        push!(pchipflist, pchipf)


        if dval == unique(svmdf[:, dvar])[end]
            xvals = range(minimum(svm_slice[:, xvar]),
                          stop=maximum(svm_slice[:, xvar]), length=10^4)
            xpos = svm_slice[end, xvar]
        end

        i += 1
    end
    # #############################################################


    # Plot CVM Curves #############################################
    for dval in unique(cvmdf[:, dvar])
        # Slice CVM DataFrame
        if !isnan(dval)
            cvm_slice = cvmdf[abs.(cvmdf[:, dvar] .- dval) .< 1e-6, :]
        else
            cvm_slice = cvmdf
        end

        # Plot Curve
        ax.axhline(cvm_slice[1, yvar],
                   color="green",
                   linewidth=1,
                   linestyle=cvmlinestyles[i],
                   marker=cvmmarkers[i])

        # Add Legend to the Curves
        cvm_label = curve_legend(cvm_slice, :iota; model="cvm")
        ax.text(xpos,
                cvm_slice[1, yvar],
                cvm_label, fontsize=10, va="bottom")


        # Plot Vertical Lines
        if plot_vlines
            for pchipf in pchipflist
                sigstar = xvals[argmin([abs(pchipf(x) - cvm_slice[1, yvar]) for x in xvals])]
                ax.axvline(sigstar, color="black", linewidth=.5, linestyle="-.")
            end
        end
    end
    # #############################################################


    # Label Axes
    ax.set_xlabel(string("\$", cvs_xlabels[xvar][1], "\$"), labelpad=10)
    ax.set_ylabel(cvs_ylabels[yvar], labelpad=10)

    # Set Title
    title_params = join([string("\$", tlabels[x][1], "= \$ ",
                                str_format_fun(tlabels[x][2], svmdf[1, x]))
                         for x in cvs_plots_title_params_order], ", ")
    plot_title = latexstring("Optimal RMP-Conditional ", cvs_ylabels[yvar], " for ", title_params)
    fig.suptitle(plot_title, fontsize=14)
    ax.set_title("(\$\\kappa\$ and \$\\iota\$ values in b.p.)", fontsize=12)

    if !isnan(figpad)
        fig.tight_layout(pad=figpad)
    end

    if save_fig
        fig_folder, fig_name = cvm_vs_svm_plot_dirs(yvar, dvar, xvar;
                                                    m=svmdf[1, :m],
                                                    main_dir_path=main_dir_path,
                                                    plots_dir=plots_dir,
                                                    cvm_vs_svm_plots_dir=cvm_vs_svm_plots_dir)
        plt.savefig(string(fig_folder, "/", fig_name), dpi=fig_dpi, bbox_inches="tight")
    end
#     display(fig)
    return fig
end


# ** Trash {{{2
# ##############################################################################
# ################################### TRASH ####################################
# ##############################################################################
function df_slicer(pt, i::Int64, combs; model::String="svm")
    tmp = combs[i]
    if model == "svm"
        svmloc = sum([abs.(pt.svm_data[:, x] .- tmp[x]) .< 1e-4
                      for x in keys(tmp)]) .== length(keys(tmp))
        return pt.svm_data[svmloc, :]
        # sort!(pt.svm_data[svmloc, :], xvar)
    else
        cvmloc = sum([abs.(pt.cvm_data[:, x] .- tmp[x]) .< 1e-4
                      for x in keys(tmp)]) .== length(keys(tmp))
        return pt.cvm_data[cvmloc, :]
    end
end


# * Risk-Management Policy Choice {{{1

# ** RMP Auxiliary {{{2
# include("_RMP/_rmp_auxiliary.jl")
function fig_size_pad_adjuster(subplots::Int64;
                               figaspect::Float64=NaN,
                               figsize::Tuple{Float64, Float64}=(.0, .0),
                               figpad::Float64=1.8)

    figaspect = isnan(figaspect) ? rmp_fig_aspect : figaspect
    if subplots == 1
        figsize = sum(figsize) == .0 ? Tuple(PyPlot.figaspect(rmp_fig_aspect)) : figsize
    else
        figsize = sum(figsize) == .0 ? rmp_multi_plot_fig_size : figsize
        figpad = rmp_multi_plot_figpad
    end

    return figsize, figpad
end


function get_cutoff_value(df::DataFrame, xvar::Symbol, yvar::Symbol, yval::Float64;
                          xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                              Base.TwicePrecision{Float64}}=range(.0, stop=.0, length=0))

    if all(df[:, yvar] .> yval)
        return Inf
    elseif all(df[:, yvar] .< yval)
        return -Inf
    else
        if size(xgrid, 1) == 0
            xgrid = range(minimum(df[:, xvar]), stop=maximum(df[:, xvar]), length=10^5)
        end

        yinterp = Dierckx.Spline1D(df[:, xvar], df[:, yvar];
                               k=3, bc="extrapolate")
        return xgrid[argmin(abs.(yinterp(xgrid) .- yval))]
    end
end


function get_misrep_cutoff_value(xvar::Symbol, yvar::Symbol,
                                 cvmdf::DataFrame,
                                 svmdf::DataFrame,
                                 misrepdf::DataFrame;
                                 xgrid::StepRangeLen{Float64,
                                                     Base.TwicePrecision{Float64},
                                                     Base.TwicePrecision{Float64}}=range(.0,
                                                                                         stop=.0,
                                                                                         length=0))

    xvals = svmdf[:, xvar]
    if size(xgrid, 1) == 0
        xgrid = range(minimum(xvals), stop=maximum(xvals), length=10^5)
    end

    misrep_yval_interp = Dierckx.Spline1D(misrepdf[:, Symbol(:r_, xvar)],
                                          misrepdf[:, Symbol(:r_, yvar)]; k=3, bc="extrapolate")
    if xvar != :sigmah
        fi_yvals = [maximum([x, cvmdf[1, yvar]]) for x in svmdf[:, yvar]]
        if xvar == :lambda
            fi_yvals = svmdf[:, yvar]
        end
        fi_yval_interp = Dierckx.Spline1D(xvals, fi_yvals; k=3, bc="extrapolate")

        return xgrid[argmin(abs.(fi_yval_interp(xgrid) .- misrep_yval_interp(xgrid)))]
    else
        fi_yval_interp = Dierckx.Spline1D(xvals, svmdf[:, yvar]; k=3, bc="extrapolate")

        svm_cv = NaN
        svm_diff = abs.(fi_yval_interp(xgrid) .- misrep_yval_interp(xgrid))
        if !isempty(svm_diff .< 1e-5)
            svm_cv = xgrid[argmin(svm_diff)]
        end

        cvm_cv = NaN
        cvm_diff = abs.(misrep_yval_interp(xgrid) .- cvmdf[1, yvar])
        if !isempty(cvm_diff .< 1e-5)
            cvm_cv = xgrid[argmin(cvm_diff)]
        end
        return cvm_cv, svm_cv
    end

end


function rmp_plot_dirs(yvars::Array{Symbol, 1}, xvar::Symbol;
                       nrmp::Bool=false,
                       m::Float64=NaN,
                       rf_model::String="svm",
                       main_dir_path::String=main_dir_path,
                       plots_dir::String=plots_dir,
                       nrmp_plots_dir::String=nrmp_plots_dir,
                       rmp_plots_dir::String=rmp_plots_dir,
                       misrep::Bool=false,
                       fname_ext::String=rmp_fname_ext)
    plots_path = string(main_dir_path, "/", plots_dir)
    if nrmp
        fig_path = string(plots_path, "/", nrmp_plots_dir)
    else
        fig_path = string(plots_path, "/", rmp_plots_dir)
    end

    m_val = "all"
    if !isnan(m)
        m_val = str_format_fun(comb_folder_dict[:m][2], m)
    end
    fig_folder = string(fig_path, "/m_", m_val)

    # Create Directories
    for fdir in [plots_path, fig_path, fig_folder]
        if !isdir(fdir)
            mkdir(fdir)
        end
    end

    # FileName
    yvars_fig_name = join([obj_fun_dict[y] for y in yvars], "_")

    type_prefix = full_info_prefix
    if misrep
        type_prefix = misrep_prefix
    end

    fn_prefix = nrmp ? nrmp_fn_prefix : rmp_fn_prefix
    if rf_model == "svm"
        fig_name = string(fn_prefix, "_", type_prefix, "_",
                          yvars_fig_name, "_", xvar,".", fname_ext)
    elseif rf_model == "cvm"
        fig_name = string("cvm_", fn_prefix, "_", type_prefix, "_",
                          yvars_fig_name, "_", xvar,".", fname_ext)
    else
        println("Wrong Risky Firm Model. Please enter 'cvm' or 'svm'. Exiting...")
        return
    end

    return fig_folder, fig_name
end


# ** RMP Curves {{{2
# include("_RMP/_rmp_curves.jl")
function plot_cvm_curve(ax, xvar::Symbol, yvar::Symbol,
                        sfdf::DataFrame,
                        text_xloc::Float64, text_yloc::Float64;
                        nrmp::Bool=false,
                        xgrid::StepRangeLen{Float64,
                                            Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}}=range(.0, stop=.0, length=0))

    if size(sfdf, 1) == 1
        ax.axhline(sfdf[1, yvar],
                   color=cvm_curve_color,
                   linewidth=1,
                   linestyle=cvmlinestyles[1],
                   marker=cvmmarkers[1])
    elseif size(xgrid, 1) == 0
        ax.plot(sfdf[:, xvar], sfdf[:, yvar];
                color=cvm_curve_color,
                linewidth=1,
                linestyle=cvmlinestyles[1],
                marker=cvmmarkers[1],
                markersize=3)
    else
        y_interp = Dierckx.Spline1D(sfdf[:, xvar], sfdf[:, yvar]; k=3, bc="extrapolate")
        ax.plot(xgrid, y_interp(xgrid);
                color=cvm_curve_color,
                linewidth=1,
                linestyle=cvmlinestyles[1])
    end


    # Add Legend to the Curve
    if nrmp
        if xvar == :sigmah
            cvm_label = latexstring("\$", cvs_xlabels[:lambda][1], "= 0.\$")
        else
            text_xloc = .95 * text_xloc
            cvm_label = latexstring("\$", cvs_xlabels[:sigmah][1], "= \\sigma_l = " , sfdf[1, :sigmal], "\$")
        end
    elseif xvar == :iota
        cvm_label = latexstring("\$", cvs_xlabels[:iota][1], "=", cvm_curve_label, "\$")
    else #if xvar == :sigmah
        cvm_label = latexstring("\$", cvs_xlabels[:iota][1], "=",
                                str_format_fun(cvs_xlabels[:iota][2], sfdf[1, :iota]), "\$ (b.p.)")
    end

    ax.text(text_xloc, text_yloc,
            cvm_label, fontsize=10, va="bottom")

    return ax
end


function plot_svm_curve(ax, xvar::Symbol, yvar::Symbol,
                        rfdf::DataFrame,
                        text_xloc::Float64, text_yloc::Float64,
                        xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}};
                        nrmp::Bool=false)

    if size(rfdf, 1) == 1
        ax.axhline(rfdf[1, yvar],
                   color=svm_curve_color,
                   linewidth=1,
                   linestyle=svmlinestyles[1],
                   marker=svmmarkers[1])
    elseif any([size(xgrid, 1) == 0, size(rfdf, 1) <= 3])
        ax.plot(rfdf[:, xvar], rfdf[:, yvar];
                color=svm_curve_color,
                linewidth=1,
                linestyle=svmlinestyles[1],
                marker=svmmarkers[1],
                markersize=3)
    else
        y_interp = Dierckx.Spline1D(rfdf[:, xvar], rfdf[:, yvar]; k=3, bc="extrapolate")
        ax.plot(xgrid, y_interp(xgrid);
                color=svm_curve_color,
                linewidth=1,
                linestyle=svmlinestyles[1])
    end


    # Add Legend to the Curve
    if nrmp
        if xvar == :sigmah
            svm_label = latexstring("\$",
                                    ModelPlots.cvs_xlabels[:lambda][1],
                                    "\$ = ",
                                    [x for x in unique(rfdf[:, :lambda]) if !isnan(x)][1])
        else
            svm_label = latexstring("\$",
                                ModelPlots.cvs_xlabels[:sigmah][1],
                                "\$ = ",
                                [x for x in unique(rfdf[:, :sigmah]) if !isnan(x)][1])
        end
    elseif xvar == :iota
        if !isnan(rfdf[1, :sigmah])
            svm_label = latexstring("(\$",
                                    ModelPlots.cvs_xlabels[:iota][1], ", ",
                                    ModelPlots.cvs_xlabels[:lambda][1], ", ",
                                    ModelPlots.cvs_xlabels[:sigmah][1],
                                    "\$) = (",
                                    rfdf[1, :iota], ", ",
                                    rfdf[1, :lambda], ", ",
                                    rfdf[1, :sigmah], ")")
        else
            svm_label = latexstring("(\$",
                                    ModelPlots.cvs_xlabels[:iota][1], ", ",
                                    ModelPlots.cvs_xlabels[:sigmah][1],
                                    "\$) = (",
                                    rfdf[1, :iota], ", ",
                                    rfdf[1, :sigmal], ")")
        end
    elseif xvar == :sigmah
        if !isnan(rfdf[1, :lambda])
            svm_label = latexstring("(\$",
                                    ModelPlots.cvs_xlabels[:iota][1], ", ",
                                    ModelPlots.cvs_xlabels[:lambda][1],
                                    "\$) = (",
                                    rfdf[1, :iota], ", ",
                                    rfdf[1, :lambda], ")")
        else
            svm_label = latexstring("\$",
                                    ModelPlots.cvs_xlabels[:iota][1],
                                    "=",
                                    str_format_fun(cvs_xlabels[:iota][2],
                                                   rfdf[1, :iota]), "\$ (b.p.)")
        end
    elseif xvar == :lambda
        svm_label = latexstring("(\$",
                                ModelPlots.cvs_xlabels[:iota][1], ", ",
                                ModelPlots.cvs_xlabels[:sigmah][1],
                                "\$) = (",
                                rfdf[1, :iota], ", ",
                                rfdf[1, :sigmah], ")")
    end

    ax.text(text_xloc,  text_yloc,
            svm_label, fontsize=10, va="bottom")

    return ax
end


function plot_misrep_curve(ax, xvar::Symbol, yvar::Symbol,
                           misrepdf::DataFrame;
                           xgrid::StepRangeLen{Float64,
                                            Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}}=range(.0, stop=.0, length=0))

    #  text_xloc::Float64, text_yloc::Float64)

    if xvar == :iota
        ax.axhline(misrepdf[1, Symbol(:r_, yvar)],
                   color=misrep_curve_color,
                   linewidth=1,
                   linestyle=svmlinestyles[1],
                   marker=svmmarkers[1])
    else
        if !in(misrepdf[1, Symbol(:s_, xvar)], misrepdf[:, Symbol(:r_, xvar)])
            cols = [vcat([misrepdf[1, Symbol(:s_, xvar)], misrepdf[:, Symbol(:r_, xvar)]...]),
                    vcat([misrepdf[1, Symbol(:s_, yvar)], misrepdf[:, Symbol(:r_, yvar)]...])]
        else
            cols = [misrepdf[:, Symbol(:r_, xvar)], misrepdf[:, Symbol(:r_, yvar)]]
        end

        tmp = sort!(DataFrame(Dict(zip([xvar, yvar], cols))), xvar)
        if size(xgrid, 1) > 0
            y_interp = Dierckx.Spline1D(tmp[:, xvar], tmp[:, yvar]; k=3, bc="extrapolate")
            ax.plot(xgrid, y_interp(xgrid),
                    color=misrep_curve_color,
                    linewidth=1,
                    linestyle=svmlinestyles[1],
                    marker=svmmarkers[1],
                    markersize=3)
        else
            ax.plot(tmp[:, xvar], tmp[:, yvar];
                    color=misrep_curve_color,
                    linewidth=1,
                    linestyle=svmlinestyles[1],
                    marker=svmmarkers[1],
                    markersize=3)
        end
        # ax.plot(misrepdf[Symbol(:r_,xvar)], misrepdf[Symbol(:r_, yvar)];
        #         color=misrep_curve_color,
        #         linewidth=1,
        #         linestyle=svmlinestyles[1],
        #         marker=svmmarkers[1],
        #         markersize=3)
    end

    # Add Legend to the Curve
    # svm_label = latexstring("(\$",
    #                         ModelPlots.cvs_xlabels[:iota][1], ", ",
    #                         ModelPlots.cvs_xlabels[:lambda][1], ", ",
    #                         ModelPlots.cvs_xlabels[:sigmah][1],
    #                         "\$) = (",
    #                         rfdf[1, :iota], ", ",
    #                         rfdf[1, :lambda], ", ",
    #                         rfdf[1, :sigmah], ")")

    # ax.text(text_xloc,  text_yloc,
    #         svm_label, fontsize=10, va="bottom")

    return ax
end


function plot_vlines(ax, xvar;
                     fv_xvar::Float64=NaN,
                     fv_color::String=fv_color,
                     mbr_xvar::Float64=NaN,
                     mbr_color::String=mbr_color,
                     cvm_misrep_xvar::Float64=NaN,
                     svm_misrep_xvar::Float64=NaN,
                     misrep_color::String=misrep_color)

    xvar = (xvar == :sigmal) ? :sigmah : xvar

    if abs(cvm_misrep_xvar - fv_xvar) < 1e-5
            cvm_misrep_xvar = NaN
    end

    # Form Dictionary with labels and values:
    vldict = vlines_labels_dict(xvar; fv_xvar=fv_xvar,
                                fv_color=fv_color,
                                mbr_xvar=mbr_xvar,
                                mbr_color=mbr_color,
                                cvm_misrep_xvar=cvm_misrep_xvar,
                                svm_misrep_xvar=svm_misrep_xvar,
                                misrep_color=misrep_color)


    xkeys = [x for x in keys(vldict) if .&(!isnan(vldict[x][:value]), !isinf(vldict[x][:value]))]
    minor_ticks = [vldict[x][:value] for x in xkeys]
    minor_labels = [vldict[x][:xsym] for x in xkeys]
    for x in xkeys
        ax.axvline(vldict[x][:value],
                   color=vldict[x][:color],
                   linewidth=.6,
                   linestyle="--",
                   marker=svmmarkers[1])
    end

    ax.set_xticks(minor_ticks, minor=true)
    ax.set_xticklabels(minor_labels, minor=true)

    return ax
end


# ** RMP Color Regions {{{2
# include("_RMP/_rmp_color_regions.jl")
function color_optimal_rm_region_fun(ax, fv_xvar::Float64,
                                     xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                         Base.TwicePrecision{Float64}},
                                     text_xloc::Float64,
                                     text_yloc::Float64)

    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .>= fv_xvar,
                    facecolor=rm_region_color, alpha=0.15)

    if !isinf(text_xloc)
        ax.text(text_xloc, text_yloc,
                "Risk-Management \n is optimal",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
                style="italic",
                bbox=Dict("facecolor" => box_color, "alpha" => 0.5, "pad" => 10))
    end


    return ax
end


function color_optimal_nrm_region_fun(ax, fv_xvar::Float64,
                                      xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                          Base.TwicePrecision{Float64}},
                                      text_xloc::Float64,
                                      text_yloc::Float64)
    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .<= fv_xvar,
                    facecolor=nrm_region_color, alpha=0.25)

    if !isinf(text_xloc)
        ax.text(text_xloc, text_yloc,
                "No Risk-Management is optimal",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
                style="italic",
                bbox=Dict("facecolor" => box_color, "alpha" => 0.5, "pad" => 10))
    end

    return ax
end


function color_conflict_region_fun(ax, xvar::Symbol, fv_xvar::Float64, mbr_xvar::Float64,
                                   xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                       Base.TwicePrecision{Float64}})

    if xvar == :iota
        region_cond = .&(xgrid .>= fv_xvar, xgrid .<= mbr_xvar)
    else #if xvar == :sigmah
        region_cond = .&(xgrid .>= minimum([fv_xvar, mbr_xvar]),
                         xgrid .<= maximum([fv_xvar, mbr_xvar]))
    end


    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData,
                                                                   ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=region_cond,
                    facecolor=conflict_region_color, alpha=0.25)

    return ax
end


function color_misrep_region_fun(ax,
                                 cvm_misrep_xvar::Float64,
                                 svm_misrep_xvar::Float64,
                                 xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                     Base.TwicePrecision{Float64}},
                                 text_xloc::Float64,
                                 text_yloc::Float64)
    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)

    if any([isnan(cvm_misrep_xvar), isnan(svm_misrep_xvar)]) || abs(cvm_misrep_xvar - svm_misrep_xvar) < 1e-5
        misrep_xvar = isnan(cvm_misrep_xvar) ? svm_misrep_xvar : cvm_misrep_xvar
        ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .>= misrep_xvar,
                        facecolor=misrep_region_color, alpha=0.25)
    else
        ax.fill_between(xgrid, 0, 1, transform=trans,
                        where= .&(xgrid .>= minimum([cvm_misrep_xvar, svm_misrep_xvar]),
                                  xgrid .<= maximum([cvm_misrep_xvar, svm_misrep_xvar])),
                        facecolor=misrep_region_color, alpha=0.25)
    end

    ax.text(text_xloc, text_yloc,
            "Misrepresentation is optimal",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
            style="italic",
            bbox=Dict("facecolor" => box_color, "alpha" => 0.5, "pad" => 10))

    return ax
end


function color_regions_fun(ax, xvar::Symbol,
                           xmin::Float64, xmax::Float64,
                           ymin::Float64, ymax::Float64,
                           xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                               Base.TwicePrecision{Float64}};
                           fv_xvar::Float64=NaN,
                           mbr_xvar::Float64=NaN,
                           cvm_misrep_xvar::Float64=NaN,
                           svm_misrep_xvar::Float64=NaN,
                           color_rm_region::Bool=true,
                           color_nrm_region::Bool=true,
                           color_conflict_region::Bool=false,
                           color_misrep_region::Bool=false)

    # # Axes Limits ##################################################
    # xmin, xmax = ax.get_xlim()
    # ymin, ymax = ax.get_ylim()

    # println(string("xmin: ", xmin))
    # println(string("ymin: ", ymin))
    # println(string("xmax: ", xmax))
    # println(string("ymax: ", ymax))
    # # ##############################################################


    # Optimal Risk-Management Region under Full Information
    if .&(!isnan(fv_xvar), color_rm_region)
        if xvar == :iota
            xloc =  fv_xvar / 2.
            yloc = .8 * ymin + .2 * ymax
        elseif xvar == :sigmah
            xloc = (xmax + fv_xvar) / 2
            yloc = .2 * ymin + .8 * ymax
        elseif xvar == :lambda
            xloc = (xmax + fv_xvar) / 2
            yloc = .4 * ymin + .6 * ymax
        end

        ax = color_optimal_rm_region_fun(ax, fv_xvar, xgrid, xloc, yloc)
    end

    # Optimal No-Risk-Management Region under Full Information
    if .&(!isnan(fv_xvar), color_nrm_region)
        if xvar == :iota
            xloc = .5 * fv_xvar + .5 * xmax
            yloc = .1 * ymin + .9 * ymax
        elseif xvar == :sigmah
            xloc = .5 * fv_xvar + .5 * xmin
            yloc = .85 * ymin + .15 * ymax
        elseif xvar == :lambda
            xpos = !isinf(fv_xvar) ? fv_xvar : xmax

            xloc = .5 * xmax + .5 * xmin
            yloc = .4 * ymin + .6 * ymax
        end

        ax = color_optimal_nrm_region_fun(ax, fv_xvar, xgrid, xloc, yloc)
    end

    # Conflict Region under Full Information:
    # RM Firm Value >= NRM Firm Value, but RM MBR <= NRM MBR
    if .&(!isnan(fv_xvar), !isnan(mbr_xvar), color_conflict_region)
        ax = color_conflict_region_fun(ax, xvar, fv_xvar, mbr_xvar, xgrid)
    end

    # Misrepresentation Region
    # Misrepresentation MBR > FI MBR
    if .&(any([!isnan(cvm_misrep_xvar), !isnan(svm_misrep_xvar)]), color_misrep_region)
        if any([isnan(cvm_misrep_xvar), isnan(svm_misrep_xvar)]) || abs(cvm_misrep_xvar - svm_misrep_xvar) < 1e-5
            misrep_xvar = isnan(cvm_misrep_xvar) ? svm_misrep_xvar : cvm_misrep_xvar
            xloc = .5 * misrep_xvar + .5 * xmax
        else
            xloc = .5 * ( cvm_misrep_xvar + svm_misrep_xvar)
        end

        yloc = .1 * ymin + .9 * ymax
        ax = color_misrep_region_fun(ax, cvm_misrep_xvar, svm_misrep_xvar, xgrid, xloc, yloc)
    end

    return ax
end


# ** RMP Plot Methods {{{2
# include("_RMP/_rmp_fi_plot_methods.jl")
function rmp_subplotfun(fig::Figure, xvar::Symbol,
                        yvar::Symbol,
                        sfdf::DataFrame,
                        rfdf::DataFrame,
                        xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}};
                        nrmp::Bool=false,
                        ax_subplot::Int64=111,
                        interp_yvar::Bool=false,
                        misrepdf::DataFrame=DataFrame(),
                        fv_xvar::Float64=NaN,
                        mbr_xvar::Float64=NaN,
                        cvm_misrep_xvar::Float64=NaN,
                        svm_misrep_xvar::Float64=NaN,
                        color_rm_region::Bool=true,
                        color_nrm_region::Bool=true,
                        color_conflict_region::Bool=false,
                        color_misrep_region::Bool=false,
                        cvmlinestyles::Array{String,1}=cvmlinestyles,
                        cvmmarkers::Array{String,1}=cvmmarkers,
                        svmlinestyles::Array{String,1}=svmlinestyles,
                        svmmarkers::Array{String,1}=svmmarkers)

    ax = fig.add_subplot(ax_subplot)

    # Plot Non-Horizontal Curve ####################################
    if xvar == :iota
        # Plot CVM Curve
        xloc = (sfdf[end - 1, xvar] + sfdf[end, xvar])/2
        yloc = (sfdf[end - 1, yvar] + sfdf[end, yvar])/2

        if interp_yvar
            ax = plot_cvm_curve(ax, xvar, yvar, sfdf, xloc, yloc;
                                xgrid=xgrid, nrmp=nrmp)
        else
            ax = plot_cvm_curve(ax, xvar, yvar, sfdf, xloc, yloc; nrmp=nrmp)
        end
    else
        # Plot SVM Curve
        xloc = .95 * rfdf[end, xvar]
        yloc = .1 * (maximum(rfdf[:, yvar]) - minimum(rfdf[:, yvar])) + minimum(rfdf[:, yvar])
        ax = plot_svm_curve(ax, xvar, yvar, rfdf, xloc, yloc, xgrid; nrmp=nrmp)
        # ##############################################################
    end
    # ##############################################################


    # Axes Limits ##################################################
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # ##############################################################


    # Plot Horizontal Curve ############################################
    if xvar == :iota
        # Plot SVM Curve ###############################################
        # Because SVM is now a horizontal curve, it must come after
        # the CVM curve for legend placement purposes (x axis)

        # Plot Curve
        xloc = .725 * xmax
        yloc = rfdf[end, yvar]
        if .&(!isempty(misrepdf), yvar == :MBR)
            yloc = yloc - .125 * (rfdf[end, yvar] - ymin)
        end
        ax = plot_svm_curve(ax, xvar, yvar, rfdf, xloc, yloc, xgrid; nrmp=nrmp)
        # ##############################################################
    else
        xloc = .95 * xmax
        yloc = sfdf[1, yvar]
        ax = plot_cvm_curve(ax, xvar, yvar, sfdf, xloc, yloc; nrmp=nrmp)
    end
    # ###################################################################

    # Plot Misrepresentation Curve #################################
    if .&(!isempty(misrepdf), (Symbol(:r_, yvar) in Symbol.(names(misrepdf))))
        ax = plot_misrep_curve(ax, xvar, yvar, misrepdf; xgrid=xgrid)
    end
    # ##############################################################

    # Axes Limits ##################################################
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # ##############################################################


    # Vertical Lines ###############################################
    ax = plot_vlines(ax, xvar; fv_xvar=fv_xvar,
                     mbr_xvar=mbr_xvar,
                     cvm_misrep_xvar=cvm_misrep_xvar,
                     svm_misrep_xvar=svm_misrep_xvar)
    # ##############################################################


    # Color Regions ################################################
    ax = color_regions_fun(ax, xvar,
                           xmin, xmax,
                           ymin, ymax,
                           xgrid;
                           fv_xvar=fv_xvar,
                           mbr_xvar=mbr_xvar,
                           cvm_misrep_xvar=cvm_misrep_xvar,
                           svm_misrep_xvar=svm_misrep_xvar,
                           color_rm_region=color_rm_region,
                           color_nrm_region=color_nrm_region,
                           color_conflict_region=color_conflict_region,
                           color_misrep_region=color_misrep_region)
    # ##############################################################


   # Axes Limits ##################################################
    #xmin, xmax = ax.get_xlim()
    #ymin, ymax = ax.get_ylim()
    # ##############################################################


    # For some reason, y limits are not matching axis ylim.
    # Force reset:
    ax.set_xlim([minimum(xgrid), maximum(xgrid)])
    ax.set_ylim([ymin, ymax])

    return ax
end


function rmp_core_plot(fig, xvar::Symbol, yvars::Array{Symbol,1},
                       sfdf::DataFrame, rfdf::DataFrame,
                       xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                              Base.TwicePrecision{Float64}},
                       subplots::Int64;
                       nrmp::Bool=false,
                       interp_yvar::Bool=false,
                       misrepdf::DataFrame=DataFrame(),
                       fv_xvar::Float64=NaN,
                       mbr_xvar::Float64=NaN,
                       cvm_misrep_xvar::Float64=NaN,
                       svm_misrep_xvar::Float64=NaN,
                       cvmlinestyles::Array{String,1}=cvmlinestyles,
                       cvmmarkers::Array{String,1}=cvmmarkers,
                       svmlinestyles::Array{String,1}=svmlinestyles,
                       svmmarkers::Array{String,1}=svmmarkers,
                       color_rm_region::Bool=true,
                       color_nrm_region::Bool=true,
                       color_conflict_region::Bool=false,
                       color_misrep_region::Bool=false)

    axes = []
    count = 1
    if xvar == :iota
        xvar_xlabel = string("Risk Management Cost \$",
                             cvs_xlabels[xvar][1], "\$ (b.p.)")
    elseif xvar == :sigmah
        xvar_xlabel = string("Post-Shock Volatility \$",
                             cvs_xlabels[xvar][1], "\$")
    elseif xvar == :lambda
        xvar_xlabel = string("Shock Intensity \$",
                             cvs_xlabels[xvar][1], "\$")
    end

    for yvar in yvars
        ax_subplot = ax_subplots[subplots][count]

        ax = rmp_subplotfun(fig, xvar, yvar,
                            sfdf, rfdf, xgrid;
                            nrmp=nrmp,
                            ax_subplot=ax_subplot,
                            interp_yvar=interp_yvar,
                            misrepdf=misrepdf,
                            fv_xvar=fv_xvar,
                            mbr_xvar=mbr_xvar,
                            cvm_misrep_xvar=cvm_misrep_xvar,
                            svm_misrep_xvar=svm_misrep_xvar,
                            cvmlinestyles=cvmlinestyles,
                            cvmmarkers=cvmmarkers,
                            svmlinestyles=svmlinestyles,
                            svmmarkers=svmmarkers,
                            color_rm_region=color_rm_region,
                            color_nrm_region=color_nrm_region,
                            color_conflict_region=color_conflict_region,
                            color_misrep_region=color_misrep_region)


        # ##############################################################
        # ################### Plot Labels and Titles ###################
        # ##############################################################
        if subplots == 1
            ax.set_xlabel(xvar_xlabel, labelpad=10)
            ax.set_title(" ")
        elseif count == 1
            ax.set_title(" \n ")
        end
        ax.set_ylabel(cvs_ylabels[yvar], labelpad=10)
        # ##############################################################

        # ##############################################################
        # ########################## Safe Type #########################
        # ##############################################################
        if !isempty(misrepdf)
            misrep_xval = misrepdf[1, Symbol(:s_, xvar)]
            misrep_yval = misrepdf[1, Symbol(:s_, yvar)]
            ax.scatter(misrep_xval,
                       misrep_yval;
                       s=25,color="purple", alpha=0.8)
            ax.text(misrep_xval, 1.0025 * misrep_yval, "\$S\$")
        end
        # ##############################################################

        push!(axes, ax)
        count += 1
    end

    if subplots > 1
        PyPlot.matplotlib.pyplot.xlabel(xvar_xlabel,
                                        labelpad=10)
    end


    return fig
end


# ATTENTION: kappa and iota values should be in basis points!
function rmp_fi_plotfun(xvar::Symbol, yvars::Array{Symbol,1},
                        sfdf::DataFrame,
                        rfdf::DataFrame;
                        nrmp::Bool=false,
                        xgrid::StepRangeLen{Float64,
                                            Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}}=range(.0,
                                                                                stop=0.,
                                                                                length=0),
                        interp_yvar::Bool=false,
                        misrepdf::DataFrame=DataFrame(),
                        fv_xvar::Float64=NaN,
                        mbr_xvar::Float64=NaN,
                        cvm_misrep_xvar::Float64=NaN,
                        svm_misrep_xvar::Float64=NaN,
                        color_rm_region::Bool=true,
                        color_nrm_region::Bool=true,
                        color_conflict_region::Bool=false,
                        color_misrep_region::Bool=false,
                        cvmlinestyles::Array{String,1}=cvmlinestyles,
                        cvmmarkers::Array{String,1}=cvmmarkers,
                        svmlinestyles::Array{String,1}=svmlinestyles,
                        svmmarkers::Array{String,1}=svmmarkers,
                        figaspect::Float64=NaN,
                        figsize::Tuple{Float64, Float64}=(.0, .0),
                        figpad::Float64=1.8,
                        save_fig::Bool=true,
                        fig_dpi::Int64=300,
                        main_dir_path::String=main_dir_path,
                        plots_dir::String=plots_dir,
                        rmp_plots_dir::String=rmp_plots_dir)

    # Figure Size and Layout Aspect
    subplots = size(yvars, 1)
    figsize, figpad = fig_size_pad_adjuster(subplots;
                                            figaspect=figaspect,
                                            figsize=figsize,
                                            figpad=figpad)


    if size(xgrid, 1) == 0
        if size(sfdf, 1) > 1
            xgrid = range(minimum(sfdf[:, xvar]), stop=maximum(sfdf[:, xvar]), length=10^5)
        else
            xgrid = range(minimum(rfdf[:, xvar]), stop=maximum(rfdf[:, xvar]), length=10^5)
        end
    end

    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=figsize)
    fig = rmp_core_plot(fig, xvar, yvars,
                        sfdf, rfdf, xgrid,
                        subplots;
                        nrmp=nrmp,
                        interp_yvar=interp_yvar,
                        misrepdf=misrepdf,
                        fv_xvar=fv_xvar,
                        mbr_xvar=mbr_xvar,
                        cvm_misrep_xvar=cvm_misrep_xvar,
                        svm_misrep_xvar=svm_misrep_xvar,
                        cvmlinestyles=cvmlinestyles,
                        cvmmarkers=cvmmarkers,
                        svmlinestyles=svmlinestyles,
                        svmmarkers=svmmarkers,
                        color_rm_region=color_rm_region,
                        color_nrm_region=color_nrm_region,
                        color_conflict_region=color_conflict_region,
                        color_misrep_region=color_misrep_region)


    # Set Sup Title
    suptitle_yvars = join([cvs_ylabels[yvar] for yvar in yvars], " and ")
    if .&(!isnan(rfdf[1, :sigmah]), !isnan(rfdf[1, :lambda]))
        suptitle_params = join([string("\$", tlabels[x][1], "= \$ ",
                                       str_format_fun(ModelPlots.tlabels[x][2], rfdf[1, x]))
                                for x in rmp_plots_title_params_order], ", ")
    else
        suptitle_params = join([string("\$", tlabels[x][1], "= \$ ",
                                       str_format_fun(ModelPlots.tlabels[x][2], x!=:sigmal ? rfdf[1, x] : sfdf[1, x]))
                                for x in rmp_plots_title_params_order],
                                ", ")
    end

    plot_suptitle = latexstring(suptitle_yvars, " under Full Information \n",
                             " for ", suptitle_params)
    fig.suptitle(plot_suptitle, fontsize=14)


    if !isnan(figpad)
        fig.tight_layout(pad=figpad)
    end


    if save_fig
        rf_model = isnan(rfdf[end, :lambda]) ? "cvm" : "svm"
        fig_folder, fig_name = rmp_plot_dirs(yvars, xvar;
                                             nrmp=nrmp,
                                             m=rfdf[1, :m],
                                             rf_model=rf_model,
                                             main_dir_path=main_dir_path,
                                             plots_dir=plots_dir,
                                             rmp_plots_dir=rmp_plots_dir,
                                             misrep=!isempty(misrepdf))
        PyPlot.savefig(string(fig_folder, "/", fig_name), dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end


# * Joint Equilibrium {{{1

# ** JEQ Auxiliary {{{2
# include("_JEQ/_jeq_auxiliary.jl")
function get_otc_values(firm_type::String, fidf::DataFrame,
                        yvar::Symbol; iota::Float64=NaN, kappa::Float64=NaN)


    # Transaction Costs and Volatility Risk Parameters
    cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [fidf[1, :sigmal]],
                                            :m  => [parse(Float64, string(fidf[1, :m]))],
                                            :gross_delta => [parse(Float64,
                                                                   string(fidf[1, :gross_delta]))],
                                            :mu_b => [1.0],
                                            :xi => [parse(Float64, string(fidf[1, :xi]))])

    xvar = :iota
    if !isnan(iota)
        cvmdict[:kappa] = [fidf[1, :kappa] * 1e-4]
    elseif !isnan(kappa)
        cvmdict[:iota] = [fidf[1, :iota] * 1e-4]
        xvar = :kappa
    else
        println("Please enter a value for iota or for kappa. Exiting...")
        return
    end

    svmdict = deepcopy(cvmdict)
    svmdict[:lambda] = [fidf[2, :lambda]]
    svmdict[:iota] = [.0]
    svmdict[:sigmah] = [fidf[2, :sigmah]]
    # #########################################################

    # Get Safe and Risky Firms' Full Info Optimal Results #####
    firm_obj_fun = :firm_value
    cvmdf, svmdf, _ = get_cvm_svm_dfs(cvmdict, svmdict;
                                      firm_obj_fun=firm_obj_fun)
    # #########################################################

    df = (firm_type == "safe") ? cvmdf : svmdf
    yinterp = Dierckx.Spline1D(df[xvar], df[yvar], k=3, bc="extrapolate")
    xval = (isnan(iota)) ? kappa : iota
    return yinterp(xval)
end


function get_otc_cut_off(y_otc::Float64,
                         yfun::Spline1D,
                         xgrid::StepRangeLen{Float64,
                                             Base.TwicePrecision{Float64},
                                             Base.TwicePrecision{Float64}})

    otc_diffs = yfun(xgrid) .- y_otc

    if all(otc_diffs .> .0)
        return Inf
    elseif all(otc_diffs .< .0)
        return -Inf
    else
        return (minimum(abs.(otc_diffs)) < 1e-4) ?  xgrid[argmin(abs.(otc_diffs))] : NaN
    end
end


function get_otc_cut_off_values(y_otc::Float64,
                                sep_yinterp::Spline1D, pool_yinterp::Spline1D,
                                xgrid::StepRangeLen{Float64,
                                             Base.TwicePrecision{Float64},
                                             Base.TwicePrecision{Float64}};
                                sep_otc::Float64=NaN,
                                pool_otc::Float64=NaN)

    if !isnan(y_otc)
        sep_otc = get_otc_cut_off(y_otc, sep_yinterp, xgrid)
        pool_otc = get_otc_cut_off(y_otc, pool_yinterp, xgrid)
    end

    if any([.&(isinf(sep_otc), sep_otc > 0), .&(isinf(pool_otc), pool_otc > 0)])
        return NaN
    elseif any(isnan.([sep_otc, pool_otc]) .== false)
        return maximum([x for x in [sep_otc, pool_otc] if !isnan(x)])
    else
        return NaN
    end
end


# ** JEQ Curves {{{2
# include("_JEQ/_jeq_curves.jl")
function plot_fi_curve(ax, fi_val::Float64; market::String="EP", kappa_val::Float64=NaN)
    # yvar::Symbol, fidf::DataFrame)

    if !(market in ["EP", "OTC"])
        println("Please enter 'EP' or 'OTC' for market type. Exiting... ")
        return
    end

    # Curve Style
    fi_curve_color = (market == "EP") ? fi_ep_curve_color : fi_otc_curve_color
    fi_linewidth = .8
    fi_linestyle = (market == "EP") ? jeq_linestyles[1] : "--"

    # Plot Curve
    ax.axhline(fi_val,
               color=fi_curve_color,
               linewidth=fi_linewidth,
               linestyle=fi_linestyle) #jeq_linestyles[1])

    # Set Label
    if !isnan(kappa_val)
        xlabel = (market == "EP") ? :kappa_ep : :kappa_otc
        fi_label = latexstring("\$", cvs_xlabels[xlabel][1], "=",
                               str_format_fun(cvs_xlabels[xlabel][2], kappa_val), "\$ (b.p.)")

    else
        fi_label = market
    end
    ax.text(.9, fi_val, fi_label, fontsize=10, va="bottom")

    return ax
end


function plot_pool_curve(ax, firm_type::String,
                         xvar::Symbol, yvar::Symbol,
                         pooldf::DataFrame;
                         xgrid::StepRangeLen{Float64,
                                             Base.TwicePrecision{Float64},
                                             Base.TwicePrecision{Float64}}=range(1., stop=1., length=10^5),
                         spline_k::Int64=3,
                         spline_bc::String="extrapolate")


    if !(firm_type in ["safe", "risky"])
        println("Please enter 'safe' or 'risky' for firm_type. Exiting...")
        return
    end
    fsym = (firm_type == "safe") ? :s_ : :r_


    pool_yinterp = Dierckx.Spline1D(pooldf[xvar], pooldf[Symbol(fsym, yvar)];
                                    k=spline_k, bc=spline_bc)

    ax.plot(pooldf[:mu_s], pooldf[Symbol(fsym, yvar)];
            color=pool_curve_color,
            linewidth=1,
            linestyle=jeq_linestyles[2],
            marker=jeq_markers[2],
            markersize=3)


    pool_label = "Pooling"
    xloc = .8 * maximum(pooldf[xvar])
    yloc =  pool_yinterp(xloc) - .6 * (pool_yinterp(1.05 * xloc) - pool_yinterp(.95 * xloc))
    ax.text(xloc, yloc, pool_label, fontsize=10, va="bottom")

    return pool_yinterp, ax
end


function plot_sep_curve(ax, firm_type::String,
                        xvar::Symbol, yvar::Symbol,
                        sepdf::DataFrame;
                        fi_val::Float64=NaN,
                        xgrid::StepRangeLen{Float64,
                                             Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}}=range(1., stop=1., length=10^5),
                        interp_yvar::Bool=true,
                        spline_k::Int64=3,
                        spline_bc::String="extrapolate")

    if !(firm_type in ["safe", "risky"])
        println("Please enter 'safe' or 'risky' for firm_type. Exiting...")
        return
    end
    fsym = (firm_type == "safe") ? :s_ : :r_


    if size(sepdf, 1) > 1
        # if !isnan(fi_val)
        #     sep_yinterp = Dierckx.Spline1D(vcat(sepdf[xvar], 1.),
        #                                    vcat(sepdf[Symbol(:s_, yvar)], fi_val);
        #                                    k=3, bc="extrapolate")
        # else
        #     sep_yinterp = Dierckx.Spline1D(sepdf[xvar], sepdf[Symbol(:s_, yvar)];
        #                                    k=3, bc="extrapolate")
        # end
        sep_yinterp = Dierckx.Spline1D(sepdf[xvar], sepdf[Symbol(fsym, yvar)];
                                       k=spline_k, bc=spline_bc)


        if interp_yvar
            ax.plot(xgrid[xgrid .> 0], sep_yinterp(xgrid[xgrid .> 0]),
                    color=sep_curve_color,
                    linewidth=1,
                    linestyle=jeq_linestyles[3])
        else
            ax.plot(sepdf[xvar], sepdf[Symbol(fsym, yvar)],
                    color=sep_curve_color,
                    linewidth=1,
                    linestyle=jeq_linestyles[3],
                    marker=jeq_markers[3],
                    markersize=3)
        end
    else
        sep_yinterp = Dierckx.Spline1D(xgrid, fill(sepdf[1, Symbol(fsym, yvar)], size(xgrid, 1));
                                       k=spline_k, bc=spline_bc)

        ax.axhline(sepdf[1, Symbol(fsym, yvar)],
                   color=sep_curve_color,
                   linewidth=1,
                   linestyle=jeq_linestyles[3])
    end

    sep_label = "Separating"
    xloc = .5 * maximum(xgrid)
    yloc =  sep_yinterp(xloc) #+ 1.5 * (sep_yinterp(1.05 * xloc) - sep_yinterp(.95 * xloc))
    ax.text(xloc, yloc, sep_label, fontsize=10, va="top")

    return sep_yinterp, ax
end


function jeq_plot_vlines(ax, xvar::Symbol;
                         fv_xvar::Float64=NaN,
                         fv_color::String=fv_color,
                         mbr_xvar::Float64=NaN,
                         mbr_color::String=mbr_color)

    xval = (!isnan(fv_xvar)) ? fv_xvar : mbr_xvar
    ax.axvline(xval, color="black", linewidth=.6, linestyle="--")

    # Form Dictionary with labels and values:
    vldict = vlines_labels_dict(xvar; fv_xvar=fv_xvar,
                                fv_color=fv_color,
                                mbr_xvar=mbr_xvar,
                                mbr_color=mbr_color)

    xkeys = [x for x in keys(vldict) if .&(!isnan(vldict[x][:value]), !isinf(vldict[x][:value]))]
    minor_ticks = [vldict[x][:value] for x in xkeys]
    minor_labels = [vldict[x][:xsym] for x in xkeys]
    for x in xkeys
        ax.axvline(vldict[x][:value],
                   color=vldict[x][:color],
                   linewidth=.6,
                   linestyle="--",
                   marker=svmmarkers[1])
    end

    ax.set_xticks(minor_ticks, minor=true)
    ax.set_xticklabels(minor_labels, minor=true)

    return ax
end


# ** JEQ Color Regions {{{2
# include("_JEQ/_jeq_color_regions.jl")
function color_otc_region_fun(ax, fv_xvar::Float64,
                              xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                  Base.TwicePrecision{Float64}},
                              text_xloc::Float64,
                              text_yloc::Float64)
    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .<= fv_xvar,
                    facecolor=otc_region_color, alpha=0.25)

    if !isinf(text_xloc)
        ax.text(text_xloc, text_yloc,
                "OTC Trading \n is optimal",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=11,
                style="italic",
                bbox=Dict("facecolor" => box_color, "alpha" => 0.5, "pad" => 10))
    end

    return ax
end


# ** JEQ Plot Methods {{{2
# include("_JEQ/_jeq_plot_methods.jl")
function jeq_subplotfun(fig::Figure, firm_type::String,
                        xvar::Symbol, yvar::Symbol,
                        fidf::DataFrame,
                        misrepdf::DataFrame,
                        pooldf::DataFrame,
                        sepdf::DataFrame,
                        xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}};
                        kappa_otc_bp::Float64=NaN,
                        ax_subplot::Int64=111,
                        interp_yvar::Bool=false,
                        fv_xvar::Float64=NaN,
                        mbr_xvar::Float64=NaN,
                        misrep_xvar::Float64=NaN,
                        color_rm_region::Bool=true,
                        color_nrm_region::Bool=true,
                        color_conflict_region::Bool=false,
                        color_misrep_region::Bool=false,
                        cvmlinestyles::Array{String,1}=cvmlinestyles,
                        cvmmarkers::Array{String,1}=cvmmarkers,
                        svmlinestyles::Array{String,1}=svmlinestyles,
                        svmmarkers::Array{String,1}=svmmarkers)

    if !(firm_type in ["safe", "risky"])
        println("Please enter 'safe' or 'risky' for firm_type. Exiting...")
        return
    end

    ax = fig.add_subplot(ax_subplot)

    frow = (firm_type == "safe") ? 1 : 2
    # Plot Full Information ################################################
    # Electronic Market
    ax = plot_fi_curve(ax, fidf[frow, yvar]; market="EP", kappa_val=fidf[frow, :kappa])

    y_otc = NaN
    if !isnan(kappa_otc_bp)
        y_otc = get_otc_values(firm_type, fidf, yvar; kappa=kappa_otc_bp * 1e-4)
        ax = plot_fi_curve(ax, y_otc; market="OTC", kappa_val=kappa_otc_bp)
    end


    # Separating Equilibrium ###############################################
    sep_yinterp, ax = plot_sep_curve(ax, firm_type,
                                     xvar, yvar, sepdf;
                                     fi_val=fidf[frow, yvar], xgrid=xgrid)


    # Pooling Equilibrium ###################################################
    pool_yinterp, ax = plot_pool_curve(ax, firm_type,
                                       xvar, yvar, pooldf, xgrid=xgrid)


    #
    # sep_otc = NaN
    # pool_otc = NaN
    # if !isnan(y_otc)
    #     sep_otc = get_otc_cut_off_value(y_otc, sep_yinterp, xgrid)
    #     pool_otc = get_otc_cut_off_value(y_otc, pool_yinterp, xgrid)
    # end

    mu_star = get_otc_cut_off_values(y_otc,
                                     sep_yinterp, pool_yinterp,
                                     xgrid)
    fv_xvar = NaN
    mbr_xvar = NaN
    # if any(isnan.([sep_otc, pool_otc]) .== false)
    #     mu_star = minimum([x for x in [sep_otc, pool_otc] if !isnan(x)])
    if !isnan(mu_star)
        if yvar == :firm_value
            fv_xvar = mu_star
        elseif yvar == :MBR
            mbr_xvar = mu_star
        end

        ax = jeq_plot_vlines(ax, xvar; fv_xvar=fv_xvar, mbr_xvar=mbr_xvar)
    end


    # Axes Limits ##################################################
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # ##############################################################


    if !isnan(fv_xvar)
        xloc = fv_xvar/2
        yloc = (.65 * ymax + .35 * ymin)
        ax = color_otc_region_fun(ax, fv_xvar, xgrid, xloc, yloc)
    end


        # sep_otc_diff = minimum(abs.(sep_yinterp(xgrid) .- y_otc))
        # pool_otc_diff = minimum(abs.(pool_yinterp(xgrid) .- y_otc))

        # mu_star = NaN
        # if .&(sep_otc_diff < 1e-4, pool_otc_diff < 1e-4)
        #     sep_otc = xgrid[argmin(abs.(sep_yinterp(xgrid) .- y_otc))]
        #     pool_otc = xgrid[argmin(abs.(pool_yinterp(xgrid) .- y_otc))]
        #     mu_star = minimum([sep_otc, pool_otc])
        # elseif sep_otc_diff < 1e-4
        #     mu_star = xgrid[argmin(abs.(sep_yinterp(xgrid) .- y_otc))]
        # elseif pool_otc_diff < 1e-4
        #     mu_star = xgrid[argmin(abs.(pool_yinterp(xgrid) .- y_otc))]
        # end

        # ax.axvline(mu_star)



    # For some reason, y limits are not matching axis ylim.
    # Force reset:
    ax.set_xlim([minimum(xgrid), maximum(xgrid)])
    ax.set_ylim([ymin, ymax])

    return ax
end


function jeq_core_plot(fig, firm_type::String,
                       yvars::Array{Symbol,1},
                      fidf::DataFrame,
                      misrepdf::DataFrame,
                      pooldf::DataFrame,
                      sepdf::DataFrame,
                      xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                              Base.TwicePrecision{Float64}},
                      subplots::Int64;
                      kappa_otc_bp::Float64=NaN,
                       interp_yvar::Bool=false,
                       fv_xvar::Float64=NaN,
                       mbr_xvar::Float64=NaN,
                       misrep_xvar::Float64=NaN,
                       cvmlinestyles::Array{String,1}=cvmlinestyles,
                       cvmmarkers::Array{String,1}=cvmmarkers,
                       svmlinestyles::Array{String,1}=svmlinestyles,
                       svmmarkers::Array{String,1}=svmmarkers,
                       color_rm_region::Bool=true,
                       color_nrm_region::Bool=true,
                       color_conflict_region::Bool=false,
                       color_misrep_region::Bool=false)




    xvar = :mu_s
    xvar_xlabel = latexstring("Measure of safe firms \$", jeq_xlabels[:mu_s][1], "\$")

    axes = []
    count = 1
    for yvar in yvars
        ax_subplot = ax_subplots[subplots][count]

        ax = jeq_subplotfun(fig, firm_type,
                            xvar, yvar,
                            fidf, misrepdf,
                            pooldf, sepdf,
                            xgrid;
                            ax_subplot=ax_subplot,
                            kappa_otc_bp=kappa_otc_bp,
                            interp_yvar=interp_yvar,
                            fv_xvar=fv_xvar,
                            mbr_xvar=mbr_xvar,
                            misrep_xvar=misrep_xvar,
                            cvmlinestyles=cvmlinestyles,
                            cvmmarkers=cvmmarkers,
                            svmlinestyles=svmlinestyles,
                            svmmarkers=svmmarkers,
                            color_rm_region=color_rm_region,
                            color_nrm_region=color_nrm_region,
                            color_conflict_region=color_conflict_region,
                            color_misrep_region=color_misrep_region)


        # ##############################################################
        # ################### Plot Labels and Titles ###################
        # ##############################################################
        if subplots == 1
            ax.set_xlabel(xvar_xlabel, labelpad=10)
            ax.set_title(" ")
        elseif count == 1
            ax.set_title(" \n ")
        end
        ax.set_ylabel(cvs_ylabels[yvar], labelpad=10)
        # ##############################################################

#         # ##############################################################
#         # ########################## Safe Type #########################
#         # ##############################################################
#         if !isempty(misrepdf)
#             misrep_xval = misrepdf[1, Symbol(:s_, xvar)]
#             misrep_yval = misrepdf[1, Symbol(:s_, yvar)]
#             ax.scatter(misrep_xval,
#                        misrep_yval;
#                        s=25,color="purple", alpha=0.8)
#             ax.text(misrep_xval, 1.0025 * misrep_yval, "\$S\$")
#         end
#         # ##############################################################

        push!(axes, ax)
        count += 1
    end

    if subplots > 1
        PyPlot.matplotlib.pyplot.xlabel(xvar_xlabel,
                                        labelpad=10)
    end


    return fig
end


function jeq_plotfun(firm_type::String, yvars::Array{Symbol,1},
                     fidf::DataFrame,
                     misrepdf::DataFrame,
                     pooldf::DataFrame,
                     sepdf::DataFrame;
                     xgrid::StepRangeLen{Float64,
                                         Base.TwicePrecision{Float64},
                                         Base.TwicePrecision{Float64}}=range(.0, stop=0., length=0),
                     interp_yvar::Bool=false,
                     kappa_otc_bp::Float64=NaN,
                     fv_xvar::Float64=NaN,
                     mbr_xvar::Float64=NaN,
                     misrep_xvar::Float64=NaN,
                     color_rm_region::Bool=true,
                     color_nrm_region::Bool=true,
                     color_conflict_region::Bool=false,
                     color_misrep_region::Bool=false,
                     cvmlinestyles::Array{String,1}=cvmlinestyles,
                     cvmmarkers::Array{String,1}=cvmmarkers,
                     svmlinestyles::Array{String,1}=svmlinestyles,
                     svmmarkers::Array{String,1}=svmmarkers,
                     figaspect::Float64=NaN,
                     figsize::Tuple{Float64, Float64}=(.0, .0),
                     figpad::Float64=1.8,
                     save_fig::Bool=true,
                     fig_dpi::Int64=300,
                     main_dir_path::String=main_dir_path,
                     plots_dir::String=plots_dir,
                     rmp_plots_dir::String=rmp_plots_dir)


    if !(firm_type in ["safe", "risky"])
        println("Please enter 'safe' or 'risky' for firm_type. Exiting...")
        return
    end

    # Figure Size and Layout Aspect
    subplots = size(yvars, 1)
    figsize, figpad = ModelPlots.fig_size_pad_adjuster(subplots;
                                            figaspect=figaspect,
                                            figsize=figsize,
                                            figpad=figpad)


    if size(xgrid, 1) == 0
        xgrid = range(.0, stop=1., length=10^5)
    end

    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=figsize)
    fig = jeq_core_plot(fig, firm_type, yvars,
                        fidf, misrepdf,
                        pooldf, sepdf,
                        xgrid, subplots;
                        kappa_otc_bp=kappa_otc_bp,
                        interp_yvar=interp_yvar,
                        fv_xvar=fv_xvar,
                        mbr_xvar=mbr_xvar,
                        misrep_xvar=misrep_xvar,
                        cvmlinestyles=cvmlinestyles,
                        cvmmarkers=cvmmarkers,
                        svmlinestyles=svmlinestyles,
                        svmmarkers=svmmarkers,
                        color_rm_region=color_rm_region,
                        color_nrm_region=color_nrm_region,
                        color_conflict_region=color_conflict_region,
                        color_misrep_region=color_misrep_region)


    # Set Sup Title
    suptitle_yvars = join([cvs_ylabels[yvar] for yvar in yvars], " and ")
    suptitle_params1 = latexstring("\$\\overline{", cvs_xlabels[:iota][1], "}=\$",
                                   str_format_fun(cvs_xlabels[:iota][2],
                                                  parse(Float64, string(fidf[1, :iota]))),
                                   " (b.p.)" )

    suptitle_params2 = join([string("\$", tlabels[x][1], "= \$ ",
                                    str_format_fun(tlabels[x][2], parse(Float64, string(sepdf[1, x]))))
                             for x in jeq_plots_title_params_order], ", ")
    plot_suptitle = latexstring("Safe Type's ", suptitle_yvars, "\n",
                                " for ", suptitle_params1, ", ", suptitle_params2)
    fig.suptitle(plot_suptitle, fontsize=14)


    if !isnan(figpad)
        fig.tight_layout(pad=figpad)
    end

    return fig
end


# * Contour Plots {{{1

# ** Contour Auxiliary {{{2
# include("_Contour/_contour_auxiliary.jl")

function get_eq_type_df(eq_type::String,
                        fidf::DataFrame, misrepdf::DataFrame,
                        sepdf::DataFrame, pooldf::DataFrame)
    if eq_type == "full_info"
        return fidf
    elseif eq_type == "misrep"
        return misrepdf
    elseif eq_type == "separating"
        return sepdf
    elseif eq_type == "pooling"
        return pooldf
    end
end


function slice_df(df::DataFrame, svar::Symbol; tol::Float64=1e-5)
    if any(isnan.(df[:, svar]))
        return df[isnan.(df[:, svar]) .==false, :]
    elseif any(abs.(df[:, svar] .- .0) .< tol)
        return df[abs.(df[:, svar] .- .0) .> tol, :]
    end
end


function interp_z_values(df::DataFrame;
                         xvar::Symbol=contour_xvar,
                         yvar::Symbol=contour_yvar,
                         zvars::Array{Symbol, 1}=contour_zvars,
                         ft_xy::Symbol=Symbol(""),
                         ft_z::Array{Symbol, 1}=[:s_, :r_],
                         spline_k::Int64=3,
                         spline_bc::String="extrapolate")
  if !(xvar in Symbol.(names(df)))
        ft_xy = :r_
    end

    # Separate DataFrames
    xdf = slice_df(df, Symbol(ft_xy, xvar))
    ydf = slice_df(df, Symbol(ft_xy, yvar))


    # Form Dictionary to store results
    tmpdict = Dict{Symbol, Any}(zip([xvar, yvar],
                                    [Spline1D(1:5, 1:5), Spline1D(1:5, 1:5)]))
    # fd = Dict(zip(zvars, repeat([deepcopy(tmpdict)], 1, size(zvars, 1))))

    println(Symbol(ft_xy, yvar))
    fd = Dict{Symbol, Any}(:xvar => xvar,
                           :yvar => yvar,
                           :xvals => xdf[:, Symbol(ft_xy, xvar)],
                           :yvals => ydf[:, Symbol(ft_xy, yvar)])

    if !(zvars[1] in Symbol.(names(df)))
        zvars = vcat([Symbol(prefix_z, zvar) for prefix_z in ft_z, zvar in zvars]...)
    end

    # Interpolate Functions
    for zvar in zvars
        fd[zvar] = deepcopy(tmpdict)
        fd[zvar][xvar] = Dierckx.Spline1D(fd[:xvals], xdf[:, zvar],
                                          k=spline_k, bc=spline_bc)
        fd[zvar][yvar] = Dierckx.Spline1D(fd[:yvals], ydf[:, zvar],
                                          k=spline_k, bc=spline_bc)
    end

    return fd
end


function form_mesh_grid(xvals::Array{Float64,1},
                        yvals::Array{Float64,1},
                        zfun; N::Int64=200)
    xgrid = range(minimum(xvals), stop=maximum(xvals), length=N)
    ygrid = range(minimum(yvals), stop=maximum(yvals), length=N)

    X = Array(repeat(xgrid, 1, N)')
    Y = repeat(ygrid, 1, N)
    Z = Array([zfun(x,y) for x in xgrid, y in ygrid]')

    return X, Y, Z
end


function get_contour_plot_path_name(df::DataFrame, zfun_name::Symbol;
                                    firm_type::Symbol=Symbol(""), fname_eq_type::String="",
                                    fname_ext::String=contour_fname_ext,
                                    nrmp::Bool=false)
    eq_type = ""
    if isempty(fname_eq_type)
        eq_type = String(df[1, :eq_type])

        if eq_type == "pool"
            eq_type = "pooling"
        elseif eq_type == "sep"
            eq_type = "separating"
        end

        fname_eq_type = eq_type_title[eq_type][1]
    end

    if .&(eq_type != "full_info", any([isempty(string(firm_type)),
                                      !(firm_type in [:safe, :risky])]))
        println("Please enter a firm type: :safe or :risky. Exiting...")
        return
    end


    ft = ""
    iota_s = NaN
    lambda_r = NaN
    mu_s = NaN
    if !nrmp
        if eq_type == "full_info"
            iota_s = minimum([x for x in df[:, :rm_iota] if x > 0.])
            lambda_r = unique([x for x in df[:, :nrm_lambda] if !isnan(x)])[1]
        else
            ft = (firm_type == :safe) ? :s_ : :r_
            mu_s = df[1, :mu_s]
            iota_s = df[1, :s_iota]
            lambda_r = unique([x for x in df[:, :r_lambda] if !isnan(x)])[1]
        end
    end

    # Inverse Coupon Rate
    pcr = df[1, :p]/df[1, :c]

    if !nrmp
        fname = string(fname_eq_type, "_mu_s_", mu_s, "_pcr_", pcr,
                       "_iota_s_", iota_s, "__lambda_r_", lambda_r,
                       "__", ft,  zfun_name)
    else
        fname = string("nrmp_", fname_eq_type, "_mu_s_", mu_s, "_pcr_", pcr,
                       "__", ft,  zfun_name)
    end

    return string(contour_plots_path, "/", fname, ".", fname_ext)
end


function get_region_grids(xgrid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},                          ygrid::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},
               eqfun; N::Int64=10^5)
    ymin = minimum(ygrid)
    ymax = maximum(ygrid)

    yvals = x -> [y for y in ygrid if eqfun(x, y)]
    min_y = x -> !isempty(yvals(x)) ? minimum(yvals(x)) : .0
    max_y = x -> !isempty(yvals(x)) ? maximum(yvals(x)) : .0

    min_y_grid = fetch(@spawn [minimum([maximum([ymin, min_y(x)]), ymax]) for x in xgrid])
    max_y_grid = fetch(@spawn [minimum([maximum([ymin, max_y(x)]), ymax]) for x in xgrid])
    min_y_fun = Dierckx.Spline1D(xgrid, min_y_grid; k=3, bc="extrapolate")
    max_y_fun = Dierckx.Spline1D(xgrid, max_y_grid; k=3, bc="extrapolate")

    ref_xgrid = range(minimum(xgrid), stop=maximum(xgrid), length=N)
    min_y_ref_grid = [minimum([maximum([ymin, min_y_fun(x)]), ymax]) for x in ref_xgrid]
    max_y_ref_grid = [minimum([maximum([ymin, max_y_fun(x)]), ymax]) for x in ref_xgrid]

    # return min_y, max_y, ref_xgrid, min_y_ref_grid, max_y_ref_grid
    return ref_xgrid, min_y_ref_grid, max_y_ref_grid
end



function get_eq_contour_mesh_grid(xvals::Array{Float64,1}, yvals::Array{Float64,1},
                                  fun_dict; N::Int64=10^3)

    # eq_bool = (x, y) -> fun_dict[:fi_ind](x, y) + 2 * fun_dict[:sep_ind](x, y) + 3 * fun_dict[:pool_ind](x, y)
    # eq_vals = (x, y) -> (fun_dict[:fi_ind](x,y) * fun_dict[:mbr][:fi](x, y) +
    #                      fun_dict[:sep_ind](x,y) * fun_dict[:mbr][:sep](x, y) +
    #                      fun_dict[:pool_ind](x,y) * fun_dict[:mbr][:pool](x, y))

    X, Y, bool_Z = form_mesh_grid(xvals, yvals, fun_dict[:eq_bool], N=N)
    _, _, bool_OTC_EP = form_mesh_grid(xvals, yvals, fun_dict[:bool_otc_ep], N=N)
    _, _, r_mbr = fetch(@spawn form_mesh_grid(xvals, yvals, fun_dict[:r_mbr], N=N))
    _, _, s_fv = fetch(@spawn form_mesh_grid(xvals, yvals, fun_dict[:s_fv], N=N))

    return Dict{Symbol, Any}(:X => X, :Y => Y,
                             :bool_Z => bool_Z,
                             :bool_OTC_EP => bool_OTC_EP,
                             :r_mbr => r_mbr,
                             :s_fv => s_fv)
end


# ** Contour Plot Methods {{{2
# include("_Contour/_contour_plot_methods.jl")
function get_title_value(df, var)
  if !(var in Symbol.(names(df)))
        var = Symbol(:r_, var)
    end

    return [parse(Float64, string(x)) for x in unique(df[:, var]) if !isnan(x)][1]
end


function get_final_contour_plot_title(df::DataFrame, zvar::Symbol,
                                      ft::Symbol;
                                      k_otc::Float64=NaN,
                                      params_list::Array{Symbol,1}=vcat(:mu_s,
                                                                        contour_plots_title_params_order))

    df[!, :pcr] .= df[1, :p]/df[1, :c]
    title_params = join([string("\$", contour_tlabels[x][1], "= \$ ",
                                str_format_fun(contour_tlabels[x][2],
                                               get_title_value(df, x)))
                         for x in params_list], ", ")

    firm_type_title = ""
    if ft in [:st, :safe]
        firm_type_title = "Safe Type's "
    elseif ft in [:rt, :risky]
        firm_type_title = "Risky Type's "
    else
        println("Error! Firm type not recognized. Exiting...")
        return
    end

    k_otc_title = ""
    if isnan(k_otc)
        eq_type_title = " in the Prevailing EP Market Equilibria "
    else
        eq_type_title = " in the Prevailing Dual Market Equilibria "
        k_otc_title = string(", \$", contour_tlabels[:kappa_otc][1], "= \$ ",
                             str_format_fun(contour_tlabels[:kappa_otc][2], k_otc))
    end

    plot_title = latexstring(firm_type_title, "Optimal ",
                             contour_tlabels[zvar][1],
                             eq_type_title,
                             "\n for ", title_params, k_otc_title)

    return plot_title
end


function get_contour_plot_title(df::DataFrame,
                                # eqfuns::Dict{Symbol, Any},
                                zvar::Symbol;
                                nrmp::Bool=false,
                                ft::Symbol=Symbol(""),
                                diff_fun::Bool=false,
                                rm_prefix::Symbol=Symbol(""),
                                params_list::Array{Symbol,1}=contour_plots_title_params_order,
                                s_fi_fv::Float64=NaN)

    eq_type = String(df[1, :eq_type])

    if nrmp
        params_list = [x for x in params_list if !(x == :lambda)]
    end

    if eq_type == "pool"
        eq_type = "pooling"
    elseif eq_type == "sep"
        eq_type = "separating"
    end


    mu_s = NaN
    if eq_type in ["pooling", "separating"]
        mu_s =  df[1, :mu_s]
        params_list = vcat(:mu_s, params_list)
    end
    df[!, :pcr] .= df[1, :p]/df[1, :c]

    title_params = join([string("\$", contour_tlabels[x][1], "= \$ ",
                                str_format_fun(contour_tlabels[x][2],
                                               get_title_value(df, x)))
                         for x in params_list], ", ")


    firm_type_title = ""
    if ft == :safe
        firm_type_title = "Safe Type's "
    else
        firm_type_title = "Risky Type's "
    end

    if diff_fun
        plot_title = latexstring(firm_type_title, eq_type_title[eq_type][2],
                                 " v.s. Full Information Eq. ",
                                 contour_tlabels[zvar][1], " Differential ",
                                 "\n for ", title_params)
    else
        if !isempty(string(rm_prefix))
            rm_policy = (rm_prefix == :rm) ? " Risk Management " : "No Risk Management "
            firm_type_title = string(firm_type_title, rm_policy)
        end

        if eq_type == "misrep"
            plot_title = latexstring(firm_type_title,
                                     contour_tlabels[zvar][1],
                                     " in case of Misrepresentation",
                                     "\n for ", title_params)
        else
            plot_title = latexstring(firm_type_title, "Optimal ",
                                     contour_tlabels[zvar][1], " in a ",
                                     eq_type_title[eq_type][2], " Equilibrium ",
                                     "\n for ", title_params)
        end
    end

    if !isnan(s_fi_fv)
        plot_title = latexstring(plot_title,
                                 "\n (Safe Type's Full Information Firm Value = ",
                                 str_format_fun("%.2f", s_fi_fv), ") \n")
    end

    return plot_title
end


function plot_iso_curves(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2};
                         iso_xlabel::Symbol=:iota,
                         iso_ylabel::Symbol=:sigmah,
                         rmp_Z::Array{Float64, 2}=Array{Float64,2}(undef, 0, 2),
                         seaborn_style=iso_plt_inputs[:seaborn_style],
                         iso_levels=iso_plt_inputs[:iso_levels],
                         heat_levels=iso_plt_inputs[:heat_levels],
                         iso_cmap=iso_cmaps["full_info"],
                         heat_cmap::String="",
                         fig_aspect=iso_plt_inputs[:fig_aspect],
                         iso_fontsize=iso_plt_inputs[:iso_fontsize],
                         use_subgrid=iso_plt_inputs[:use_subgrid],
                         subgrid_rows=iso_plt_inputs[:subgrid_rows],
                         iso_cols=iso_plt_inputs[:iso_cols],
                         heat_cols=iso_plt_inputs[:heat_cols],
                         cat_Z=[],
                         cat_cmap="GnBu",
                         cat_alpha::Float64=.25,
                         rmp_cmap=rmp_cmaps["misrep"],
                         rmp_alpha::Float64=.15)
                         # cat_Z::Array{Int64, 2}=Array{Int64, 2}[])

    if isempty(heat_cmap)
        heat_cmap = iso_cmap
    end

    if !isempty(seaborn_style)
        Seaborn.set(style=seaborn_style)
    end

    w, h = figaspect(fig_aspect)
    fig = PyPlot.figure(figsize=(w, h))

    # Choose between subgrids or subplots ##################################
    if use_subgrid
        fig = PyPlot.figure(figsize=(w, h))
        ax1 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols),
                                  (0, 0), rowspan=subgrid_rows, colspan=iso_cols)
        ax2 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols),
                                  (0, iso_cols), rowspan=subgrid_rows, colspan=heat_cols)
    else
        fig, axs = PyPlot.subplots(1, 2, figsize=(w, h), sharey=true)
        ax1 = axs[1] # fig.add_subplot(121)
        ax2 = axs[2] # fig.add_subplot(122)
    end
    # ######################################################################

    CS = ax1.contour(X, Y, Z, levels=iso_levels, cmap=iso_cmap)
    ax1.clabel(CS, inline=5, fontsize=iso_fontsize)
    ax1.set_xlabel(latexstring("\$", xylabels[iso_xlabel][1], "\$"), labelpad=10)
    ax1.set_ylabel(latexstring("\$", xylabels[iso_ylabel][1], "\$"), labelpad=10)

    # Color Region where Firm's change their RMP relative to
    # their Full Information Eq. RMP:
    if !isempty(rmp_Z)
        CS3 = ax1.contourf(X, Y, rmp_Z, cmap=rmp_cmap) #, alpha=rmp_alpha)
    end

    CS2 = ax2.contourf(X, Y, Z, levels=heat_levels, cmap=heat_cmap)
    if use_subgrid
        ax2.tick_params(
            axis="y",          # changes apply to the x-axis
            which="both",      # both major and minor ticks are affected
            bottom=false,      # ticks along the bottom edge are off
            top=false,         # ticks along the top edge are off
            left=false,
            right=false,
            labelleft=false,
            labelbottom=false)
    end
    ax2.set_xlabel(latexstring("\$", xylabels[iso_xlabel][1], "\$"), labelpad=10)

    # Add Colorbar
    cb2 = fig.colorbar(CS2)

    if !isempty(cat_Z)
        cats = sort(unique(cat_Z))
        cat_tick_labels = [eq_cat_dict[x][2] for x in [:fi, :sep, :pool, :otc]
                           if eq_cat_dict[x][1] in cats]

        if size(cats, 1) < size([x for x in keys(eq_cat_dict)], 1)
            cat_Z = cat_Z .- 1
            cats = cats .- 1
        end

        cat_levels = size(cats, 1) - 1
        CS1 = ax1.contourf(X, Y, cat_Z,
                           cmap=cat_cmap, levels=cat_levels) #, alpha=cat_alpha)
        cb1 = fig.colorbar(CS1, ax=ax1, ticks=reverse(cats))#, orientation="horizontal")
        cb1.set_ticklabels(cat_tick_labels)
        cb1.set_clim(1, cat_levels + 1)
    end


    return fig, ax1, ax2
end


function plot_iso_contour_curves(fd::Dict{Symbol, Any},
                                 zfun;
                                 rmp_diff_fun=nothing,
                                 fig_title::LaTeXString=LaTeXString(""),
                                 file_path_name::String="",
                                 seaborn_style=iso_plt_inputs[:seaborn_style],
                                 iso_levels=iso_plt_inputs[:iso_levels],
                                 heat_levels=iso_plt_inputs[:heat_levels],
                                 iso_cmap=iso_cmaps["full_info"],
                                 heat_cmap::String="",
                                 fig_aspect=iso_plt_inputs[:fig_aspect],
                                 iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                 use_subgrid=iso_plt_inputs[:use_subgrid],
                                 subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                 iso_cols=iso_plt_inputs[:iso_cols],
                                 heat_cols=iso_plt_inputs[:heat_cols],
                                 title_font_size=iso_plt_inputs[:title_font_size],
                                 fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                 tight_pad=iso_plt_inputs[:tight_pad],
                                 h_pad=iso_plt_inputs[:h_pad],
                                 w_pad=iso_plt_inputs[:w_pad],
                                 rmp_cmap=rmp_cmaps["misrep"],
                                 rmp_alpha::Float64=.15)

    X, Y, Z = form_mesh_grid(fd[:xvals], fd[:yvals], zfun)

    rmp_Z = Array{Float64, 2}(undef, 0, 2)
    if rmp_diff_fun != nothing
        _, _, rmp_Z = ModelPlots.form_mesh_grid(fd[:xvals], fd[:yvals], rmp_diff_fun)
    end

    fig, ax1, ax2 = plot_iso_curves(X, Y, Z;
                                    rmp_Z=rmp_Z,
                                    seaborn_style=seaborn_style,
                                    iso_levels=iso_levels,
                                    heat_levels=heat_levels,
                                    iso_cmap=iso_cmap,
                                    heat_cmap=heat_cmap,
                                    fig_aspect=fig_aspect,
                                    iso_fontsize=iso_fontsize,
                                    use_subgrid=use_subgrid,
                                    subgrid_rows=subgrid_rows,
                                    iso_cols=iso_cols,
                                    heat_cols=heat_cols,
                                    rmp_cmap=rmp_cmap,
                                    rmp_alpha=rmp_alpha)

    if !isempty(fig_title)
        fig.suptitle(fig_title, fontsize=title_font_size)
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end


function plot_equilibria_iso_contour_curves(X, Y, Z, eq_type_Z;
                                            sup_title_sep::Bool=false,
                                            fig_title::LaTeXString=LaTeXString(""),
                                            file_path_name::String="",
                                            iso_xlabel::Symbol=:iota,
                                            iso_ylabel::Symbol=:sigmah,
                                            seaborn_style=iso_plt_inputs[:seaborn_style],
                                            iso_levels=iso_plt_inputs[:iso_levels],
                                            heat_levels=iso_plt_inputs[:heat_levels],
                                            iso_cmap=iso_cmaps["full_info"],
                                            heat_cmap::String="",
                                            fig_aspect=iso_plt_inputs[:fig_aspect],
                                            iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                            use_subgrid=iso_plt_inputs[:use_subgrid],
                                            subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                            iso_cols=iso_plt_inputs[:iso_cols],
                                            heat_cols=iso_plt_inputs[:heat_cols],
                                            title_font_size=iso_plt_inputs[:title_font_size],
                                            fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                            tight_pad=iso_plt_inputs[:tight_pad],
                                            h_pad=iso_plt_inputs[:h_pad],
                                            w_pad=iso_plt_inputs[:w_pad],
                                            cat_cmap="GnBu",
                                            rmp_cmap=rmp_cmaps["misrep"])

    fig, ax1, ax2 = ModelPlots.plot_iso_curves(X, Y, Z;
                                               iso_xlabel=iso_xlabel,
                                               iso_ylabel=iso_ylabel,
                                               iso_cmap=iso_cmap,
                                               iso_levels=15,
                                               cat_Z=eq_type_Z,
                                               cat_cmap=cat_cmap,
                                               rmp_cmap=rmp_cmap)
#    ax1.contourf(X, Y, eq_type_Z, cmap="GnBu_r", alpha=.4)

    if !isempty(fig_title)
        ttl = fig.suptitle(fig_title, fontsize=title_font_size)

        if sup_title_sep
            ttl.set_position([.5, 1.05])
        end
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end


# ** Contour Payoff Functions {{{2
# include("_Contour/_contour_payoff_functions.jl")
function check_rm_vars(xvar::Symbol, yvar::Symbol)
    if .&(xvar != :iota, yvar != :iota)
        println("Either xvar or yvar must be :iota. Exiting...")
        return false
    end

    return true
end


function get_rm_payoff_funs(fd::Dict, xvar::Symbol, yvar::Symbol, obj_fun::Symbol)
    rm_vars_present = check_rm_vars(xvar, yvar)
    if !rm_vars_present
        return
    end

    if xvar == :iota
        rm_fun = (x, y) -> fd[obj_fun][xvar](x)
        nrm_fun = (x, y) -> fd[obj_fun][yvar](y)
    else
        nrm_fun = (x, y) -> fd[obj_fun][xvar](x)
        rm_fun = (x, y) -> fd[obj_fun][yvar](y)
    end

    return rm_fun, nrm_fun
end



function fi_payoff_functions(fi_fd::Dict;
                             xvar::Symbol=:iota, yvar::Symbol=:sigmah)

    fi_funs = Dict{Symbol, Any}(:xvar => xvar,
                                :yvar => yvar)

    # Firm Value ##################################################
    # Get RMP-Specific Payoffs
    fi_funs[:rm_fv], fi_funs[:nrm_fv] = get_rm_payoff_funs(fi_fd, xvar, yvar, :firm_value)


    # RM Payoff > NRM Payoff ?
    fi_funs[:rm_cond] = (x, y) -> fi_funs[:rm_fv](x, y) >= fi_funs[:nrm_fv](x, y)


    # Payoff is the Maximum RMP-Conditional Firm Value
    fi_funs[:fv] = (x, y) -> maximum([fi_fd[:firm_value][xvar](x),
                                      fi_fd[:firm_value][yvar](y)])
    # #############################################################

    for zvar in [z for z in contour_zvars if z != :firm_value]
        zsym = contour_zvars_sym[zvar]
        fi_funs[Symbol(:rm_, zsym)], fi_funs[Symbol(:nrm_, zsym)] = get_rm_payoff_funs(fi_fd, fi_fd[:xvar], fi_fd[:yvar], zvar)
        fi_funs[zsym] = (x, y) -> fi_funs[:rm_cond](x, y) ? fi_funs[Symbol(:rm_, zsym)](x, y) : fi_funs[Symbol(:nrm_, zsym)](x, y)
    end

    return fi_funs
end


function misrep_payoff_functions(fi_funs::Dict{Symbol, Any}, mp_fd::Dict{Symbol, Any};
                                 xvar::Symbol=:iota, yvar::Symbol=:sigmah)

    if any([fi_funs[:xvar] != xvar, fi_funs[:yvar] != yvar])
        println("Full Info and Misrep x and y variables do not coincide. Exiting...")
        return
    end

    mp_funs = Dict{Symbol, Any}(:xvar => xvar, :yvar => yvar)

    # Compute Type-Specific MBR under Full Information ################
    mp_funs[:rm_mbr], mp_funs[:nrm_mbr] = get_rm_payoff_funs(mp_fd, mp_fd[:xvar], mp_fd[:yvar], :r_MBR)
    mp_funs[:rm_cond] = (x, y) -> mp_funs[:rm_mbr](x, y) >= mp_funs[:nrm_mbr](x, y)

    # Choose RMP that maximizes the MBR in case of Misrepresentation
    mp_funs[:mbr] = (x, y) -> maximum([mp_fd[:r_MBR][xvar](x),
                                       mp_fd[:r_MBR][yvar](y)])

    # Difference between Misrep MBR and FI MBR
    mp_funs[:mp_fi_mbr_diff] = (x, y) -> mp_funs[:mbr](x, y) - fi_funs[:mbr](x, y)

    # Does the Misrep RMP differ from the optimal RMP in Full Info Eq?
    mp_funs[:mp_fi_rm_diff] = (x, y) -> (mp_funs[:rm_cond](x, y) != fi_funs[:rm_cond](x, y)) ? 1. : 0.

    for zvar in [z for z in contour_zvars if z != :MBR]
        zsym = contour_zvars_sym[zvar]
        mp_funs[Symbol(:rm_, zsym)], mp_funs[Symbol(:nrm_, zsym)] = get_rm_payoff_funs(mp_fd, mp_fd[:xvar],
                                                                                       mp_fd[:yvar], zvar)
        mp_funs[zsym] = (x, y) -> mp_funs[:rm_cond](x, y) ? mp_funs[Symbol(:rm_, zsym)](x, y) : mp_funs[Symbol(:nrm_, zsym)](x, y)
        mp_funs[Symbol(:mp_fi_, zsym, :_diff)] = (x, y) -> mp_funs[zsym](x, y) - fi_funs[zsym](x, y)
    end

    return mp_funs
end


function jeq_payoff_functions(fi_funs::Dict{Symbol, Any}, jfd::Dict;
                              eq_type::String="pooling",
                              xvar::Symbol=:iota, yvar::Symbol=:sigmah)


    if any([fi_funs[:xvar] != xvar, fi_funs[:yvar] != yvar])
        println("Full Info and Misrep x and y variables do not coincide. Exiting...")
        return
    end

    jeq_funs = Dict{Symbol, Any}(:xvar => xvar, :yvar => yvar,
                                 :safe => Dict{Symbol, Any}(),
                                 :risky => Dict{Symbol, Any}())

    r_obj_fun = :MBR
    if eq_type == "separating"
        r_obj_fun = :firm_value
    end


    # Set objective function symbol
    zsym = contour_zvars_sym[r_obj_fun]

    # Risky Firm's Objective Function Payoff in case of Risk-Management v.s. No Risk-Management
    jeq_funs[:risky][Symbol(:rm_, zsym)], jeq_funs[:risky][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar, Symbol(:r_, r_obj_fun))

    # Choose Risk-Management if it maximizes Payoff
    jeq_funs[:risky][:rm_cond] = (x, y) -> jeq_funs[:risky][Symbol(:rm_, zsym)](x, y) >= jeq_funs[:risky][Symbol(:nrm_, zsym)](x, y)

    # Does this RMP differ from the optimal RMP in Full Info Eq?
    jeq_funs[:risky][:jeq_fi_rm_diff] = (x, y) -> (jeq_funs[:risky][:rm_cond](x, y) != fi_funs[:rm_cond](x, y)) ? 1. : 0.

    # Risky Firm's Objective Function Payoff
    jeq_funs[:risky][zsym] = (x, y) -> maximum([jfd[Symbol(:r_, r_obj_fun)][xvar](x),
                                                jfd[Symbol(:r_, r_obj_fun)][yvar](y)])

    # Difference in RMP between Joint and Full Information Equilibrium
    jeq_funs[:risky][Symbol(:jeq_fi_, zsym, :_diff)] = (x, y) -> jeq_funs[:risky][zsym](x, y) - fi_funs[zsym](x, y)

    # Safe Firm's Payoff Depends on what Risky Firm chooses
    jeq_funs[:safe][Symbol(:rm_, zsym)], jeq_funs[:safe][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar,
                                                                                                   Symbol(:s_, r_obj_fun))
    jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[:safe][Symbol(:rm_, zsym)](x, y) : jeq_funs[:safe][Symbol(:nrm_, zsym)](x, y)

    # jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:safe][zsym](x, y) - Need FI payoff #fi_funs[zsym](x, y)

    for zvar in [z for z in contour_zvars if z != r_obj_fun]
        zsym2 = contour_zvars_sym[zvar]

        for ft in keys(contour_firm_types)
            ft_z = contour_firm_types[ft]
            jeq_funs[ft][Symbol(:rm_, zsym2)], jeq_funs[ft][Symbol(:nrm_, zsym2)] = get_rm_payoff_funs(jfd, jfd[:xvar],
                                                                                                     jfd[:yvar],
                                                                                                     Symbol(ft_z, zvar))
            jeq_funs[ft][zsym2] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[ft][Symbol(:rm_, zsym2)](x, y) : jeq_funs[ft][Symbol(:nrm_, zsym2)](x, y)
        end

        jeq_funs[:risky][Symbol(:jeq_fi_, zsym2, :_diff)] = (x, y) -> jeq_funs[:risky][zsym2](x, y) - fi_funs[zsym2](x, y)
    end


    return deepcopy(jeq_funs)
end


function get_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs,
                                     fi_fv::Float64, fi_fv_fun, k_otc::Float64)
    jeq_ind = (x, y) -> mp_funs[:mbr](x, y) >= fi_funs[:mbr](x, y)
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:safe][:fv](x, y) >= sep_funs[:safe][:fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false

    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :mbr => Dict{Symbol, Any}())
    for zvar in [:fv, :mbr, :lev]
        fun_dict[:mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[zvar](x, y) : .0
        fun_dict[:mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[:risky][zvar](x, y) : .0
        fun_dict[:mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[:risky][zvar](x, y) : .0
    end

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) +
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) +
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:risky][:mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:risky][:mbr](x, y))
    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv + #fi_funs[:fv](x, y) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:safe][:fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:safe][:fv](x, y))
    fun_dict[:bool_otc] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : 0
    fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : fun_dict[:eq_bool](x, y)

    catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)],
                    [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
    fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]

    return fun_dict
end

# ** NRMP Functions {{{2

# *** Check Conditions {{{3
function nrmp_check_cvm_cond(df::DataFrame, var::Symbol, value::Float64; tol::Float64=1e-5)
    # CVM lambda is either NaN or zero
    cvm_lambda_cond = .&(any([isnan(value), abs.(value) < tol]),
                         occursin("lambda", string(var)))

    # CVM sigmah is either NaN or equal to sigmal
    cvm_sigmah_cond = .&(any([isnan(value), abs(minimum(df[:, :sigmal]) - value) < tol]),
                         occursin("sigmah", string(var)))

    return any([cvm_lambda_cond, cvm_sigmah_cond])
end

# *** Identify Equilibrium and Firm Types {{{3
function get_eq_ft_types(fd, ft::Symbol)
    df_eq_type = Symbol(fd[:df][1, :eq_type])

    # Equilibrium Type
    eq_type = Symbol("")
    if df_eq_type in [:fi, :full_info, :full_information]
        eq_type = :fi
    elseif df_eq_type in [:mp, :misrep, :misrepresentation]
        eq_type = :mp
    elseif df_eq_type in [:pool, :pooling]
        eq_type = :pool
    elseif df_eq_type in [:sep, :separating]
        eq_type = :sep
    else
        println("Error! Could not identify equilibrium type. Returning...")
        return
    end

    # Check firm type:
    if ft in [:st, :safe]
        ft = :st
    elseif ft in [:rt, :risky]
        ft = :rt
    elseif !isempty(String(ft))
        println("wrong ft argument. Please enter :st (:safe), :rt (:risky) or empty Symbol.")
        println("Exiting...")
        return
    end

    return eq_type, ft
end

# *** Extract Z Variable {{{3
function nrmp_return_zval(df::DataFrame, xvar::Symbol, xval::Float64,
                          yvar::Symbol, yval::Float64,
                          zvar::Symbol, eq_type::Symbol;
                          misrep_type::Symbol=:rt, tol::Float64=1e-5)

    cvm_cond1 = nrmp_check_cvm_cond(df, xvar, xval; tol=tol)
    cvm_cond2 = nrmp_check_cvm_cond(df, yvar, yval; tol=tol)
    cvm_cond = cvm_cond1 || cvm_cond2

    eq_cond = (eq_type in [:misrep, :mp, :separating, :sep])

    if cvm_cond
        return unique(nrmp_cvm_df_slicer(df; tol=tol)[:, zvar])[1]
    else
        tmp(var, val) = isnan(val) ? isnan.(df[:, var]) : abs.(df[:, var] .- val) .< tol

        println(string(xvar, ": ", xval, "; ", yvar, ": ", yval))
        if .&(xvar == :mu_s, eq_cond)
            cond = .&(tmp(yvar, yval), df[:, :misrep_type] .== misrep_type)
        elseif .&(yvar == :mu_s, eq_cond)
            cond = .&(tmp(xvar, xval), df[:, :misrep_type] .== misrep_type)
        elseif .&(xvar in Symbol.(names(df)), yvar in Symbol.(names(df)))
            cond  = .&(tmp(xvar, xval), tmp(yvar, yval))
        elseif xvar in Symbol.(names(df))
            cond = tmp(xvar, xval)
        else
            cond = tmp(yvar, yval)
        end

        return df[cond, zvar][1]
    end
end


function nrmp_contour_names_zvar(zvar::Symbol, eq_type::Symbol;
                                 ft::Symbol=Symbol(""),
                                 diff::Bool=false,
                                 diff_percent::Bool=false)

    ft_prefix = !isempty(string(ft)) ? Symbol(string(ft)[1], :_) : Symbol("")

    zvar_adj = Symbol("")
    if zvar in [:fv, :firm_value]
        zvar_adj = :fv
    elseif zvar in [:mbr, :MBR]
        zvar_adj = :mbr
    elseif zvar in [:yield, :yield_spd]
        zvar_adj = zvar
    else
        println("Error! Could not identify zvar. Returning...")
        return
    end

    fname_eq_type = any([diff, diff_percent]) ? Symbol(eq_type, :_fi) : eq_type
    fname_diff = diff ? :_diff : Symbol("")
    fname_diff = diff_percent ? :_perc_diff : fname_diff

    fname_zvar = Symbol(ft_prefix, fname_eq_type, :_, zvar_adj, fname_diff)


    plot_zvar = fname_zvar
    if !occursin("diff", string(fname_zvar))
        plot_zvar =  Symbol(ft_prefix, zvar_adj)
    end

    return fname_zvar, plot_zvar
end

# *** Additional Z Vars and Surface Interpolation {{{3
function nrmp_surface_interpolator(sd;
                                   zvars::Array{Symbol,1}=Array{Symbol, 1}(),
                                   non_zvars::Array{Symbol,1}=[:df, :xvar, :xvals, :X, :yvar, :yvals, :Y, :interp],
                                   kx::Int64=3, ky::Int64=3) #, s::Float64=.25)
    if isempty(zvars)
        zvars = [x for x in keys(sd) if !(x in non_zvars)]
    end

    sd[:interp] = Dict{Symbol, Any}([zvar => nothing for zvar in zvars])

    for zvar in zvars
        sd[:interp][zvar] = Dierckx.Spline2D(sd[:xvals], sd[:yvals],
                                             sd[zvar][:zvals], kx=kx, ky=ky) #, s=s)
    end

    return sd
end


function nrmp_contour_plot_additional_zvars(fid, resd, eq_type::Symbol)

    if eq_type in [:misrep, :misrepresentation]
        zvar = :mbr
        resd[:r_mp_fi_mbr_diff] = Dict(:zvals => resd[Symbol(:r_, zvar)][:zvals]
                                      .- fid[zvar][:zvals][2:end, 2:end])
        resd[:r_mp_fi_mbr_perc_diff] = Dict(:zvals => (resd[:r_mp_fi_mbr_diff][:zvals]
                                                      ./ fid[zvar][:zvals][2:end, 2:end]) .* 100 )

        zvar = :fv
        resd[:r_mp_fi_fv_diff] = Dict(:zvals => resd[Symbol(:r_, zvar)][:zvals]
                                       .- fid[zvar][:zvals][2:end, 2:end])
        resd[:r_mp_fi_fv_perc_diff] = Dict(:zvals => (resd[:r_mp_fi_fv_diff][:zvals]
                                                       ./ fid[zvar][:zvals][2:end, 2:end]) .* 100 )
    else
        eqt = Symbol("")
        if eq_type in [:pool, :pooling]
            eqt = :pool
        elseif  eq_type in [:sep, :separating]
            eqt = :sep
        else
            println("Error! Equilibrium type not identified. Exiting...")
            return
        end


        zvar = :mbr
        zvar1 = Symbol(:r_, eqt, :_fi_mbr_diff)
        resd[zvar1] = Dict(:zvals => resd[Symbol(:r_, zvar)][:zvals]
                           .- fid[zvar][:zvals][2:end, 2:end])

        zvar2 = Symbol(:r_, eqt, :_fi_mbr_perc_diff)
        resd[zvar2] = Dict(:zvals => (resd[zvar1][:zvals]
                                      ./ fid[zvar][:zvals][2:end, 2:end]) .* 100)

        zvar3 = Symbol(:s_, eqt, :_fi_mbr_diff)
        resd[zvar3] = Dict(:zvals => resd[Symbol(:s_, zvar)][:zvals]
                           .- fid[zvar][:zvals][1, 1])

        zvar4 = Symbol(:s_, eqt, :_fi_mbr_perc_diff)
        resd[zvar4] = Dict(:zvals => (resd[zvar3][:zvals]
                                      ./ fid[zvar][:zvals][1, 1]) .* 100)

        zvar = :fv
        zvar5 = Symbol(:s_, eqt, :_fi_fv_diff)
        resd[zvar5] = Dict(:zvals => resd[Symbol(:s_, zvar)][:zvals]
                           .- fid[zvar][:zvals][1, 1])

        zvar6 = Symbol(:s_, eqt, :_fi_fv_perc_diff)
        resd[zvar6] = Dict(:zvals => (resd[zvar5][:zvals]
                                      ./ fid[zvar][:zvals][1, 1]) .* 100)
    end

    return nrmp_surface_interpolator(resd)
end


function get_nrmp_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs,
                                          k_ep::Float64, k_otc::Float64, fi_fv_fun)
    jeq_ind = (x, y) -> mp_funs[:r_mbr](x, y) >= fi_funs[:mbr](x, y)
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:s_fv](x, y) >= sep_funs[:s_fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false

    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :r_mbr => Dict{Symbol, Any}())

    fun_dict[:r_mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[:mbr](x, y) : .0
    fun_dict[:r_mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[:r_mbr](x, y) : .0
    fun_dict[:r_mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[:r_mbr](x, y) : .0

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) +
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) +
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:r_mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:r_mbr](x, y))

    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv_fun(k_ep) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:s_fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:s_fv](x, y))

    fun_dict[:bool_otc] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : 0
    fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : fun_dict[:eq_bool](x, y)

    catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)],
                    [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
    fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]

    return fun_dict
end

# *** Prepare and Slice DataFrames {{{3
function nrmp_cvm_df_slicer(df::DataFrame; tol::Float64=1e-5)
    # CVM sigmah is either NaN or equal to sigmal
    sigmah_var = :nrm_sigmah
    cvm_sigmah_cond = isnan.(df[:, sigmah_var]) .| (abs.(df[:, sigmah_var] .-
                                                 minimum(df[:, :sigmal])) .< tol)

    # CVM lambda is either NaN or zero
    lambda_var = :nrm_lambda
    cvm_lambda_cond = isnan.(df[:, sigmah_var]) .| (abs.(df[:, lambda_var]) .< tol)

    return df[cvm_sigmah_cond .| cvm_lambda_cond, :]
end

function nrmp_contour_plot_prepare_df(resdf::DataFrame;
                                      zvars::Array{Symbol,1}=contour_zvars,
                                      misrep_type::Symbol=:rt,
                                      min_lambda::Float64=NaN,
                                      max_lambda::Float64=NaN,
                                      min_sigmah::Float64=NaN,
                                      xvar::Symbol=:nrm_sigmah,
                                      xvals::Array{Float64, 1}=Array{Float64, 1}(),
                                      yvar::Symbol=:nrm_lambda,
                                      yvals::Array{Float64, 1}=Array{Float64, 1}(),
                                      tol::Float64=1e-5)

    df = deepcopy(resdf)
    eq_type = Symbol(unique(df[:, :eq_type])[1])

    type_prefix = Symbol("")
    ft_vec = [Symbol("")]

  if !(eq_type in [:full_info, :fi])
      xvar = xvar != :mu_s ? Symbol(:r_, xvar) : xvar
      yvar = yvar != :mu_s ? Symbol(:r_, yvar) : yvar

      ft_vec = [:s_, :r_]

      df = df[df[:, :misrep_type] .== misrep_type, :]
  end

    # Form Z variable names
    zvar_vec = vcat([Symbol(ft, x) for ft in ft_vec, x in zvars]...)
    zf(x) = x == :firm_value ? :fv : Symbol(lowercase(string(x)))
    zvar_name_vec = vcat([Symbol(ft, zf(x)) for ft in ft_vec, x in zvars]...)
    zvardf = DataFrame( :zvar => zvar_vec, :zvar_name=> zvar_name_vec)

    for ft in ft_vec
        if !isnan(min_lambda)
            col = Symbol(ft, :nrm_lambda)
            df[isnan.(df[:, col]), col] .= min_lambda
        end
        if !isnan(max_lambda)
            col = Symbol(ft, :nrm_lambda)
            df = df[df[:, col] .<= max_lambda, :]
        end
        if !isnan(min_sigmah)
            col = Symbol(ft, :nrm_sigmah)
            df[isnan.(df[:, col]), col] .= min_sigmah
        end
    end

    eq_cond = eq_type in [:misrep, :mp, :separating, :sep]
    if xvar in Symbol.(names(df))
        xvals = .&(xvar == :mu_s, eq_cond) ? xvals : unique(df[:, xvar])
    end
    if yvar in Symbol.(names(df))
        yvals = .&(yvar == :mu_s, eq_cond) ? yvals : unique(df[:, yvar])
    end

    # Adjust x and y values
    xvals = sort([x for x in xvals if !isnan(x)])
    yvals = sort([y for y in yvals if !isnan(y)])

    # Store Results
    resd = Dict{Symbol, Any}(:df => df,
                             :xvar => xvar,
                             :xvals => sort(xvals),
                             :X => repeat(xvals, 1, size(yvals, 1)),
                             :yvar => yvar,
                             :yvals => sort(yvals),
                             :Y => Array(repeat(yvals, 1, size(xvals, 1))'))

    # Function to extract zvar value:
    # zfun(xval, yval) = df[.&(abs.(df[:, xvar] .- xval) .< tol,
    #                          abs.(df[:, yvar] .- yval) .< tol), zvar][1]

    for row in 1:size(zvardf, 1)
        zvar = zvardf[row, :zvar]
        zvar_name = zvardf[row, :zvar_name]

        resd[zvar_name] = Dict{Symbol, Any}(:zvals => Array{Float64,2}[],
                                            :fun => nothing)

        resd[zvar_name][:zvals] = [nrmp_return_zval(df, xvar, xval, yvar, yval, zvar, eq_type;
                                                    misrep_type=misrep_type, tol=tol)
                                   for xval in resd[:xvals], yval in resd[:yvals]]

        resd[zvar_name][:fun] = Dierckx.Spline2D(resd[:xvals], resd[:yvals], resd[zvar_name][:zvals])
    end

    return resd
end

function nrmp_remove_extra_vals(df::DataFrame)
    eq_type = df[1, :eq_type]
    fi_eq_cond = eq_type in [:fi, :full_info]
    prefix = fi_eq_cond ? Symbol("") : :r_

    for xy in [:sigmah, :lambda]
        xy_col = Symbol(prefix, :nrm_, xy)
        val_vec = svm_param_values_dict[xy]
        if .&(fi_eq_cond, xy == :sigmah)
            val_vec = vcat(svm_param_values_dict[:sigmal], val_vec)
        end

        missing_vals = [x for x in unique(df[:, xy_col]) if .&(!(x in val_vec), !isnan.(x))]
        println(missing_vals)
        for val in missing_vals
            cond = isnan.(df[:, xy_col]) .|  (abs.(df[:, xy_col] .- val) .> 1e-4)
            df = df[cond, :]
        end
    end

    return df
end

# *** File Path {{{3
function get_nrmp_contour_plot_path_name(df::DataFrame, fname_zvar::Symbol;
                                         xvar::Symbol=:sigmah, yvar::Symbol=:lambda,
                                         eq_type::Symbol=Symbol(""),
                                         fname_ext::String=contour_fname_ext)


    fixed_vars_vec = [:mu_s, :lambda]
    fixed_var = [x for x in fixed_vars_vec if !(x in [xvar, yvar])][1]

    # Equilibrium Type
    if isempty(string(eq_type))
        str_eq_type = string(df[1, :eq_type])
        if str_eq_type == "pool"
            str_eq_type = "pooling"
        elseif str_eq_type == "sep"
            str_eq_type = "separating"
        end
        eq_type = Symbol(eq_type_title[str_eq_type][1])
    end

    # Inverse Coupon Rate
    pcr = df[1, :p]/df[1, :c]

    # Misrepresenting Type
    mft=Symbol("")
    mu_s = NaN
    lambda_val = NaN
    if (eq_type in [:fi, :full_info])
        lambda_val = df[1, :nrm_lambda]
    else
        mft = df[1, :misrep_type]
        lambda_val = df[1, :r_nrm_lambda]

        # Measure of Safe Firms
        if !(:mu_s in [xvar, yvar])
            mu_s = df[1, :mu_s] # (eq_type == :pool) ? df[1, :mu_s] : mu_s

            # Misrepresentation and Separating cases
            if (eq_type in [:misrep, :misrepresentation, :sep, :separating]) # [:pool, :pooling])
                st_misrep_cond = .&(mft == :st, abs(mu_s) < 1e-5)
                rt_misrep_cond = .&(mft == :rt, abs(mu_s - 1.) < 1e-5)
                if .&(!st_misrep_cond, !rt_misrep_cond)
                    println("Error! Misrepresentation type and mu_s don't match!")
                    println("Exiting...")
                    return
                end
            end
        end
    end
    mft_name = !isempty(string(mft)) ? string(:_misrep_, mft) : ""
    mu_s_name = !isnan(mu_s) ? string(:_mu_s_, mu_s) : ""
    lambda_name = !isnan(lambda_val) ? string(:_lambda_, lambda_val) : ""

    fixed_var_name = fixed_var == :lambda ? lambda_name : mu_s_name

    # File Name
    fname = string("nrmp_", eq_type, mft_name, fixed_var_name, "_pcr_", pcr,
                   "_x_", xvar, "_y_", yvar,
                   "__", fname_zvar)

    return string(contour_plots_path, "/", fname, ".", fname_ext)
end

# *** Plot Title {{{3
function get_nrmp_contour_plot_title(df::DataFrame,
                                     eq_type::Symbol,
                                     zvar::Symbol;
                                     k_otc::Float64=NaN,
                                     xvar::Symbol=:sigmah,
                                     yvar::Symbol=:lambda,
                                     ft::Symbol=Symbol(""),
                                     diff::Bool=false,
                                     diff_percent::Bool=false,
                                     params_list::Array{Symbol,1}=contour_plots_title_params_order,
                                     s_fi_fv::Float64=NaN,
                                     s_otc_fv::Float64=NaN)

    title_eq_type = [eq_type_title[k][2]  for k in keys(eq_type_title) if (eq_type_title[k][1] == eq_type)][1]

    params_list = [x for x in params_list if !(x in[xvar, yvar])]

    mu_s = NaN
    if .&(eq_type in [:pool, :pooling, :ep_market, :dual_market], !(:mu_s in [xvar, yvar])) #, "separating"]
        mu_s =  df[1, :mu_s]
        params_list = vcat(:mu_s, params_list)
    end
    df[:, :pcr] .= df[1, :p]/df[1, :c]

    title_params = join([string("\$", contour_tlabels[x][1], "= \$ ",
                                str_format_fun(contour_tlabels[x][2],
                                               get_title_value(df, x)))
                         for x in params_list], ", ")

    firm_type_title = ""

    if ft == :st
        firm_type_title = "Safe Type's "
    else
        firm_type_title = "Risky Type's "
    end

    if diff | diff_percent
        differential = diff_percent ? " Differential (%) " : " Differential "

        plot_title = latexstring(firm_type_title, title_eq_type,
                                 " v.s. Full Information Eq. ",
                                 contour_tlabels[zvar][1], differential,
                                 "\n for ", title_params)
    else
        if eq_type in [:mp, :misrep, :misrepresentation]
            plot_title = latexstring(firm_type_title,
                                     contour_tlabels[zvar][1],
                                     " in case of Misrepresentation",
                                     "\n for ", title_params)
        elseif eq_type in [:pool, :pooling, :sep, :separating]
            plot_title = latexstring(firm_type_title, "Optimal ",
                                     contour_tlabels[zvar][1], " in a ",
                                     title_eq_type, " Equilibrium ",
                                     "\n for ", title_params)
        else
            k_otc_title = ""
            if !isnan(k_otc)
                k_otc_title = string(", \$", contour_tlabels[:kappa_otc][1], "= \$ ",
                                     str_format_fun(contour_tlabels[:kappa_otc][2], k_otc))
            end

            plot_title = latexstring(firm_type_title, "Optimal ",
                                     contour_tlabels[zvar][1], " in the ",
                                     title_eq_type,
                                     "\n for ", title_params, k_otc_title)
        end
    end


    if .&(!isnan(s_fi_fv), !isnan(s_otc_fv))
        plot_title = latexstring(plot_title,
                                 "\n Safe Type's Full Information Firm Value = ",
                                 str_format_fun("%.2f", s_fi_fv),
                                 ", Safe Type's OTC Firm Value = ",
                                 str_format_fun("%.2f", s_otc_fv), " \n")

    elseif !isnan(s_fi_fv)
        plot_title = latexstring(plot_title,
                                 "\n (Safe Type's Full Information Firm Value = ",
                                 str_format_fun("%.2f", s_fi_fv), ") \n")
    elseif !isnan(s_otc_fv)
        plot_title = latexstring(plot_title,
                                 "\n (Safe Type's OTC Firm Value = ",
                                 str_format_fun("%.2f", s_otc_fv), ") \n")
    end

    return plot_title
end

# *** Plot Functions {{{3
function plot_nrmp_iso_contour_curves(fd::Dict{Symbol, Any},
                                      zvar;
                                      rfd::Dict=Dict(),
                                      iso_xlabel::Symbol=:sigmah,
                                      iso_ylabel::Symbol=:lambda,
                                      rmp_Z::Array{Float64, 2}=Array{Float64,2}(undef, 0, 2),
                                      diff_fun=nothing,
                                      fig_title::LaTeXString=LaTeXString(""),
                                      file_path_name::String="",
                                      seaborn_style=iso_plt_inputs[:seaborn_style],
                                      iso_levels=iso_plt_inputs[:iso_levels],
                                      heat_levels=iso_plt_inputs[:heat_levels],
                                      iso_cmap=iso_cmaps["full_info"],
                                      heat_cmap::String="",
                                      fig_aspect=iso_plt_inputs[:fig_aspect],
                                      iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                      use_subgrid=iso_plt_inputs[:use_subgrid],
                                      subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                      iso_cols=iso_plt_inputs[:iso_cols],
                                      heat_cols=iso_plt_inputs[:heat_cols],
                                      title_font_size=iso_plt_inputs[:title_font_size],
                                      fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                      tight_pad=iso_plt_inputs[:tight_pad],
                                      h_pad=iso_plt_inputs[:h_pad],
                                      w_pad=iso_plt_inputs[:w_pad],
                                      cat_Z=[],
                                      cat_cmap="GnBu",
                                      cat_alpha::Float64=.25,
                                      rmp_cmap=rmp_cmaps["misrep"],
                                      rmp_alpha::Float64=.15)


    fig, ax1, ax2 = plot_iso_curves(fd[:X], fd[:Y], fd[zvar][:zvals],
                                    iso_xlabel=iso_xlabel,
                                    iso_ylabel=iso_ylabel,
                                    rmp_Z=rmp_Z,
                                    seaborn_style=seaborn_style,
                                    iso_levels=iso_levels,
                                    heat_levels=heat_levels,
                                    iso_cmap=iso_cmap,
                                    heat_cmap=heat_cmap,
                                    fig_aspect=fig_aspect,
                                    iso_fontsize=iso_fontsize,
                                    use_subgrid=use_subgrid,
                                    subgrid_rows=subgrid_rows,
                                    iso_cols=iso_cols,
                                    heat_cols=heat_cols,
                                    cat_Z=cat_Z,
                                    cat_cmap=cat_cmap,
                                    cat_alpha=cat_alpha,
                                    rmp_cmap=rmp_cmap,
                                    rmp_alpha=rmp_alpha)

    if !isempty(rfd)
        for rf in keys(rfd)
            ax1.plot(rfd[rf][:x], rfd[rf][:y], color=rfd[rf][:color],
                     marker=rfd[rf][:marker], markersize=rfd[rf][:ms])
            ax1.annotate(rfd[rf][:name],
                         fontsize=rfd[rf][:fontsize],
                         color=rfd[rf][:color],
                         xy=(rfd[rf][:x], rfd[rf][:y]),
                         xytext=(rfd[rf][:xf] * rfd[rf][:x],
                                 rfd[rf][:yf] * rfd[rf][:y]))
        end
    end

    if !isempty(fig_title)
        fig.suptitle(fig_title, fontsize=title_font_size)
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig, ax1, ax2
end


function plot_nrmp_contour(fd::Dict{Symbol, Any}, zvar::Symbol;
                           ft::Symbol=Symbol(""),
                           rfd::Dict=Dict(),
                           iso_xlabel::Symbol=:sigmah,
                           iso_ylabel::Symbol=:lambda,
                           rmp_Z::Array{Float64, 2}=Array{Float64,2}(undef, 0, 2),
                           diff::Bool=false,
                           diff_percent::Bool=false,
                           fig_title::LaTeXString=LaTeXString(""),
                           file_path_name::String="",
                           fname_ext::String=contour_fname_ext,
                           seaborn_style=iso_plt_inputs[:seaborn_style],
                           iso_levels=iso_plt_inputs[:iso_levels],
                           heat_levels=iso_plt_inputs[:heat_levels],
                           iso_cmap=iso_cmaps["full_info"],
                           heat_cmap::String="",
                           fig_aspect=iso_plt_inputs[:fig_aspect],
                           iso_fontsize=iso_plt_inputs[:iso_fontsize],
                           use_subgrid=iso_plt_inputs[:use_subgrid],
                           subgrid_rows=iso_plt_inputs[:subgrid_rows],
                           iso_cols=iso_plt_inputs[:iso_cols],
                           heat_cols=iso_plt_inputs[:heat_cols],
                           title_font_size=iso_plt_inputs[:title_font_size],
                           fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                           tight_pad=iso_plt_inputs[:tight_pad],
                           h_pad=iso_plt_inputs[:h_pad],
                           w_pad=iso_plt_inputs[:w_pad],
                           cat_Z=[],
                           cat_cmap="GnBu",
                           cat_alpha::Float64=.25,
                           rmp_cmap=rmp_cmaps["misrep"],
                           rmp_alpha::Float64=.15)

    if zvar == :firm_value
        zvar = :fv
    elseif zvar == :MBR
        zvar = :mbr
    end

    # Equilibrium and firm type:
    eq_type, ft = get_eq_ft_types(fd, ft)

    # Figure Title
    if isempty(fig_title)
        fig_title = get_nrmp_contour_plot_title(fd[:df], eq_type, zvar; ft=ft,
                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                diff=diff, diff_percent=diff_percent)
    end

    # Z Variable Names
    fname_zvar, plot_zvar = nrmp_contour_names_zvar(zvar, eq_type;
                                                    ft=ft, diff=diff,
                                                    diff_percent=diff_percent)


    if isempty(file_path_name)
        file_zvar = isempty(rfd) ? fname_zvar : Symbol(fname_zvar, :_slide)

        file_path_name = get_nrmp_contour_plot_path_name(fd[:df], file_zvar;
                                                         xvar=iso_xlabel, yvar=iso_ylabel,
                                                         eq_type=eq_type,
                                                         fname_ext=fname_ext)
    end


    println(file_path_name)
    fig, ax1, ax2 = plot_nrmp_iso_contour_curves(fd, plot_zvar;
                                                 rfd=rfd,
                                                 iso_xlabel=iso_xlabel,
                                                 iso_ylabel=iso_ylabel,
                                                 fig_title=fig_title,
                                                 file_path_name=file_path_name,
                                                 iso_cmap=iso_cmap,
                                                 seaborn_style=seaborn_style,
                                                 iso_levels=iso_levels,
                                                 heat_levels=heat_levels,
                                                 heat_cmap=heat_cmap,
                                                 fig_aspect=fig_aspect,
                                                 iso_fontsize=iso_fontsize,
                                                 use_subgrid=use_subgrid,
                                                 subgrid_rows=subgrid_rows,
                                                 iso_cols=iso_cols,
                                                 heat_cols=heat_cols,
                                                 title_font_size=title_font_size,
                                                 fig_dpi=fig_dpi,
                                                 tight_pad=tight_pad,
                                                 h_pad=h_pad,
                                                 w_pad=w_pad,
                                                 cat_Z=cat_Z,
                                                 cat_cmap=cat_cmap,
                                                 cat_alpha=cat_alpha,
                                                 rmp_cmap=rmp_cmap,
                                                 rmp_alpha=rmp_alpha)

    return fig, ax1, ax2
end


function plot_nrmp_equilibria_contour(pltd, eq_type_Z, zvar,
                                      eq_type::Symbol, df::DataFrame;
                                      ft::Symbol=Symbol(""),
                                      iso_xlabel::Symbol=:sigmah,
                                      iso_ylabel::Symbol=:lambda,
                                      rmp_Z::Array{Float64, 2}=Array{Float64,2}(undef, 0, 2),
                                      diff::Bool=false,
                                      diff_percent::Bool=false,
                                      fig_title::LaTeXString=LaTeXString(""),
                                      file_path_name::String="",
                                      fname_ext::String=contour_fname_ext,
                                      seaborn_style=iso_plt_inputs[:seaborn_style],
                                      iso_levels=iso_plt_inputs[:iso_levels],
                                      heat_levels=iso_plt_inputs[:heat_levels],
                                      iso_cmap=iso_cmaps["full_info"],
                                      heat_cmap::String="",
                                      fig_aspect=iso_plt_inputs[:fig_aspect],
                                      iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                      use_subgrid=iso_plt_inputs[:use_subgrid],
                                      subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                      iso_cols=iso_plt_inputs[:iso_cols],
                                      heat_cols=iso_plt_inputs[:heat_cols],
                                      title_font_size=iso_plt_inputs[:title_font_size],
                                      fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                      tight_pad=iso_plt_inputs[:tight_pad],
                                      h_pad=iso_plt_inputs[:h_pad],
                                      w_pad=iso_plt_inputs[:w_pad],
                                      cat_Z=[],
                                      cat_cmap="GnBu",
                                      cat_alpha::Float64=.25,
                                      rmp_cmap=rmp_cmaps["misrep"],
                                      rmp_alpha::Float64=.15)

    if zvar == :fv
        zvar = :firm_value
    elseif zvar == :mbr
        zvar = :MBR
    end

    # Equilibrium and firm type:
    # eq_type, ft = get_eq_ft_types(fd, ft)

    # Figure Title
    if isempty(fig_title)
        fig_title = get_nrmp_contour_plot_title(df, eq_type, zvar; ft=ft,
                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                diff=diff, diff_percent=diff_percent)
    end

    # Z Variable Names
    fname_zvar, plot_zvar = nrmp_contour_names_zvar(zvar, eq_type;
                                                    ft=ft, diff=diff,
                                                    diff_percent=diff_percent)


    if isempty(file_path_name)
        file_path_name = get_nrmp_contour_plot_path_name(df, fname_zvar;
                                                         xvar=iso_xlabel, yvar=iso_ylabel,
                                                         eq_type=eq_type,
                                                         fname_ext=fname_ext)
    end

    fig, ax1, ax2 = plot_iso_curves(pltd[:X], pltd[:Y], pltd[zvar_name],
                                    cat_Z=pltd[eq_type_Z],
                                    iso_xlabel=iso_xlabel,
                                    iso_ylabel=iso_ylabel,
                                    rmp_Z=rmp_Z,
                                    seaborn_style=seaborn_style,
                                    iso_levels=iso_levels,
                                    heat_levels=heat_levels,
                                    iso_cmap=iso_cmap,
                                    heat_cmap=heat_cmap,
                                    fig_aspect=fig_aspect,
                                    iso_fontsize=iso_fontsize,
                                    use_subgrid=use_subgrid,
                                    subgrid_rows=subgrid_rows,
                                    iso_cols=iso_cols,
                                    heat_cols=heat_cols,
                                    cat_cmap=cat_cmap,
                                    cat_alpha=cat_alpha,
                                    rmp_cmap=rmp_cmap,
                                    rmp_alpha=rmp_alpha)

    if !isempty(fig_title)
        fig.suptitle(fig_title, fontsize=title_font_size)
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end




end

# ** OTC Payoff Functions {{{2
function interpolate_svmdf(df::DataFrame,
                           xvar::Symbol, yvar::Symbol;
                           obj_fun::Symbol=:firm_value,
                           zvars::Array{Symbol, 1}=[:MBR, :firm_value],
                           m::Float64=1.,
                           fixed_var::Symbol=Symbol(""), fixed_val::Float64=NaN,
                           tol::Float64=1e-5,
                           s_interp::Float64=1.)


    # Filter DataFrame ========================
    # Objective Function
    cond = Symbol.(df[:, :obj_fun]) .== obj_fun

    # Maturity
    cond = .&(cond, abs.(df[:, :m] .- m) .< tol)

    # Additional variable
    if .&(fixed_var in Symbol.(names(df)), !isnan(fixed_val))
        cond = .&(cond, abs.(df[:, fixed_var] .- fixed_val) .< tol)
   end

    # Filter
    df = df[cond, :]
    # ==========================================

    # Check values for common parameters
    for col in [:gross_delta, :alpha, :V0, :r, :pi, :xi]
        if size(unique(df[:, col]), 1) > 1
            println(string("Error! More than one value for parameter ", col))
            return
        end
    end


    zd = Dict{Symbol, Any}(:xvar => xvar,
                           :xvals => unique(df[:, xvar]),
                           :yvar => yvar,
                           :yvals => unique(df[:, yvar]),
                           :kvals => unique(df[:, :kappa]))

    for k in zd[:kvals]
        cond = abs.(df[:, :kappa] .- k) .< tol
        tmp = df[cond, vcat(xvar, yvar, zvars...)]

        try
            try
                zfuns =fetch(@spawn [Dierckx.Spline2D(tmp[:, xvar], tmp[:, yvar],
                                                      tmp[:, zvar]; s=s_interp)
                                     for zvar in zvars])
                zd[Symbol(k)] = Dict{Symbol, Any}(zip(zvars, zfuns))
            catch
                zfuns =fetch(@spawn [Dierckx.Spline2D(tmp[:, xvar], tmp[:, yvar],
                                                      tmp[:, zvar]; s=1.5*s_interp)
                                     for zvar in zvars])
                zd[Symbol(k)] = Dict{Symbol, Any}(zip(zvars, zfuns))
            end
        catch
           zfuns =fetch(@spawn [Dierckx.Spline2D(tmp[:, xvar], tmp[:, yvar],
                                                      tmp[:, zvar]; s=size(tmp, 1))
                                     for zvar in zvars])
           zd[Symbol(k)] = Dict{Symbol, Any}(zip(zvars, zfuns))
        end
   end


#     Z = Array{Float64, 2}(unstack(df, xvar, yvar, :firm_value)[:, 2:end])
#     itp = Interpolations.interpolate(Z, BSpline(Cubic(Line(OnGrid()))))
#     sitp = Interpolations.scale(itp, xvals, yvals)

    return zd
 end

function get_svm_interp_vals(zd, zvar::Symbol,
                             xval::Float64, yval::Float64,
                             kappa::Float64)

    sym(k) = Symbol(zd[:kvals][k])
    zvals = fetch(@spawn [zd[sym(k)][zvar](xval, yval) for k in 1:size(zd[:kvals],1)])
    zfun = Dierckx.Spline1D(zd[:kvals], zvals; k=3, bc="nearest", s=0.0)

    return zfun(kappa)
end

function get_cvm_interp_vals(df, zvar::Symbol, kappa::Float64)

    # Safe Type's Firm Value as a function of Transaction Costs Kappa
    cond = abs.(df[:, :iota]) .< 1e-5
    zfun = Dierckx.Spline1D(df[cond, :kappa], df[cond, zvar];
                                 k=3, bc="extrapolate")

    return zfun(kappa)
end


function get_fi_interp_vals(cvmdf, svmdf,
                            xvar::Symbol, yvar::Symbol,
                            zvar::Symbol, kappa::Float64)

    zd = interpolate_svmdf(svmdf, xvar, yvar)
    svmZ = [get_svm_interp_vals(zd, zvar, x, y, kappa) for x in zd[:xvals], y in zd[:yvals]]
    cvmZ = get_cvm_interp_vals(cvmdf, zvar, kappa)

    # Add row for sigmah = sigmal
    return vcat([cvmZ for x in 1:size(zd[:yvals], 1)]',  svmZ)
end


function get_otc_payoff_funs(cvmdf::DataFrame, svmdf::DataFrame, kotc::Float64;
                             xvar::Symbol=:sigmah, yvar::Symbol=:lambda,
                             avals::Array{Float64,1}=Array{Float64}(undef, 0),
                             zvars::Array{Symbol, 1}=[:MBR, :firm_value],
                             m::Float64=1.)


    X = Array{Float64}(undef, 0, 2)
    Y = Array{Float64}(undef, 0, 2)

    sigvals = vcat(unique(cvmdf[:, :sigmal]), unique(svmdf[:, :sigmah]))

    # The other variable is either lambda or mu_s, but
    # mu_s is not in svmdf
    avals = isempty(avals) ? unique(svmdf[:, :lambda]) : avals

    # Sigmah grid
    S = Array{Float64,2}(repeat(sigvals, 1, size(avals, 1)))

    # Additional Variable grid
    L = Array{Float64,2}(repeat(avals', size(sigvals, 1), 1))

    X = xvar == :sigmah ? S : L
    Y = xvar == :sigmah ? L : S


    otcd = Dict{Symbol, Any}(:kotc => kotc,
                             :xvar => xvar,
                             :yvar => yvar,
                             :X => X,
                             :Y => Y)

    for zvar in zvars

        # Adjust name
        szvar = Symbol(lowercase(string(zvar)))
        szvar = szvar == :firm_value ? :fv : szvar

        otcd[szvar] = Dict(:zvals => get_fi_interp_vals(cvmdf, svmdf,
                                                       xvar, yvar, zvar, kotc))
    end

    return otcd
end

# * END MODULE {{{1
end
