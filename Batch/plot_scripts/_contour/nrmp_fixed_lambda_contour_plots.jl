main_path = "/home/artur/BondPricing"
run_modularity = false
run_fi = true
run_misrep = true
run_pool = true
run_sep = true
run_jeq = true

# ####################################################################
# ATTENTION!!! #######################################################
# ####################################################################
#  Z Variables name on the dataframes: [:firm_value, :MBR]
#  Z variable names on the dictionaries: [:fv, :mbr]
# ####################################################################

# * Plot Inputs
xvar = :sigmah
yvar = :mu_s
plot_xvar = :nrm_sigmah
plot_yvar = :mu_s
title_font_size = ModelPlots.iso_plt_inputs[:title_font_size]

# * Load Modules and Prepare Results
if run_modularity
    plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
    plot_script_name = "nrmp_fixed_lambda_contour_plots_init.jl"
    include(string(plot_script_path, "/", plot_script_name))
else
    module_path = string(main_path, "/", "Julia/modules/")
    modls = ["Batch", "ModelObj",
             "AnalyticFunctions", "BondPrInterp",
             "EqFinDiff", "JointEqStructs", "FullInfoEq",
             "ModelPlots", "JointEq"]
    for modl in modls
        include(string(joinpath(module_path, modl), ".jl"))
    end
end


iso_levels=12
# * Full Information Plots
if run_fi
    iso_cmap = ModelPlots.iso_cmaps["full_info"]

# ** FI - [Risky] Type's Firm Value
    zvar = :fv
    ft=Symbol("")
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(fid, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** FI - [Risky] Type's MBR
    zvar = :MBR
    ft=Symbol("")
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(fid, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)
end

# * Misrepresentation
if run_misrep
    iso_cmap = ModelPlots.iso_cmaps["misrep"]

# ** MP - Risky Type's MBR
    zvar = :MBR
    ft=:rt
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(mpd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** MP v.s. FI - Risky Type's MBR Differential
    zvar = :MBR
    ft = :rt
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(mpd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** MP v.s. FI - Risky Type's MBR Differential (%)
    zvar = :MBR
    ft = :rt
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(mpd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)

end

# * Pooling Payoffs
if run_pool
    iso_cmap = ModelPlots.iso_cmaps["pooling"]

# ** POOL - Risky Type's MBR
    zvar = :mbr
    ft=:rt
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** POOL - Pool v.s. Full Info - Risky Type's MBR Differential
    zvar = :MBR
    ft=:rt
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** POOL - Pool v.s. Full Info - Risky Type's MBR Differential (%)
    zvar = :MBR
    ft=:rt
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)

# ** POOL - Pool v.s. Full Info - Safe Type's MBR Differential
    zvar = :MBR
    ft=:st
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** POOL - Pool v.s. Full Info - Safe Type's MBR Differential (%)
    zvar = :MBR
    ft=:st
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)


# ** POOL - Safe Type's Firm Value
    zvar = :fv
    ft=:st
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** POOL - Pool v.s. Full Info - Safe Type's Firm Value Differential
    zvar = :fv
    ft=:st
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** POOL - Pool v.s. Full Info - Safe Type's Firm Value Differential (%)
    zvar = :fv
    ft=:st
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(poold, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)
end

# * Separating Payoffs
if run_sep
    iso_cmap = ModelPlots.iso_cmaps["separating"]

# ** SEP - Risky Type's MBR
    zvar = :MBR
    ft=:rt
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** SEP - Pool v.s. Full Info - Risky Type's MBR Differential
    zvar = :MBR
    ft=:rt
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** SEP - Pool v.s. Full Info - Risky Type's MBR Differential (%)
    zvar = :MBR
    ft=:rt
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)

# ** SEP - Pool v.s. Full Info - Safe Type's MBR Differential
    zvar = :MBR
    ft=:st
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** SEP - Pool v.s. Full Info - Safe Type's MBR Differential (%)
    zvar = :MBR
    ft=:st
    zvar_diff=false
    diff_percent=true
    iso_levels=iso_levels

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)


# ** SEP - Safe Type's Firm Value
    zvar = :fv
    ft=:st
    zvar_diff=false
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)


# ** SEP - Pool v.s. Full Info - Safe Type's Firm Value Differential
    zvar = :fv
    ft=:st
    zvar_diff=true
    diff_percent=false

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap)

# ** SEP - Pool v.s. Full Info - Safe Type's Firm Value Differential (%)
    zvar = :fv
    ft=:st
    zvar_diff=false
    diff_percent=true
    iso_levels=12

    fig, ax1, ax2 = ModelPlots.plot_nrmp_contour(sepd, zvar;
                                                 ft=ft,
                                                 diff=zvar_diff,
                                                 diff_percent=diff_percent,
                                                 iso_xlabel=xvar,
                                                 iso_ylabel=yvar,
                                                 iso_cmap=iso_cmap,
                                                 iso_levels=iso_levels)
end


# * Joint Equilibrium
if run_jeq
    iso_xlabel = :sigmah
    iso_ylabel = :mu_s
    iso_levels=15
    iso_cmap = "cool" # ModelPlots.iso_cmaps["pooling"]
    cat_cmap = "bone_r"
    k_otc = 32.5
    sup_title_sep=true

# ** Electronic Platforms
# *** Risky Type's Market to Book Ratio
    zvar = :mbr
    ft = :rt
    ft_zvar = ft == :st ? Symbol(:s_, zvar) : Symbol(:r_, zvar)
    eq_type = :ep_market
    zvar_diff = false
    diff_percent=false

    # s_otc_fv generated in the init file.
    # plot_title = ModelPlots.get_final_contour_plot_title(poold[:df], zvar, ft)
    plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                                        xvar=iso_xlabel, yvar=iso_ylabel,
                                                        diff=zvar_diff, diff_percent=diff_percent)


    # function will back out the inverse coupon ratio and the misrepresenting type from the dataframe.
    file_path_name = ModelPlots.get_nrmp_contour_plot_path_name(poold[:df], ft_zvar;
                                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                                eq_type=eq_type)

    ep_r_mbr_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[ft_zvar],
                                                                 pltd[:bool_Z];
                                                                 iso_xlabel=iso_xlabel,
                                                                 iso_ylabel=iso_ylabel,
                                                                 fig_title=plot_title,
                                                                 file_path_name=file_path_name,
                                                                 iso_cmap=iso_cmap,
                                                                 iso_levels=iso_levels,
                                                                 cat_cmap=cat_cmap)

# *** Safe Type's Firm Value
    zvar = :fv
    ft = :st
    ft_zvar = ft == :st ? Symbol(:s_, zvar) : Symbol(:r_, zvar)
    eq_type = :ep_market
    zvar_diff = false
    diff_percent=false

    # s_fi_fv generated in the init file.
    # plot_title = ModelPlots.get_final_contour_plot_title(poold[:df], zvar, ft)
    plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                                        xvar=iso_xlabel, yvar=iso_ylabel,
                                                        diff=zvar_diff, diff_percent=diff_percent,
                                                        s_fi_fv=s_fi_fv)

    # function will back out the inverse coupon ratio and the misrepresenting type from the dataframe.
    file_path_name = ModelPlots.get_nrmp_contour_plot_path_name(poold[:df], ft_zvar;
                                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                                eq_type=eq_type)

    ep_s_fv_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[ft_zvar],
                                                                pltd[:bool_Z];
                                                                sup_title_sep=sup_title_sep,
                                                                iso_xlabel=iso_xlabel,
                                                                iso_ylabel=iso_ylabel,
                                                                fig_title=plot_title,
                                                                file_path_name=file_path_name,
                                                                iso_cmap=iso_cmap,
                                                                iso_levels=iso_levels,
                                                                cat_cmap=cat_cmap)


# ** Dual Market Equilibria
    zvar = :fv
    ft = :st
    ft_zvar = ft == :st ? Symbol(:s_, zvar) : Symbol(:r_, zvar)
    eq_type = :dual_market
    zvar_diff = false
    diff_percent=false

    # s_fi_fv and s_otc_fv generated in the init file.

    # plot_title = ModelPlots.get_final_contour_plot_title(poold[:df], zvar, ft; k_otc=k_otc)
    plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                                        xvar=iso_xlabel, yvar=iso_ylabel,
                                                        diff=zvar_diff, diff_percent=diff_percent,
                                                        k_otc=k_otc,
                                                        s_fi_fv=s_fi_fv,
                                                        s_otc_fv=s_otc_fv)


    # function will back out the inverse coupon ratio and the misrepresenting type from the dataframe.
    file_path_name = ModelPlots.get_nrmp_contour_plot_path_name(poold[:df], ft_zvar;
                                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                                eq_type=eq_type)


    ep_s_fv_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[ft_zvar],
                                                                pltd[:bool_OTC_EP];
                                                                sup_title_sep=sup_title_sep,
                                                                iso_xlabel=iso_xlabel,
                                                                iso_ylabel=iso_ylabel,
                                                                fig_title=plot_title,
                                                                file_path_name=file_path_name,
                                                                iso_cmap=iso_cmap,
                                                                iso_levels=iso_levels,
                                                                cat_cmap=cat_cmap)
end
