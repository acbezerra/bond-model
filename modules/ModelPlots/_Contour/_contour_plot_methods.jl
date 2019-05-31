
function get_title_value(df, var)
    if !(var in names(df))
        var = Symbol(:r_, var)
    end
    
    return [parse(Float64, string(x)) for x in unique(df[var]) if !isnan(x)][1]
end


function iso_curves_title(df::DataFrame, zvar::Symbol;
                          ft_prefix::Symbol=Symbol(""))

    firm_type = "" 
    if ft_prefix == :s_
        firm_type = "Safe Type's "
    else
        firm_type = "Risky Type's "
    end
    
    title_params = join([string("\$", contour_tlabels[x][1], "= \$ ",
                                str_format_fun(contour_tlabels[x][2], 
                                               get_title_value(df, x)))
                         for x in contour_plots_title_params_order], ", ")
    
    eq_type = df[1, :eq_type]
    if eq_type == "misrep"
        plot_title = latexstring(firm_type, contour_tlabels[zvar][1], " in case of Misrepresentation",
                                 "\n for ", title_params)
    else
        plot_title = latexstring(firm_type, "Optimal ", contour_tlabels[zvar][1], " in a ",
                                 eq_type_title[eq_type], " Equilibrium ",
                                 "\n for ", title_params)
    end   
    
    
    return plot_title
end


function plot_iso_curves(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2};
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
                         heat_cols=iso_plt_inputs[:heat_cols])
    
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
        ax1 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols), (0, 0),
                                  rowspan=subgrid_rows, colspan=iso_cols)
        ax2 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols), (0, iso_cols),
                                  rowspan=subgrid_rows, colspan=heat_cols)
    else
        fig, axs = PyPlot.subplots(1, 2, figsize=(w, h), sharey=true)
        ax1 = axs[1] # fig.add_subplot(121)
        ax2 = axs[2] # fig.add_subplot(122)
    end
    # ######################################################################
    
    CS = ax1.contour(X, Y, Z, levels=iso_levels, cmap=iso_cmap)
    ax1.clabel(CS, inline=5, fontsize=iso_fontsize)
    ax1.set_ylabel(latexstring("\$", xylabels[:sigmah][1], "\$"), labelpad=10)
    ax1.set_xlabel(latexstring("\$", xylabels[:iota][1], "\$"), labelpad=10)

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
    ax2.set_xlabel(latexstring("\$", xylabels[:iota][1], "\$"), labelpad=10)
    
    # Add Colorbar
    cbar = fig.colorbar(CS2)

    return fig, ax1, ax2
end


function plot_fi_iso_curves(fidf::DataFrame; xvar::Symbol=:iota,
                            yvar::Symbol=:sigmah, 
                            zvars::Array{Symbol,1}=[:firm_value, :MBR],
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
                            tight_pad=iso_plt_inputs[:tight_pad],
                            h_pad=iso_plt_inputs[:h_pad],
                            w_pad=iso_plt_inputs[:w_pad])
    
    cond = (fidf[:iota] .!= 2.5)
    fi_fd = interp_z_values(fidf[cond, :], xvar, yvar, zvars)

    fi_rm_cond_fun, fi_fv_fun, fi_mbr_fun = fi_payoff_functions(fi_fd; 
                                                                xvar=xvar, 
                                                                yvar=yvar)

    fi_X, fi_Y, fi_FV = form_mesh_grid(fi_fd[:xvals], 
                                       fi_fd[:yvals],
                                       fi_fv_fun)

    fi_fig, ax1, ax2 = plot_iso_curves(fi_X, fi_Y, fi_FV; 
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
                                       heat_cols=heat_cols)

    fi_title = iso_curves_title(fidf, :firm_value)
    fi_fig.suptitle(fi_title, fontsize=title_font_size)
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    return fi_fig
end


function plot_misrep_iso_curve(fidf::DataFrame, misrepdf::DataFrame; 
                               xvar::Symbol=:iota, 
                               yvar::Symbol=:sigmah,
                               zvars::Array{Symbol,1}=[:MBR],
                               seaborn_style=iso_plt_inputs[:seaborn_style], 
                               iso_levels=iso_plt_inputs[:iso_levels],
                               heat_levels=iso_plt_inputs[:heat_levels],
                               iso_cmap=iso_cmaps["misrep"],
                               heat_cmap::String="",
                               fig_aspect=iso_plt_inputs[:fig_aspect],
                               iso_fontsize=iso_plt_inputs[:iso_fontsize],
                               use_subgrid=iso_plt_inputs[:use_subgrid],
                               subgrid_rows=iso_plt_inputs[:subgrid_rows],
                               iso_cols=iso_plt_inputs[:iso_cols],
                               heat_cols=iso_plt_inputs[:heat_cols],
                               title_font_size=iso_plt_inputs[:title_font_size],
                               tight_pad=iso_plt_inputs[:tight_pad],
                               h_pad=iso_plt_inputs[:h_pad],
                               w_pad=iso_plt_inputs[:w_pad])

    fi_fd = interp_z_values(fidf, xvar, yvar, vcat(:firm_value, zvars))
    mp_fd = interp_z_values(misrepdf, xvar, yvar,
                            [Symbol(:r_, zvar) for zvar in zvars];
                            ft_xy=:r_) 

    r_fi_mbr_fun, r_mp_mbr_fun, r_mbr_diff_fun = misrep_payoff_functions(fi_fd, mp_fd; 
                                                                         xvar=xvar, 
                                                                         yvar=yvar)

    mp_X, mp_Y, mp_XY_diff_MBR = form_mesh_grid(mp_fd[:xvals], mp_fd[:yvals],
                                                r_mbr_diff_fun)

    mp_fig, ax1, ax2 = plot_iso_curves(mp_X, mp_Y, mp_XY_diff_MBR; 
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
                                       heat_cols=heat_cols)
    
    # Set Title and Padding
    mp_title = iso_curves_title(misrepdf, :MBR; ft_prefix=:r_)
    mp_fig.suptitle(mp_title, fontsize=title_font_size)
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    return mp_fig
end



function plot_pool_iso_curve(pooldf::DataFrame; pool_fd::Dict=Dict(),
                             xvar::Symbol=:iota, yvar::Symbol=:sigmah,
                             zvars::Array{Symbol, 1}=[:r_MBR, :s_firm_value],
                             seaborn_style=iso_plt_inputs[:seaborn_style], 
                             iso_levels=iso_plt_inputs[:iso_levels],
                             heat_levels=iso_plt_inputs[:heat_levels],
                             iso_cmap=iso_cmaps["pooling"],
                             heat_cmap::String="",
                             fig_aspect=iso_plt_inputs[:fig_aspect],
                             iso_fontsize=iso_plt_inputs[:iso_fontsize],
                             use_subgrid=iso_plt_inputs[:use_subgrid],
                             subgrid_rows=iso_plt_inputs[:subgrid_rows],
                             iso_cols=iso_plt_inputs[:iso_cols],
                             heat_cols=iso_plt_inputs[:heat_cols],
                             title_font_size=iso_plt_inputs[:title_font_size],
                             tight_pad=iso_plt_inputs[:tight_pad],
                             h_pad=iso_plt_inputs[:h_pad],
                             w_pad=iso_plt_inputs[:w_pad])
    
    if isempty(pool_fd)
        pool_fd = ModelPlots.interp_z_values(pooldf, xvar, yvar, 
                                             zvars; ft_xy=:r_) 
    end
    
    r_rm_cond_fun, r_pool_mbr_fun, s_pool_fv_fun = joint_eq_payoff_functions(pool_fd;
                                                                             eq_type="pooling",
                                                                             xvar=xvar, 
                                                                             yvar=yvar)
    
    pool_X, pool_Y, pool_mbr = form_mesh_grid(pool_fd[:xvals], 
                                              pool_fd[:yvals],
                                              r_pool_mbr_fun)
    
    pool_fig, ax1, ax2 = plot_iso_curves(pool_X, pool_Y, pool_mbr; 
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
                                         heat_cols=heat_cols)

    # Set Title and Padding
    pool_title = iso_curves_title(pooldf, :MBR; ft_prefix=:r_)
    pool_fig.suptitle(pool_title, fontsize=title_font_size)
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    return pool_fig
end
