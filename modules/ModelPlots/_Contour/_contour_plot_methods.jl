

function get_title_value(df, var)
    if !(var in names(df))
        var = Symbol(:r_, var)
    end
    
    return [parse(Float64, string(x)) for x in unique(df[var]) if !isnan(x)][1]
end


function get_final_contour_plot_title(df::DataFrame, zvar::Symbol,
                                      ft::Symbol;
                                      k_otc::Float64=NaN,
                                      params_list::Array{Symbol,1}=vcat(:mu_s,
                                                                        contour_plots_title_params_order))

    df[:pcr] = df[1, :p]/df[1, :c]
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
                                eqfuns::Dict{Symbol, Any}, 
                                zvar::Symbol;
                                ft::Symbol=Symbol(""),
                                diff_fun::Bool=false,
                                rm_prefix::Symbol=Symbol(""),
                                params_list::Array{Symbol,1}=contour_plots_title_params_order)
    
    eq_type = df[1, :eq_type]

    mu_s = NaN
    if eq_type in ["pooling", "separating"]
        mu_s =  df[1, :mu_s]
        params_list = vcat(:mu_s, params_list)
    end
    df[:pcr] = df[1, :p]/df[1, :c]

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
                         heat_cols=iso_plt_inputs[:heat_cols],
                         cat_Z=[],
                         cat_cmap="GnBu",
                         cat_alpha=.25)
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
                           cmap=cat_cmap, levels=cat_levels, alpha=cat_alpha)
        cb1 = fig.colorbar(CS1, ax=ax1, ticks=reverse(cats))#, orientation="horizontal")
        cb1.set_ticklabels(cat_tick_labels)
        cb1.set_clim(1, cat_levels + 1)
    end
    

    return fig, ax1, ax2
end


function plot_iso_contour_curves(fd::Dict{Symbol, Any},
                                 zfun;
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
                                 w_pad=iso_plt_inputs[:w_pad])
    
    X, Y, Z = form_mesh_grid(fd[:xvals], fd[:yvals], zfun)

    fig, ax1, ax2 = plot_iso_curves(X, Y, Z; 
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
                                            cat_cmap="GnBu")

    fig, ax1, ax2 = ModelPlots.plot_iso_curves(X, Y, Z;
                                               iso_cmap=iso_cmap, 
                                               iso_levels=15,
                                               cat_Z=eq_type_Z,
                                               cat_cmap=cat_cmap)
#    ax1.contourf(X, Y, eq_type_Z, cmap="GnBu_r", alpha=.4)
    
    if !isempty(fig_title)
        fig.suptitle(fig_title, fontsize=title_font_size)
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end

 

# ################################################################
# TRASH ##########################################################
# ################################################################
    # if eq_type in ["full_info", "misrep"]
    #     fun_candidates =  [x for x in keys(eqfuns) if
    #                        occursin(string(contour_zvars_sym[zvar]), string(x))]
    # else
    #     fun_candidates =  [x for x in keys(eqfuns[ft]) if
    #                    occursin(string(contour_zvars_sym[zvar]), string(x))]
    # end

        # zfun  = [x for x in fun_candidates if .&(occursin(string(eq_type_title[eq_type][1]),
        #                                                   string(x)),
        #                                          occursin("diff", string(x)))][1]

      # zfun  = [x for x in fun_candidates if .&(occursin(string(rm_prefix), string(x)),
        #                                          !occursin("diff", string(x)))][1]


# function iso_curves_title(df::DataFrame, zvar::Symbol;
#                           ft_prefix::Symbol=Symbol(""))

#     firm_type = "" 
#     if ft_prefix == :s_
#         firm_type = "Safe Type's "
#     else
#         firm_type = "Risky Type's "
#     end
    
#     title_params = join([string("\$", contour_tlabels[x][1], "= \$ ",
#                                 str_format_fun(contour_tlabels[x][2], 
#                                                get_title_value(df, x)))
#                          for x in contour_plots_title_params_order], ", ")
    
#     eq_type = df[1, :eq_type]
#     if eq_type == "misrep"
#         plot_title = latexstring(firm_type, contour_tlabels[zvar][1], " in case of Misrepresentation",
#                                  "\n for ", title_params)
#     else
#         plot_title = latexstring(firm_type, "Optimal ", contour_tlabels[zvar][1], " in a ",
#                                  eq_type_title[eq_type][2], " Equilibrium ",
#                                  "\n for ", title_params)
#     end   
    
    
#     return plot_title
# end


# function plot_fi_iso_curves(fidf::DataFrame;
#                             xvar::Symbol=:iota,
#                             yvar::Symbol=:sigmah, 
#                             zvars::Array{Symbol,1}=[:firm_value, :MBR],
#                             seaborn_style=iso_plt_inputs[:seaborn_style], 
#                             iso_levels=iso_plt_inputs[:iso_levels],
#                             heat_levels=iso_plt_inputs[:heat_levels],
#                             iso_cmap=iso_cmaps["full_info"],
#                             heat_cmap::String="",
#                             fig_aspect=iso_plt_inputs[:fig_aspect],
#                             iso_fontsize=iso_plt_inputs[:iso_fontsize],
#                             use_subgrid=iso_plt_inputs[:use_subgrid],
#                             subgrid_rows=iso_plt_inputs[:subgrid_rows],
#                             iso_cols=iso_plt_inputs[:iso_cols],
#                             heat_cols=iso_plt_inputs[:heat_cols],
#                             title_font_size=iso_plt_inputs[:title_font_size],
#                             tight_pad=iso_plt_inputs[:tight_pad],
#                             h_pad=iso_plt_inputs[:h_pad],
#                             w_pad=iso_plt_inputs[:w_pad])
    
#     cond = (fidf[:iota] .!= 2.5)
#     fi_fd = interp_z_values(fidf[cond, :], xvar, yvar, zvars)

#     # fi_rm_cond_fun, fi_fv_fun, fi_mbr_fun
#     fi_funs = fi_payoff_functions(fi_fd; xvar=xvar, yvar=yvar)

#     fi_X, fi_Y, fi_FV = form_mesh_grid(fi_fd[:xvals], 
#                                        fi_fd[:yvals],
#                                        fi_funs[:fv])

#     fi_fig, ax1, ax2 = plot_iso_curves(fi_X, fi_Y, fi_FV; 
#                                        seaborn_style=seaborn_style, 
#                                        iso_levels=iso_levels,
#                                        heat_levels=heat_levels, 
#                                        iso_cmap=iso_cmap,
#                                        heat_cmap=heat_cmap,
#                                        fig_aspect=fig_aspect,
#                                        iso_fontsize=iso_fontsize,
#                                        use_subgrid=use_subgrid,
#                                        subgrid_rows=subgrid_rows,
#                                        iso_cols=iso_cols,
#                                        heat_cols=heat_cols)

#     fi_title = iso_curves_title(fidf, :firm_value)
#     fi_fig.suptitle(fi_title, fontsize=title_font_size)
#     PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

#     return fi_fig
# end


# function plot_misrep_iso_curve(mp_fd::Dict{Symbol, Any}, mp_funs::Dict{Symbol, Any},
#     #fidf::DataFrame, misrepdf::DataFrame; 
#                                xvar::Symbol=:iota, 
#                                yvar::Symbol=:sigmah,
#                                zvars::Array{Symbol,1}=[:MBR],
#                                seaborn_style=iso_plt_inputs[:seaborn_style], 
#                                iso_levels=iso_plt_inputs[:iso_levels],
#                                heat_levels=iso_plt_inputs[:heat_levels],
#                                iso_cmap=iso_cmaps["misrep"],
#                                heat_cmap::String="",
#                                fig_aspect=iso_plt_inputs[:fig_aspect],
#                                iso_fontsize=iso_plt_inputs[:iso_fontsize],
#                                use_subgrid=iso_plt_inputs[:use_subgrid],
#                                subgrid_rows=iso_plt_inputs[:subgrid_rows],
#                                iso_cols=iso_plt_inputs[:iso_cols],
#                                heat_cols=iso_plt_inputs[:heat_cols],
#                                title_font_size=iso_plt_inputs[:title_font_size],
#                                tight_pad=iso_plt_inputs[:tight_pad],
#                                h_pad=iso_plt_inputs[:h_pad],
#                                w_pad=iso_plt_inputs[:w_pad])

#     # fi_fd = interp_z_values(fidf, xvar, yvar, unique(vcat(:firm_value, zvars)))
#     # mp_fd = interp_z_values(misrepdf, xvar, yvar,
#     #                         [Symbol(:r_, zvar) for zvar in zvars];
#     #                         ft_xy=:r_) 

#     # fi_funs, mp_funs = misrep_payoff_functions(fi_fd, mp_fd; xvar=xvar, yvar=yvar)

#     mp_X, mp_Y, mp_XY_diff_MBR = form_mesh_grid(mp_fd[:xvals], mp_fd[:yvals],
#                                                 mp_funs[:mp_fi_mbr_diff])

#     mp_fig, ax1, ax2 = plot_iso_curves(mp_X, mp_Y, mp_XY_diff_MBR;
#                                        xvar=mp_funs[:xvar], yvar=mp_funs[:yvar],
#                                        seaborn_style=seaborn_style, 
#                                        iso_levels=iso_levels, 
#                                        heat_levels=heat_levels, 
#                                        iso_cmap=iso_cmap, 
#                                        heat_cmap=heat_cmap,
#                                        fig_aspect=fig_aspect, 
#                                        iso_fontsize=iso_fontsize,
#                                        use_subgrid=use_subgrid,
#                                        subgrid_rows=subgrid_rows,
#                                        iso_cols=iso_cols,
#                                        heat_cols=heat_cols)
    
#     # Set Title and Padding
#     # mp_title = iso_curves_title(misrepdf, :MBR; ft_prefix=:r_)
#     # mp_fig.suptitle(mp_title, fontsize=title_font_size)
#     PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

#     return mp_fig
# end



# function plot_pool_iso_curve(fidf::DataFrame, pooldf::DataFrame, firm_type::Symbol,
#                              plt_fun::Symbol; pool_fd::Dict=Dict(),
#                              xvar::Symbol=:iota, yvar::Symbol=:sigmah,
#                              zvars::Array{Symbol, 1}=[:r_MBR, :s_firm_value],
#                              seaborn_style=iso_plt_inputs[:seaborn_style], 
#                              iso_levels=iso_plt_inputs[:iso_levels],
#                              heat_levels=iso_plt_inputs[:heat_levels],
#                              iso_cmap=iso_cmaps["pooling"],
#                              heat_cmap::String="",
#                              fig_aspect=iso_plt_inputs[:fig_aspect],
#                              iso_fontsize=iso_plt_inputs[:iso_fontsize],
#                              use_subgrid=iso_plt_inputs[:use_subgrid],
#                              subgrid_rows=iso_plt_inputs[:subgrid_rows],
#                              iso_cols=iso_plt_inputs[:iso_cols],
#                              heat_cols=iso_plt_inputs[:heat_cols],
#                              title_font_size=iso_plt_inputs[:title_font_size],
#                              tight_pad=iso_plt_inputs[:tight_pad],
#                              h_pad=iso_plt_inputs[:h_pad],
#                              w_pad=iso_plt_inputs[:w_pad])

    
#     fi_zvars = unique(vcat(:firm_value,
#                            [Symbol(strip(string(z), ['r', '_', 's'])) for z in zvars]))
#     fi_fd = interp_z_values(fidf, xvar, yvar, fi_zvars)
#     if isempty(pool_fd)
#         pool_fd = interp_z_values(pooldf, xvar, yvar, 
#                                   zvars; ft_xy=:r_) 
#     end
    
#     # r_rm_cond_fun, r_pool_mbr_fun, s_pool_fv_fun
#     pool_funs = joint_eq_payoff_functions(fi_fd, pool_fd; eq_type="pooling",
#                                           xvar=xvar, yvar=yvar)

    
#     pool_X, pool_Y, pool_mbr = form_mesh_grid(pool_fd[:xvals], 
#                                               pool_fd[:yvals],
#                                               pool_funs[firm_type][plt_fun])
    
#     pool_fig, ax1, ax2 = plot_iso_curves(pool_X, pool_Y, pool_mbr; 
#                                          seaborn_style=seaborn_style, 
#                                          iso_levels=iso_levels, 
#                                          heat_levels=heat_levels, 
#                                          iso_cmap=iso_cmap, 
#                                          heat_cmap=heat_cmap,
#                                          fig_aspect=fig_aspect, 
#                                          iso_fontsize=iso_fontsize,
#                                          use_subgrid=use_subgrid,
#                                          subgrid_rows=subgrid_rows,
#                                          iso_cols=iso_cols,
#                                          heat_cols=heat_cols)

#     # Set Title and Padding
#     pool_title = iso_curves_title(pooldf, :MBR; ft_prefix=:r_)
#     pool_fig.suptitle(pool_title, fontsize=title_font_size)
#     PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

#     return pool_fig
# end
