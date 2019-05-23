
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