

function rmp_subplotfun(fig::Figure, xvar::Symbol,
                        yvar::Symbol,
                        sfdf::DataFrame,
                        rfdf::DataFrame,
                        xgrid::StepRangeLen{Float64,
                                              Base.TwicePrecision{Float64},
                                              Base.TwicePrecision{Float64}};
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
                                xgrid=xgrid)
        else
            ax = plot_cvm_curve(ax, xvar, yvar, sfdf, xloc, yloc)
        end
    else
        # Plot SVM Curve
        xloc = .95 * rfdf[end, xvar]
        yloc = .1 * (maximum(rfdf[yvar]) - minimum(rfdf[yvar])) + minimum(rfdf[yvar])
        ax = plot_svm_curve(ax, xvar, yvar, rfdf, xloc, yloc, xgrid)
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
        ax = plot_svm_curve(ax, xvar, yvar, rfdf, xloc, yloc, xgrid)
        # ##############################################################
    else
        xloc = .95 * xmax 
        yloc = sfdf[1, yvar] 
        ax = plot_cvm_curve(ax, xvar, yvar, sfdf, xloc, yloc)
    end
    # ###################################################################

    # Plot Misrepresentation Curve #################################
    if .&(!isempty(misrepdf), (Symbol(:r_, yvar) in names(misrepdf)))
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
            xgrid = range(minimum(sfdf[xvar]), stop=maximum(sfdf[xvar]), length=10^5)
        else
            xgrid = range(minimum(rfdf[xvar]), stop=maximum(rfdf[xvar]), length=10^5)
        end
    end
    
    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=figsize)
    fig = rmp_core_plot(fig, xvar, yvars,
                        sfdf, rfdf, xgrid,
                        subplots;
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
    suptitle_params = join([string("\$", tlabels[x][1], "= \$ ",
                                   str_format_fun(ModelPlots.tlabels[x][2], rfdf[1, x]))
                         for x in rmp_plots_title_params_order], ", ")
    plot_suptitle = latexstring(suptitle_yvars, " under Full Information \n", 
                             " for ", suptitle_params)
    fig.suptitle(plot_suptitle, fontsize=14)

    
    if !isnan(figpad)
        fig.tight_layout(pad=figpad)
    end
    

    if save_fig
        fig_folder, fig_name = rmp_plot_dirs(yvars, xvar;
                                             m=rfdf[1, :m],
                                             main_dir_path=main_dir_path,
                                             plots_dir=plots_dir,
                                             rmp_plots_dir=rmp_plots_dir,
                                             misrep=!isempty(misrepdf))
        PyPlot.savefig(string(fig_folder, "/", fig_name), dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end 
