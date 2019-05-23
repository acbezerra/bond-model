

# ATTENTION: kappa and iota values should be in basis points!
function rmp_fi_plotfun(xvar::Symbol, yvars::Array{Symbol,1},
                        sfdf::DataFrame,
                        rfdf::DataFrame;
                        xgrid::StepRangeLen{Float64,
                                            Base.TwicePrecision{Float64},
                                            Base.TwicePrecision{Float64}}=range(.0, stop=0., length=0),
                        fv_iota::Float64=NaN,
                        mbr_iota::Float64=NaN,
                        misrep_iota::Float64=NaN,
                        color_rm_region::Bool=false,
                        color_nrm_region::Bool=false,
                        color_conflict_region::Bool=false,
                        color_misrep_region::Bool=true,
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
        xgrid = range(minimum(sfdf[:iota]), stop=maximum(sfdf[:iota]), length=10^5)
    end
    
    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=figsize)
    fig = rmp_core_plot(fig, xvar, yvars,
                        sfdf, rfdf, xgrid,
                        subplots;
                        fv_iota=fv_iota,
                        mbr_iota=mbr_iota,
                        cvmlinestyles=cvmlinestyles,
                        cvmmarkers=cvmmarkers,
                        svmlinestyles=svmlinestyles,
                        svmmarkers=svmmarkers,
                        color_rm_region=color_rm_region,
                        color_nrm_region=color_nrm_region,
                        color_conflict_region=color_conflict_region)

    
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
    

    # if save_fig
    #     fig_folder, fig_name = rmp_plot_dirs(yvars, xvar;
    #                                          m=rfdf[1, :m],
    #                                          main_dir_path=main_dir_path,
    #                                          plots_dir=plots_dir,
    #                                          rmp_plots_dir=rmp_plots_dir)
    #     PyPlot.savefig(string(fig_folder, "/", fig_name), dpi=fig_dpi, bbox_inches="tight")
    # end

    return fig
end 
