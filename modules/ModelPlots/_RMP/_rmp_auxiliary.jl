

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

    if all(df[yvar] .> yval)
        return Inf
    elseif all(df[yvar] .< yval)
        return -Inf
    else
        if size(xgrid, 1) == 0
            xgrid = range(minimum(df[xvar]), stop=maximum(df[xvar]), length=10^5)
        end
        
        yinterp = Dierckx.Spline1D(df[xvar], df[yvar];
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

    xvals = svmdf[xvar]
    if size(xgrid, 1) == 0
        xgrid = range(minimum(xvals), stop=maximum(xvals), length=10^5)
    end
    
    misrep_yval_interp = Dierckx.Spline1D(misrepdf[Symbol(:r_, xvar)], 
                                          misrepdf[Symbol(:r_, yvar)]; k=3, bc="extrapolate")
    if xvar != :sigmah
        fi_yvals = [maximum([x, cvmdf[1, yvar]]) for x in svmdf[yvar]]
        if xvar == :lambda
            fi_yvals = svmdf[yvar]
        end 
        fi_yval_interp = Dierckx.Spline1D(xvals, fi_yvals; k=3, bc="extrapolate")
        
        return xgrid[argmin(abs.(fi_yval_interp(xgrid) .- misrep_yval_interp(xgrid)))]
    else
        fi_yval_interp = Dierckx.Spline1D(xvals, svmdf[yvar]; k=3, bc="extrapolate")

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
                       m::Float64=NaN,
                       main_dir_path::String=main_dir_path,
                       plots_dir::String=plots_dir,
                       rmp_plots_dir::String=rmp_plots_dir,
                       misrep::Bool=false,
                       fname_ext::String=rmp_fname_ext)
    plots_path = string(main_dir_path, "/", plots_dir)
    fig_path = string(plots_path, "/", rmp_plots_dir)
    
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

    type_prefix = rmp_full_info_prefix
    if misrep
        type_prefix = rmp_misrep_prefix
    end
    
    fig_name = string(rmp_fn_prefix, "_", type_prefix, "_", yvars_fig_name, "_", xvar,".", fname_ext)
   
    return fig_folder, fig_name
end

