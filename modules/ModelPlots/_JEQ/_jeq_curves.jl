


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

