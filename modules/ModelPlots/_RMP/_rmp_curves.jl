

function plot_cvm_curve(ax, xvar::Symbol, yvar::Symbol,
                        sfdf::DataFrame,
                        text_xloc::Float64, text_yloc::Float64;
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
        ax.plot(sfdf[xvar], sfdf[yvar];
                color=cvm_curve_color, 
                linewidth=1,
                linestyle=cvmlinestyles[1],
                marker=cvmmarkers[1], 
                markersize=3)
    else
        y_interp = Dierckx.Spline1D(sfdf[xvar], sfdf[yvar]; k=3, bc="extrapolate")
        ax.plot(xgrid, y_interp(xgrid);
                color=cvm_curve_color, 
                linewidth=1,
                linestyle=cvmlinestyles[1])
    end
    

    # Add Legend to the Curve
    if xvar == :iota
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
                                              Base.TwicePrecision{Float64}})

    if size(rfdf, 1) == 1   
        ax.axhline(rfdf[1, yvar], 
                   color=svm_curve_color,
                   linewidth=1, 
                   linestyle=svmlinestyles[1], 
                   marker=svmmarkers[1])
    elseif size(xgrid, 1) == 0
        ax.plot(rfdf[xvar], rfdf[yvar];
                color=svm_curve_color, 
                linewidth=1,
                linestyle=svmlinestyles[1],
                marker=svmmarkers[1], 
                markersize=3)
    else
        y_interp = Dierckx.Spline1D(rfdf[xvar], rfdf[yvar]; k=3, bc="extrapolate")
        ax.plot(xgrid, y_interp(xgrid);
                color=svm_curve_color, 
                linewidth=1,
                linestyle=svmlinestyles[1])
    end


    # Add Legend to the Curve
    if xvar == :iota
        svm_label = latexstring("(\$", 
                                ModelPlots.cvs_xlabels[:iota][1], ", ",
                                ModelPlots.cvs_xlabels[:lambda][1], ", ",
                                ModelPlots.cvs_xlabels[:sigmah][1],
                                "\$) = (",
                                rfdf[1, :iota], ", ",
                                rfdf[1, :lambda], ", ", 
                                rfdf[1, :sigmah], ")")
    elseif xvar == :sigmah
        svm_label = latexstring("(\$", 
                                ModelPlots.cvs_xlabels[:iota][1], ", ",
                                ModelPlots.cvs_xlabels[:lambda][1], 
                                "\$) = (",
                                rfdf[1, :iota], ", ",
                                rfdf[1, :lambda], ")")
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
        cols = [vcat([misrepdf[1, Symbol(:s_, xvar)], misrepdf[Symbol(:r_, xvar)]...]), 
                vcat([misrepdf[1, Symbol(:s_, yvar)], misrepdf[Symbol(:r_, yvar)]...])]
        
        tmp = sort!(DataFrame(Dict(zip([xvar, yvar], cols))), xvar)
        if size(xgrid, 1) > 0
            y_interp = Dierckx.Spline1D(tmp[xvar], tmp[yvar]; k=3, bc="extrapolate")
            ax.plot(xgrid, y_interp(xgrid),
                    color=misrep_curve_color, 
                    linewidth=1,
                    linestyle=svmlinestyles[1],
                    marker=svmmarkers[1], 
                    markersize=3)
        else
            ax.plot(tmp[xvar], tmp[yvar];
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

