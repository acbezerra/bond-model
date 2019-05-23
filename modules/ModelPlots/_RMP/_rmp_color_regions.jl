

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


function color_misrep_region_fun(ax, misrep_xvar::Float64,
                                 xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                     Base.TwicePrecision{Float64}},
                                 text_xloc::Float64, 
                                 text_yloc::Float64)
    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .>= misrep_xvar,
                    facecolor=misrep_region_color, alpha=0.25)
    
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
                           misrep_xvar::Float64=NaN,
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
    if .&(!isnan(misrep_xvar), color_misrep_region)
        xloc = .5 * misrep_xvar + .5 * xmax
        yloc = .1 * ymin + .9 * ymax
        ax = color_misrep_region_fun(ax, misrep_xvar, xgrid, xloc, yloc) 
    end

    return ax
end
