


function color_otc_region_fun(ax, fv_xvar::Float64,
                              xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                  Base.TwicePrecision{Float64}},
                              text_xloc::Float64, 
                              text_yloc::Float64)
    trans = PyPlot.matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.fill_between(xgrid, 0, 1, transform=trans,
                    where=xgrid .<= fv_xvar,
                    facecolor=otc_region_color, alpha=0.25)

    if !isinf(text_xloc)
        ax.text(text_xloc, text_yloc,
                "OTC Trading \n is optimal",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=11,
                style="italic",
                bbox=Dict("facecolor" => box_color, "alpha" => 0.5, "pad" => 10))
    end
    
    return ax
end

