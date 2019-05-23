
function form_cvm_svm_combinations(cvm_m_comb_nums::Array{Int64,1},
                           svm_m_comb_nums::Array{Int64,1};
                           m::Float64=1.)
    id_vars = [:comb_num, :m_comb_num]

    sbt = get_bt(; model="svm", m=m, m_comb_num=1)
    cbt = get_bt(; model="cvm", m=m, m_comb_num=1)
    
    cvmls = [ ]
    for i in cvm_m_comb_nums
        tmp = cbt.bp.df[.&(cbt.bp.df[:m].==m, cbt.bp.df[:m_comb_num].== i), :]
        cols = [x for x in names(tmp) if !(x in vcat(id_vars, [:lambda, :sigmah]))]
        push!(cvmls, Dict{Symbol, Float64}(cols .=> [tmp[1, var] for var in cols] ))
    end

    svmls = [ ]
    for i in svm_m_comb_nums
        tmp = sbt.bp.df[.&(sbt.bp.df[:m].==m, sbt.bp.df[:m_comb_num].== i), :]
        cols = [x for x in names(tmp) if !(x in vcat(id_vars, :sigmah))]
        push!(svmls, Dict{Symbol, Float64}(cols .=> [tmp[1, var] for var in cols] ))
    end

    return Dict("svm" => svmls, 
                "cvm" => cvmls)
end


function curve_legend(df_slice, var2::Symbol; model::String="svm",
                      xlabels::Dict{Symbol, Array{String,1}}=cvs_xlabels)

    var1_symbol = (model == "svm") ? Symbol(:kappa, :_ep) : Symbol(:kappa, :_otc)
    var2_val = (var2 == :iota) ? df_slice[1, var2] * 1e4 : df_slice[1, var2]
    
    return latexstring("\$\\left(", xlabels[var1_symbol][1], ", ",
                       xlabels[var2][1], "\\right) = \$ (",
                       df_slice[1, :kappa] * 1e4,
                       ", ", var2_val, ")")
end


function cvm_vs_svm_plot_dirs(yvar::Symbol, dvar::Symbol, xvar::Symbol;
                              m::Float64=NaN,
                              main_dir_path::String=main_dir_path,
                              plots_dir::String=plots_dir,
                              cvm_vs_svm_plots_dir::String=cvm_vs_svm_plots_dir)
    plots_path = string(main_dir_path, "/", plots_dir)
    fig_path = string(plots_path, "/", cvm_vs_svm_plots_dir)

    
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
    fig_name = string("cvm_vs_svm_", yvar, "_", dvar, "_", xvar,".png")
   
    return fig_folder, fig_name
end


function cvm_vs_svm_plotfun(cvmdf, svmdf,
                            xvar::Symbol, yvar::Symbol, dvar::Symbol;
                            figaspect::Float64=.5,
                            figpad::AbstractFloat=1.8, 
                            plot_vlines::Bool=false, 
                            xvars::Array{Symbol,1}=cvs_xvars,
                            cvmlinestyles::Array{String,1}=cvmlinestyles,
                            cvmmarkers::Array{String,1}=cvmmarkers,
                            svmlinestyles::Array{String,1}=svmlinestyles,
                            svmmarkers::Array{String,1}=svmmarkers,
                            save_fig::Bool=true,
                            fig_dpi::Int64=300,
                            main_dir_path::String=main_dir_path,
                            plots_dir::String=plots_dir,
                            cvm_vs_svm_plots_dir::String=cvm_vs_svm_plots_dir)


    # cvmdf = pt.cvm_data[[(x in cvm_combs) for x in pt.cvm_data[:comb_num]], :]
    # svmdf = pt.svm_data[[(x in svm_combs) for x in pt.svm_data[:comb_num]], :]

    # x-axis variable:
    # xvar =  [x for x in cvs_xvars if size(unique(svmdf[x]), 1)  == size(unique(pt.svm_data[x]), 1)][1]

    
    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=Tuple(PyPlot.figaspect(figaspect)))
    ax = fig.add_subplot(111)
    
    # Plot SVM Curves #############################################
    pchipflist = []
    xpos=0.
    i = 1
    for dval in unique(svmdf[dvar])
        # Slice SVM DataFrame
        svm_slice = svmdf[abs.(svmdf[dvar] .- dval) .< 1e-6, :]

        # Plot Curve
        ax.plot(svm_slice[xvar], svm_slice[yvar],
                color="blue", 
                linewidth=1,
                linestyle=svmlinestyles[i],
                marker=svmmarkers[i], 
                markersize=3)

        # Add Legend to the Curves
        svm_label = curve_legend(svm_slice, :lambda; model="svm")     
        ax.text(svm_slice[end, xvar], svm_slice[end, yvar],
                  svm_label, fontsize=10, va="bottom")


        # Interpolate SVM Curves and Store Them
        pchipf = Dierckx.Spline1D(svm_slice[xvar], svm_slice[yvar];
                                  k=3, bc="extrapolate")
        push!(pchipflist, pchipf)


        if dval == unique(svmdf[dvar])[end] 
            xvals = range(minimum(svm_slice[xvar]),
                          stop=maximum(svm_slice[xvar]), length=10^4)
            xpos = svm_slice[end, xvar]
        end

        i += 1
    end 
    # #############################################################
    

    # Plot CVM Curves #############################################
    for dval in unique(cvmdf[dvar])
        # Slice CVM DataFrame
        if !isnan(dval)
            cvm_slice = cvmdf[abs.(cvmdf[dvar] .- dval) .< 1e-6, :]
        else
            cvm_slice = cvmdf
        end
  
        # Plot Curve
        ax.axhline(cvm_slice[1, yvar], 
                   color="green",
                   linewidth=1, 
                   linestyle=cvmlinestyles[i], 
                   marker=cvmmarkers[i])

        # Add Legend to the Curves
        cvm_label = curve_legend(cvm_slice, :iota; model="cvm") 
        ax.text(xpos, 
                cvm_slice[1, yvar],
                cvm_label, fontsize=10, va="bottom")


        # Plot Vertical Lines
        if plot_vlines
            for pchipf in pchipflist
                sigstar = xvals[argmin([abs(pchipf(x) - cvm_slice[1, yvar]) for x in xvals])]
                ax.axvline(sigstar, color="black", linewidth=.5, linestyle="-.")
            end
        end        
    end
    # #############################################################

    
    # Label Axes 
    ax.set_xlabel(string("\$", cvs_xlabels[xvar][1], "\$"), labelpad=10)
    ax.set_ylabel(cvs_ylabels[yvar], labelpad=10)
    
    # Set Title
    title_params = join([string("\$", tlabels[x][1], "= \$ ",
                                str_format_fun(tlabels[x][2], svmdf[1, x]))
                         for x in cvs_plots_title_params_order], ", ")
    plot_title = latexstring("Optimal RMP-Conditional ", cvs_ylabels[yvar], " for ", title_params)
    fig.suptitle(plot_title, fontsize=14)
    ax.set_title("(\$\\kappa\$ and \$\\iota\$ values in b.p.)", fontsize=12)

    if !isnan(figpad)
        fig.tight_layout(pad=figpad)
    end

    if save_fig
        fig_folder, fig_name = cvm_vs_svm_plot_dirs(yvar, dvar, xvar;
                                                    m=svmdf[1, :m],
                                                    main_dir_path=main_dir_path,
                                                    plots_dir=plots_dir,
                                                    cvm_vs_svm_plots_dir=cvm_vs_svm_plots_dir)
        plt.savefig(string(fig_folder, "/", fig_name), dpi=fig_dpi, bbox_inches="tight")
    end
#     display(fig)
    return fig
end




# ##############################################################################
# ################################### TRASH ####################################
# ##############################################################################
function df_slicer(pt, i::Int64, combs; model::String="svm")
    tmp = combs[i]
    if model == "svm"
        svmloc = sum([abs.(pt.svm_data[x] .- tmp[x]) .< 1e-4 
                      for x in keys(tmp)]) .== length(keys(tmp))
        return pt.svm_data[svmloc, :]
        # sort!(pt.svm_data[svmloc, :], xvar)
    else
        cvmloc = sum([abs.(pt.cvm_data[x] .- tmp[x]) .< 1e-4 
                      for x in keys(tmp)]) .== length(keys(tmp))
        return pt.cvm_data[cvmloc, :]
    end
end


# function cvm_vs_svm_plot_path_fname(; main_dir::String=main_dir,)
#     # Paths to Plot ############################################
#     # Main Directory
#     main_dir_path = form_main_dir_path(main_dir)

    
    
# function cvm_vs_svm_plotfun(pt, yvar::Symbol, cvm_combs, svm_combs, dvar;
#                             figaspect::Float64=.5,
#                             figpad::AbstractFloat=1.8, 
#                             plot_vlines::Bool=false, 
#                             figPath::String="",
#                             xvars::Array{Symbol,1}=cvs_xvars,
#                             cvmlinestyles::Array{String,1}=cvmlinestyles,
#                             cvmmarkers::Array{String,1}=cvmmarkers,
#                             svmlinestyles::Array{String,1}=svmlinestyles,
#                             svmmarkers::Array{String,1}=svmmarkers)

#                             #pardict::Dict{String,Array{Any,1}};
    
#     cvmcombs = pardict["cvm"]
#     svmcombs = pardict["svm"]

#     xvar = [x for x in xvars if !(x in keys(svm_combs[1]))][1]
        
#     Seaborn.set(style="darkgrid")
#     fig = PyPlot.figure(figsize=Tuple(PyPlot.figaspect(figaspect)))
#     ax = fig.add_subplot(111)
    
#     # Plot SVM Curves
#     pchipflist = []
#     xpos=0.
#     i = 1
#     for i in 1:size(svmcombs, 1)
#         # Slice SVM DataFrame
#         svm_slice = df_slicer(pt, i, svmcombs; model="svm")

#         # Plot Curve
#         ax.plot(svm_slice[xvar], svm_slice[yvar],
#                 color="blue", 
#                 linewidth=1,
#                 linestyle=svmlinestyles[i],
#                 marker=svmmarkers[i], 
#                 markersize=3)

#         # Add Legend to the Curves
#         svm_label = curve_legend(svm_slice, :lambda; model="svm")     
#         ax.text(svm_slice[end, xvar], svm_slice[end, yvar],
#                   svm_label, fontsize=10, va="bottom")


#         # Interpolate SVM Curves and Store Them
#         pchipf = Dierckx.Spline1D(svm_slice[xvar], svm_slice[yvar];
#                                   k=3, bc="extrapolate")
#         push!(pchipflist, pchipf)
        
#        if i == 1
#            xvals = range(minimum(svm_slice[xvar]),
#                          stop=maximum(svm_slice[xvar]), length=10^4)
#            # ax.plot(xvals, pchipf(xvals), color='red')
#        end

#         if i == size(svmcombs,1)
#             xpos = svm_slice[end, xvar]
#         end

#         i += 1
#     end 
            
#     # Plot CVM Curves
#     for i in 1:size(cvmcombs, 1)
#         # Slice CVM DataFrame
#          cvm_slice = df_slicer(pt, i, cvmcombs; model="cvm")
#         tmp = cvmcombs[i]

#         # Plot Curve
#         ax.axhline(cvm_slice[1, yvar], 
#                    color="green",
#                    linewidth=1, 
#                    linestyle=cvmlinestyles[i], 
#                    marker=cvmmarkers[i])

#         # Add Legend to the Curves
#         cvm_label = curve_legend(cvm_slice, :iota; model="cvm") 
#         ax.text(xpos, 
#                 cvm_slice[1, yvar],
#                 cvm_label, fontsize=10, va="bottom")


#         # Plot Vertical Lines
#         if plot_vlines
#             for pchipf in pchipflist
#                 sigstar = xvals[argmin([abs(pchipf(x) - cvm_slice[1, yvar]) for x in xvals])]
#                 ax.axvline(sigstar, color="black", linewidth=.5, linestyle="-.")
#             end
#         end        
#     end

#     # Label Axes 
#     ax.set_xlabel(string("\$", cvs_xlabels[xvar][1], "\$"), labelpad=10)
#     ax.set_ylabel(cvs_ylabels[yvar], labelpad=10)
    
#     # Set Title
#     title_params = join([string("\$", tlabels[x][1], "= \$ ",
#                                 str_format_fun(tlabels[x][2], svmcombs[1][x]))
#                          for x in cvs_plots_title_params_order], ", ")
#     plot_title = latexstring("Optimal RMP-Conditional ", cvs_ylabels[yvar], " for ", title_params)
#     fig.suptitle(plot_title, fontsize=14)
#     ax.set_title("(\$\\kappa\$ and \$\\iota\$ values in b.p.)", fontsize=12)

#     if !isnan(figpad)
#         fig.tight_layout(pad=figpad)
#     end

#     if !isempty(figPath)
#         figFolder = string(figPath, "/m_", convert(Integer, svmcombs[1][:m]))
#         # figName = string("cvm_vs_svm_", cvs_ylabels[yvar], "__", 
#         #                  join([string(x, '_', svmcombs[1][x])
#         #                        for x in vcat(fixed_vars, yvars) if 
#         #                        !(x in [:sigmah, :lambda, :kappa, :iota])], "__"),
#         #                 ".png")

#         # plt.savefig(string(figFolder, "/", figName), dpi=300, bbox_inches="tight")
#     end
# #     display(fig)
#     return fig
# end

