
function cvm_data_handler(pt)
    kappa_vals = unique(pt.cvm_data[:kappa])
    sigmal_vals = unique(pt.cvm_data[:sigmal])
    iota_vals = [ i*10^4 for i in unique(pt.cvm_data[:iota]) if i < maximum(pt.cvm_data[:iota])]
    # scalarMap = Seaborn.color_palette("cool", n_colors=size(sigma_vals, 1))
    scalarMap = Seaborn.get_cmap("cool", size(sigmal_vals, 1))
    l_styles = ["-" , ":", "--", "-."]

    return kappa_vals, sigmal_vals, iota_vals, scalarMap, l_styles
end


function cvm_single_axis_plot(pt, fig, plot_sub_num::Int64, varname::Symbol;
                              kappa::Float64=NaN,
                              sigma_low::Float64=NaN,
                              legend::Bool=true)
    (kappa_vals, sigmal_vals,
     iota_vals, scalarMap, l_styles) = cvm_data_handler(pt)

    if isnan(kappa)
        kappa = kappa_vals[1]
    end
    if isnan(sigma_low)
        sigma_low = minimum(sigmal_vals)
    end

    ax = fig.add_subplot(plot_sub_num)
    for j in 1:size(sigmal_vals, 1)
        colorVal = scalarMap(j)
        pos = .&(abs.(pt.cvm_data[:kappa] .- kappa) .< 1e-6,
                 abs.(pt.cvm_data[:sigmal] .- sigmal_vals[j]) .< 1e-6)

        if (abs.(sigmal_vals[j] .- sigma_low) .< 1e-6)
            pos = .&(pos, [(x .* 10^4 in iota_vals) for x in pt.cvm_data[:iota]])
                     
            ax.plot(pt.cvm_data[pos, :iota] .* 10^4,
                    pt.cvm_data[pos, varname],
                    linewidth=1.1,
                    linestyle=l_styles[j],
                    color=colorVal,
                    label=sigmal_vals[j])
        else
            firm_val = pt.cvm_data[pos, varname][1]
            ax.axhline(y=firm_val,
                       linewidth=1.1,
                       linestyle=l_styles[j],
                       color=colorVal,
                       label=sigmal_vals[j],
                       xmin=.01, xmax=.99)
         end
    end
    
    if legend
        ax.legend(loc=0,
                   title="Volatility \$\\sigma\$",
                   ncol=1,
                   frameon=true,
                   shadow=true,
                   framealpha=.85,
                   edgecolor="white",
                   facecolor="white",
                   fancybox=true)
        # bbox_to_anchor=(.95, 0.85),

        ax.get_legend().get_title().set_color("#34495E")
    end
    # Axes' Labels:
    # if ylabel:
    #     ax.set_ylabel('Values for $\sigma = \overline{\sigma}$ and $\iota\geqslant 0$',
    #                    fontsize=12, labelpad=10)
    ax.set_xlabel("Risk Management Cost, \$\\iota\$ (b.p.)",
                   fontsize=12, labelpad=10)
    return ax
end


function cvm_double_axes_plot(pt, fig,
                              plot_sub_num::Int64, varname::Symbol;
                              kappa::Float64=NaN,
                              sigma_low::Float64=NaN,
                              y1label::Bool=true, y2label::Bool=true,
                              legend::Bool=true)
    (kappa_vals, sigmal_vals,
     iota_vals, scalarMap, l_styles) = cvm_data_handler(pt)

    if isnan(kappa)
        kappa = kappa_vals[1]
    end
    if isnan(sigma_low)
        sigma_low = minimum(sigmal_vals)
    end
    
    ax1 = fig.add_subplot(plot_sub_num)
    ax2 = ax1.twinx()
    for j in 1:size(sigmal_vals, 1)
        colorVal = scalarMap(j)
        pos = .&(abs.(pt.cvm_data[:kappa] .- kappa) .< 1e-6,
                 abs.(pt.cvm_data[:sigmal] .- sigmal_vals[j]) .< 1e-6)

        if (abs.(sigmal_vals[j] .- sigma_low) .< 1e-6)
            pos = .&(pos, [(x .* 10^4 in iota_vals) for x in pt.cvm_data[:iota]])


            ax1.plot(pt.cvm_data[pos, :iota] .* 10^4,
                     pt.cvm_data[pos, varname],
                     linewidth=1.1,
                     linestyle=l_styles[j],
                     color=colorVal,
                     label=sigmal_vals[j])
        else
            firm_val = pt.cvm_data[pos, varname][1]
            ax2.axhline(y=firm_val,
                        linewidth=1.1,
                        linestyle=l_styles[j],
                        color=colorVal,
                        label=sigmal_vals[j],
                        xmin=.01, xmax=.99)
        end
    end
    

    if legend
        # handles,labels = ax1.get_legend_handles_labels()
        handles, labels = [ ], [ ]
        # for ax in fig.axes:
        for ax in [ax1, ax2]
            for (h ,l) in zip(ax.get_legend_handles_labels()...)
                push!(handles, h)
                push!(labels, l)
            end
        end
        
        ax1.legend(handles, labels,
                   loc=0,
                   title="Volatility \$\\sigma\$",
                   ncol=1,
                   frameon=true,
                   shadow=true,
                   framealpha=.85,
                   edgecolor="white",
                   facecolor="white",
                   fancybox=true)
                   # bbox_to_anchor=(.95, 0.85),
    
        ax1.get_legend().get_title().set_color("#34495E")
    end
    
    # Axes' Labels:
    if y1label
        ax1.set_ylabel(string("Values for \$\\sigma = \\underline{\\sigma}\$",
                              " and \$\\iota \\geq 0\$"), usetex=true,
                       fontsize=12, labelpad=10)
    end
        
    if y2label
        ax2.set_ylabel(string("Values for \$\\sigma > \\underline{\\sigma}\$",
                              " and \$\\iota = 0\$"), usetex=true,
                       fontsize=12, labelpad=10)
    end
        
    ax1.set_xlabel("Risk Management Cost, \$\\iota\$ (b.p.)",
                   fontsize=12, labelpad=10)

    ax2.grid(nothing)

    return ax1

end


function cvm_plot_path_fname(pt, var::Symbol, fixed_params::Dict{Symbol, Any};
                             title_params_order::Array{Symbol,1}=cvm_plots_title_params_order,
                             main_dir::String=main_dir,
                             plots_dir::String=plots_dir,)

    # Paths to Plot ############################################
    # Main Directory
    main_dir_path = form_main_dir_path(main_dir)
    
    comb_dir_params = [x for x in title_params_order if x in keys(fixed_params)]
    comb_dir = join([str_format_fun(par_val_printer(x), 
                                    par_val_adj(x, fixed_params[x])) 
                     for x in comb_dir_params], "__")

    # Form Paths ###############################################
    dirs = [plots_dir, "CVM", comb_dir]
    graphs_path = main_dir_path
    for dir in dirs
        graphs_path = string(graphs_path, "/", dir)
        if !isdir(graphs_path)
            mkdir(graphs_path)
        end
    end
    # ##########################################################

    # Filename #################################################
        
    file_name = string("cvm_", obj_fun_dict[fixed_params[:obj_fun]],
                       "_opt_", var)
    
    return [graphs_path, file_name]
end


function plot_cvm_optimal_solutions(pt, var::Symbol;
                                    kappa::Float64=NaN,
                                    sigma_low::Float64=NaN,
                                    title_params_order::Array{Symbol,1}=cvm_plots_title_params_order,
                                    figaspect::Float64=.55,
                                    facecolor::String="w",
                                    save_fig::Bool=true,
                                    return_fig::Bool=true,
                                    fig_dpi::Int64=400,
                                    graph_format::String="png")

    if isnan(kappa)
        kappa = minimum(unique(pt.cvm_data[:kappa]))
    end
    if isnan(sigma_low)
        sigma_low = minimum(unique(pt.cvm_data[:sigmal]))
    end

    Seaborn.set_style("darkgrid")
    fig = PyPlot.figure(figsize=PyPlot.figaspect(figaspect), facecolor=facecolor)
    if !(var in [:firm_value, :MBR])
        ax1 = ModelPlots.cvm_double_axes_plot(pt, fig, 111, var;
                                              kappa=kappa, sigma_low=sigma_low)
    else
        ax1 = ModelPlots.cvm_single_axis_plot(pt, fig, 111, var;
                                              kappa=kappa, sigma_low=sigma_low)
    end
    
    ax1.set_title(string("Constant Volatility Model Optimal ", vartitles[var]), fontsize=14)
    PyPlot.tight_layout()


    if save_fig
        fixed_params = Dict{Symbol, Any}(:obj_fun => Symbol(pt.cvm_data[1, :obj_fun]),
                                         :mu_b => unique(pt.cvm_data[:mu_b])[1],
                                         :m => unique(pt.cvm_data[:m])[1],
                                         :xi => unique(pt.cvm_data[:xi])[1],
                                         :kappa => kappa,
                                         :sigmal => sigma_low)
        
        folder_path, file_name = cvm_plot_path_fname(pt, var, fixed_params;
                                                     title_params_order=title_params_order)
        PyPlot.savefig(string(folder_path, "/", file_name, ".", graph_format), dpi=fig_dpi, format=graph_format)
    end
    
    if return_fig
        return fig
    end
end
