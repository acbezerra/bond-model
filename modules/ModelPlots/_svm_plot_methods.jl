
function set_svm_surf(pt, fixed_params::Dict{Symbol, Float64})
    if isempty(pt._svm_surf)
        println("Setting SVM Slice Data")
    else
        println("Updating SVM Slice Data")
    end

    pt._svm_surf = plot_svmdf_slicer(pt, fixed_params)

    return pt
end


function svm_interp_smooth_surface(pt, fixed_params::Dict{Symbol,Float64},
                                   z_var::Symbol;
                                   interp::Bool=true,
                                   smooth::Bool=true)

    # Variables
    xy_vars = [:kappa, :lambda, :sigmah]

    kappa_vec = unique(pt._svm_surf[:kappa])
    sigmah_vec = unique(pt._svm_surf[:sigmah])
    lambda_vec = unique(pt._svm_surf[:lambda])

    # xy axes - pick 2 variables from (kappa, lambda, sigmah):
    xy_list = [var for var in xy_vars if 
                !(var in keys(fixed_params))]

    # Sigmah on Y-axis, Kappa on X-axis:
    # if Sigmah & Lambda -> lambda on X-Axis
    if .&((:kappa in xy_list), (:sigmah in xy_list))
        x_vec = kappa_vec
        y_vec = sigmah_vec
        var_axes = [:kappa, :sigmah]
    elseif .&((:kappa in xy_list), (:lambda in xy_list))
        x_vec = kappa_vec
        y_vec = lambda_vec
        var_axes = [:kappa, :lambda]
    else
        x_vec = lambda_vec
        y_vec = sigmah_vec
        var_axes = [:lambda, :sigmah]
    end

    # Form Combinations:
    comb = hcat([Array([x, y]) for x in x_vec, y in  y_vec]...)'

    # #####################################################################
    # ####################### Extract z-var values: #######################
    # #####################################################################
    z_val = Array{Float64}(undef, size(x_vec, 1) .* size(y_vec, 1), 1)

    # Boolean Indexer for Fixed Parameter Values:            
    pos_loc = [(abs.(pt._svm_surf[x] .- fixed_params[x]) .< 1e-6) 
                                for x in keys(fixed_params)]
    bool_loc =  .&(pos_loc...) # sum(pos_loc) .== size(pos_loc, 1)

    for i in 1:size(comb, 1)
        if .&((:kappa in xy_list), (:sigmah in xy_list))
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:kappa] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:sigmah] .- comb[i, 2]) .< 1e-6),
                                    bool_loc), z_var][1]
        elseif .&((:kappa in xy_list), (:lambda in xy_list))
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:kappa] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:lambda] .- comb[i, 2]) .< 1e-6),
                                     bool_loc), z_var][1]
        else
            value = pt._svm_surf[.&((abs.(pt._svm_surf[:lambda] .- comb[i, 1]) .< 1e-6),
                                    (abs.(pt._svm_surf[:sigmah] .- comb[i, 2]) .< 1e-6),
                                    bool_loc), z_var][1]
        end

        # Keep in mind that not all combinations of x_vec and
        # y_vec values may be available.
        if isempty(value)
            z_val[i] = NaN
        else
            z_val[i] = value
        end
    end

    # #####################################################################
    # ########################### Form Surface ############################
    # #####################################################################
    # Form Surface
    xyz = hcat(comb, z_val)

    # Remove combinations for which output is NaN:
    xyz[vec(sum(isnan.(xyz), dims = 2) .== 0), :]

    # Interpolate Surface
    if interp
        # Spline: Setting s to length xyz gives a SmoothBivariateSpline
        xyz_interp = Dierckx.Spline2D(xyz[:, 1], xyz[:, 2], xyz[:, 3]; 
                                      kx=3, ky=3,s=size(xyz, 1))

        # Increase Data Points

        # To spot inconsistencies and errors in the
        # manipulation of grid points, I use
        # rectangular xy grid:
        x_step_num = 100
        y_step_num = 105
        x_grid_ref = range(minimum(xyz[:, 1]), stop=maximum(xyz[:, 1]), length=x_step_num)
        y_grid_ref = range(minimum(xyz[:, 2]), stop=maximum(xyz[:, 2]), length=y_step_num)
        xy_grid_ref = hcat([Array([x, y]) for x in x_grid_ref, y in y_grid_ref]...)'
        grid_x = xy_grid_ref[:, 1]
        grid_y = xy_grid_ref[:, 2]

        # Back-out Interpolated Values:
        grid_z = evaluate(xyz_interp, grid_x, grid_y)
    else
        grid_x = xyz[:, 1]
        grid_y = xyz[:, 2]
        grid_z = xyz[:, 3]
    end

    # Kernel Surface Smoother
#     f = None
#     if smooth:
#         f = Rbf(grid_x, grid_y, grid_z, smooth=0.15, epsilon=1.)

    return [xyz, grid_x, grid_y, grid_z, var_axes]
end


function svm_plot_heatmap(pt, xy_list, z_var::Symbol,
                          fixed_params::Dict{Symbol,Float64};
                          plt_cmap::String="viridis",
                          order_ascend::Bool=true,
                          reverse_x_axis::Bool=false,
                          reverse_y_axis::Bool=false,
                          make_title::Bool=true,
                          add_title::Bool=false)
    # ############################################
    # ############### Axis Format: ###############
    # ############################################
    title_list = [string("\$", pt.xylabels[x][1], "=\$ ", pt.xylabels[x][2]) 
                        for x in keys(pt.xylabels) if !(x in xy_list)]
    title_values = [fixed_params[x] for x in keys(pt.xylabels) if !(x in xy_list)]                                        
     

    svmdf = deepcopy(pt._svm_surf)
                
    # Localize entries:
    # Localize before changing kappa to basis points:
    cond = true
    for key in keys(fixed_params)
        cond = .&(cond, (abs.(svmdf[key] .- fixed_params[key]) .< 1e-6))
    end

    # Sigmah on Y-axis, Kappa on X-axis:
    # if Sigmah & Lambda -> lambda on X-Axis
    # Adjust value of Kappa
    if .&((:sigmah in xy_list), (:kappa in xy_list))
        cols = vcat([:sigmah, :kappa], [z_var])
        # Kappa in Basis Points:
        svmdf[:kappa] .= convert.(Int64, svmdf[:kappa] .* 1e4)
    elseif .&((:sigmah in xy_list), (:lambda in xy_list))
        cols = vcat([:sigmah, :lambda], [z_var])
        # Kappa in Basis Points
        # fixed_params[:kappa] = #convert.(Int64, fixed_params[:kappa] .* 1e4)
    else
        cols = vcat([:lambda, :kappa], [z_var])
        # Kappa in Basis Points:
        svmdf[:kappa] .=  convert.(Int64, svmdf[:kappa] .* 1e4)
    end

    # Pivot Table for Plot: ##############################
    xtickvals = unique(svmdf[cols[2]])
    ytickvals = unique(svmdf[cols[1]])
    pivotdf = unstack(svmdf[cond, cols], cols[1], cols[2], cols[3])
    Z = convert(Matrix, pivotdf[:, [x for x in names(pivotdf) if !(x in cols)]])

    if reverse_x_axis
        Z = reverse(Z, dims=2)
        xtickvals = reverse(xtickvals)
    end

    if reverse_y_axis
        Z = reverse(Z, dims=1)
        ytickvals = reverse(ytickvals)
    end

    ax = Seaborn.heatmap(Z, cmap=plt_cmap,
                          xticklabels=xtickvals,
                          yticklabels=ytickvals)
    ax.set_xlabel(string("\$", pt.xylabels[cols[2]][1], "\$"))
    ax.set_ylabel(string("\$", pt.xylabels[cols[1]][1], "\$"))
    # ##################################################
                                             
                                         
    # ############### TITLE ###############
    ax_title = " "
    if make_title | add_title
        # str_format_fun(a, b) = @eval @sprintf($a, $b)
        formatted_var_string = join([str_format_fun(title_list[i], 
                                                    title_values[i]) 
                                     for i in 1:size(title_list, 1)], ", ")
        ax_title = string(pt.zlabels[z_var][1], " values for ",
                          formatted_var_string)                                
    end

    if add_title
        ax.set_title(ax_title)
    end

    return [ax, ax_title]

end


function svm_plot_surface(pt, ax, xv, yv, zv,
                          xy_list::Array{Symbol,1},
                          z_var::Symbol,
                          fixed_params::Dict{Symbol,Float64};
                          title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                          plt_cmap::String="viridis",
                          seaborn_style::String="darkgrid",
                          make_title=true, add_title=false,
                          zpad=10, view_elev=25., view_azim=210)
                          #xylabels::Array{Symbol,1}=[:mu_b, :m, :xi, :kappa, :lambda, :sigmal, :sigmah])
    # Documentation available at: https://matplotlib.org/mpl_toolkits/mplot3d/api.html

    # ATTENTION! -> only adjust kappa after computing zv:
    # Set Kappa to Basis Points:
    if xy_list[2] == :kappa
        xv = convert.(Int64, xv .* 10^4)
    elseif xy_list[1] == :kappa
        yv = convert.(Int64, yv .* 10^4)
    end

    cp_fixed_params = deepcopy(fixed_params)
    if :kappa in keys(cp_fixed_params)
        cp_fixed_params[:kappa] = convert.(Int64, cp_fixed_params[:kappa] .* 10^4)
    end
    # ############################################
    # ############### Axis Format: ###############
    # ############################################
    major_xticks = range(minimum(xv), stop=maximum(xv), length=6)
    major_yticks = range(minimum(yv), stop=maximum(yv), length=6)
    major_zticks = range(ceil(minimum(zv)), stop=ceil(maximum(zv)), length=10)
    if z_var == :c
        major_zticks = range(minimum(zv), stop=maximum(zv), length=7)
    end

    xy_format = [[string("\$",  pt.xylabels[x][1], "\$"), pt.xylabels[x][2]] 
                    for x in keys(pt.xylabels) if (x in xy_list)]

    title_list = [string("\$", pt.xylabels[x][1], "=\$ ", pt.xylabels[x][2]) 
                        for x in title_params_order if !(x in xy_list)]
    title_values = [cp_fixed_params[x] for x in title_params_order if !(x in xy_list)]

    if isempty(seaborn_style)
        Seaborn.reset_orig()
    else
        Seaborn.set(style="darkgrid")
    end
    surf = ax.plot_trisurf(xv, yv, zv,
                           linewidth=0.2,
                           cmap=plt_cmap,
                           antialiased=true)                    

    # Customize the x axis
    # ax.xaxis_inverted()
    ax.invert_xaxis()
    ax.xaxis.set_rotate_label(false)  # disable automatic rotation
    ax.set_xlabel(xy_format[1][1], labelpad=10, rotation=0)
    ax.xaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(xy_format[1][2]))
    ax.set_xticks(major_xticks)

    # Customize the y axis
    ax.yaxis.set_rotate_label(false)  # disable automatic rotation
    ax.set_ylabel(xy_format[2][1], labelpad=10, rotation=0)
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter(xy_format[2][2]))
    ax.set_yticks(major_yticks)

    # Customize the z axis
    ax.zaxis.set_rotate_label(false)  # disable automatic rotation
    ax.zaxis.set_major_formatter(PyPlot.matplotlib.ticker.FormatStrFormatter("%.2f"))
    ax.set_zlabel(pt.zlabels[z_var][1], rotation=90, labelpad=zpad)
    ax.set_zticks(major_zticks)
    ax.tick_params(axis="z", pad=zpad/2)

    ax.minorticks_on()
    ax.grid(which="minor", alpha=0.3)
    # ############### TITLE ###############
    ax_title = " "
    if make_title | add_title
        # str_format_fun(a, b) = @eval @sprintf($a, $b)
        formatted_var_string = join([str_format_fun(title_list[i], 
                                                    title_values[i]) 
                                     for i in 1:size(title_list, 1)], ", ")
        ax_title = latexstring(pt.zlabels[z_var][1], " values for ",
                          formatted_var_string)                                
    end

    if add_title
        ax.set_title(ax_title)
    end

#     # Add a color bar which maps values to colors.
#     # fig.colorbar(surf, shrink=0.6, aspect=6)

    ax.view_init(view_elev, view_azim)

    return surf, ax, ax_title
end


function svm_heat_surf_plot_path_fname(pt, xy_list::Array{Symbol,1}, z_var::Symbol, 
                                       fixed_params::Dict{Symbol,Float64};
                                       title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                                       main_dir::String=main_dir,
                                       plots_dir::String=plots_dir,
                                       graph_dir::String=heat_surf_graph_dir,
                                       mat_dir_prefix::String=mat_dir_prefix)
    
    # Functions to Format Numeric-to-String
    # str_format_fun(a::String,b::Float64) = @eval @sprintf($a, $b)
    # par_val_printer(x::Symbol) = !(x in [:iota, :kappa]) ? string(x, "_", pt.xylabels[x][2]) : string(x, "_(bp)_", pt.xylabels[x][2])
    # par_val_adj(x::Symbol, val::Float64) = !(x in [:iota, :kappa]) ? val : val * 1e4
    
    # Paths to Plot ##########################################
    
    # Main Directory
    main_dir_path = form_main_dir_path(main_dir)
    
    # Maturity Directory
    mat_dir = string(mat_dir_prefix, str_format_fun(pt.xylabels[:m][2], fixed_params[:m]))

    # Get Graph Type 
    graph_type = ""
    graph_sub_folder = ""
    if .&((:kappa in xy_list), (:sigmah in xy_list))
        graph_type = "kappa_sigmah"
        graph_sub_folder = str_format_fun(string("lambda_", pt.xylabels[:lambda][2]),
                                          fixed_params[:lambda])
    elseif .&((:kappa in xy_list), (:lambda in xy_list))
        graph_type = "kappa_lambda"
        graph_sub_folder = str_format_fun(string("sigmah_", pt.xylabels[:sigmah][2]),
                                          fixed_params[:sigmah])
    else
        graph_type = "lambda_sigmah"
        graph_sub_folder = str_format_fun(string("kappa_(bp)_", pt.xylabels[:kappa][2]),
                                          fixed_params[:kappa] * 10^4)
    end
    
    # Form Paths ##############################################
    dirs = [plots_dir, "SVM", graph_dir, mat_dir,
            graph_type, graph_sub_folder]
    graphs_path = main_dir_path
    for dir in dirs
        graphs_path = string(graphs_path, "/", dir)
        if !isdir(graphs_path)
            mkdir(graphs_path)
        end
    end
    # ##########################################################

    # Filename #################################################
    if isempty(title_params_order)
        title_params_order = keys(fixed_params)
    end

    title_params = [x for x in title_params_order if x in keys(fixed_params)]
    par_values_str = join([str_format_fun(par_val_printer(k), 
                           par_val_adj(k, fixed_params[k])) 
                      for k in title_params], "__")
    file_name = string("svm_", obj_fun_dict[Symbol(pt.svm_data[1, :obj_fun])],
                       "_", z_var, "_", par_values_str)

    return [graphs_path, file_name]
end


function svm_plot_heatmap_surf(pt, xy_list::Array{Symbol, 1}, 
                               z_var::Symbol, 
                               fixed_params::Dict{Symbol,Float64};
                               title_params_order::Array{Symbol,1}=svm_plots_title_params_order,
                               save_fig::Bool=true,
                               elev::Float64=25.,
                               azim::Float64=210.,
                               zpad::Float64=10.,
                               plt_cmap::String="viridis", # 'Spectral', 'cool',
                               seaborn_style::String="darkgrid",
                               heat_reverse_x::Bool=false,
                               heat_reverse_y::Bool=true,
                               interp_bool::Bool=false,
                               smooth_bool::Bool=false,
                               ax1_dist::Float64=8.5,
                               axs_wspace::Float64=.1,
                               cbaxes::Array{Float64,1}=[.975, 0.15, 0.015, 0.675],
                               return_fig::Bool=true,
                               fig_dpi::Int64=400,
                               graph_format::String="eps")
# sup_title_x::Float64=.575,
    
    # ###################################################
    # ################### MULTI PLOT ####################
    # ###################################################
    using3D()
    fig = PyPlot.figure("pyplot_surfaceplot", figsize=PyPlot.figaspect(.4), facecolor="w")

    # ###################################################
    # ##################### SURFACE #####################
    # ###################################################
    ax1 = fig.add_subplot(1, 2, 2, projection="3d")

    xyz, grid_x, grid_y, grid_z, var_axes = svm_interp_smooth_surface(pt, 
                                                                      fixed_params, 
                                                                      z_var,
                                                                      interp=interp_bool,
                                                                      smooth=smooth_bool)

    # Notice that variables X and Y are listed in var_axes and follow the
    # convention of sigmah on the Y-Axis and Kappa on the X-Axis. When passing
    # the list of variables to the plot_surface function, use var_axes list:
    surf, ax1, fig_title = svm_plot_surface(pt, ax1, xyz[:,1], xyz[:,2], xyz[:,3], 
                                            xy_list, z_var, fixed_params;
                                            title_params_order=title_params_order,
                                            plt_cmap=plt_cmap, 
                                            seaborn_style=seaborn_style,
                                            make_title=true, 
                                            add_title=false,
                                            zpad=zpad, view_elev=elev, view_azim=azim)

    # Set the background color of the pane YZ
    ax1.w_xaxis.set_pane_color(PyPlot.matplotlib.colors.hex2color("#d5d8dc"))

    # Add a color bar which maps values to colors.
    #    cb = fig.colorbar(surf, aspect=20, ax=ax1)
    cbax = fig.add_axes(cbaxes) 
    cb = fig.colorbar(surf, aspect=20, cax=cbax)
    cb.outline.set_visible(false)

    ax1.patch.set_facecolor("white")
    ax1.dist = ax1_dist
    # ###################################################

    # ###################################################
    # ##################### HEATMAP #####################
    # ###################################################
    ax2 = fig.add_subplot(1, 2, 1)


    ax2, _ = svm_plot_heatmap(pt, xy_list, z_var, fixed_params;
                              plt_cmap=plt_cmap,
                              reverse_x_axis=heat_reverse_x,
                              reverse_y_axis=heat_reverse_y)
    # ###################################################

    fig.suptitle(fig_title, fontsize=14) #, x=sup_title_x)
    PyPlot.subplots_adjust(wspace=axs_wspace)
    #PyPlot.show()

    if save_fig
        m = unique(pt._svm_surf[:m])[1]
        folder_path, file_name = svm_heat_surf_plot_path_fname(pt, xy_list, 
                                                               z_var, fixed_params;
                                                               title_params_order=title_params_order)
        PyPlot.savefig(string(folder_path, "/", file_name, ".", graph_format), dpi=fig_dpi, format=graph_format)
    end

    if return_fig
        return fig 
    end
end
    
