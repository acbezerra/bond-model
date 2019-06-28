

function get_eq_type_df(eq_type::String, 
                        fidf::DataFrame, misrepdf::DataFrame, 
                        sepdf::DataFrame, pooldf::DataFrame)
    if eq_type == "full_info"
        return fidf
    elseif eq_type == "misrep"
        return misrepdf
    elseif eq_type == "separating"
        return sepdf
    elseif eq_type == "pooling"
        return pooldf
    end    
end


function slice_df(df::DataFrame, svar::Symbol; tol::Float64=1e-5)
    if any(isnan.(df[svar]))
        return df[isnan.(df[svar]) .==false, :]
    elseif any(abs.(df[svar] .- .0) .< tol)
        return df[abs.(df[svar] .- .0) .> tol, :]
    end
end


function interp_z_values(df::DataFrame;
                         xvar::Symbol=contour_xvar,
                         yvar::Symbol=contour_yvar,
                         zvars::Array{Symbol, 1}=contour_zvars,
                         ft_xy::Symbol=Symbol(""),
                         ft_z::Array{Symbol, 1}=[:s_, :r_],
                         spline_k::Int64=3, 
                         spline_bc::String="extrapolate")
    if !(xvar in names(df))
        ft_xy = :r_
    end
    
    # Separate DataFrames
    xdf = slice_df(df, Symbol(ft_xy, xvar))
    ydf = slice_df(df, Symbol(ft_xy, yvar))

    # Form Dictionary to store results
    tmpdict = Dict{Symbol, Any}(zip([xvar, yvar], 
                                    [Spline1D(1:5, 1:5), Spline1D(1:5, 1:5)]))
    # fd = Dict(zip(zvars, repeat([deepcopy(tmpdict)], 1, size(zvars, 1))))

    fd = Dict{Symbol, Any}(:xvar => xvar,
                           :yvar => yvar,
                           :xvals => xdf[Symbol(ft_xy, xvar)],
                           :yvals => ydf[Symbol(ft_xy, yvar)])

    if !(zvars[1] in names(df))
        zvars = vcat([Symbol(prefix_z, zvar) for prefix_z in ft_z, zvar in zvars]...)
    end
    
    # Interpolate Functions
    for zvar in zvars
        fd[zvar] = deepcopy(tmpdict)
        fd[zvar][xvar] = Dierckx.Spline1D(fd[:xvals], xdf[zvar], 
                                          k=spline_k, bc=spline_bc)
        fd[zvar][yvar] = Dierckx.Spline1D(fd[:yvals], ydf[zvar], 
                                          k=spline_k, bc=spline_bc)
    end
        
    return fd
end


function form_mesh_grid(xvals::Array{Float64,1},
                        yvals::Array{Float64,1},
                        zfun; N::Int64=200)
    xgrid = range(minimum(xvals), stop=maximum(xvals), length=N)
    ygrid = range(minimum(yvals), stop=maximum(yvals), length=N)

    X = Array(repeat(xgrid, 1, N)')
    Y = repeat(ygrid, 1, N)
    Z = Array([zfun(x,y) for x in xgrid, y in ygrid]')

    return X, Y, Z
end


function get_contour_plot_path_name(df::DataFrame, zfun_name::Symbol;
                                    firm_type::Symbol=Symbol(""), fname_eq_type::String="",
                                    fname_ext::String=contour_fname_ext)
    eq_type = ""
    if isempty(fname_eq_type)
        eq_type = df[1, :eq_type]
        fname_eq_type = eq_type_title[eq_type][1]
    end
    
    if .&(eq_type != "full_info", any([isempty(string(firm_type)),
                                      !(firm_type in [:safe, :risky])]))
        println("Please enter a firm type: :safe or :risky. Exiting...")
        return
    end


    mu_s = NaN
    if eq_type == "full_info"
        ft = ""
        iota_s = minimum([x for x in df[:iota] if x > 0.])    
        lambda_r = unique([x for x in df[:lambda] if !isnan(x)])[1]
    else
        ft = (firm_type == :safe) ? :s_ : :r_
        mu_s = df[1, :mu_s]
        iota_s = df[1, :s_iota]
        lambda_r = unique([x for x in df[:r_lambda] if !isnan(x)])[1]
    end

    # Inverse Coupon Rate
    pcr = df[1, :p]/df[1, :c]
    
    fname = string(fname_eq_type, "_mu_s_", mu_s, "_pcr_", pcr, 
                   "_iota_s_", iota_s, "__lambda_r_", lambda_r,
                   "__", ft,  zfun_name)
    return string(contour_plots_path, "/", fname, ".", fname_ext)
end


function get_region_grids(xgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                          ygrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                          eqfun; N::Int64=10^5)
    ymin = minimum(ygrid)
    ymax = maximum(ygrid)
    
    yvals = x -> [y for y in ygrid if eqfun(x, y)]
    min_y = x -> !isempty(yvals(x)) ? minimum(yvals(x)) : .0
    max_y = x -> !isempty(yvals(x)) ? maximum(yvals(x)) : .0
    
    min_y_grid = fetch(@spawn [minimum([maximum([ymin, min_y(x)]), ymax]) for x in xgrid])
    max_y_grid = fetch(@spawn [minimum([maximum([ymin, max_y(x)]), ymax]) for x in xgrid])
    min_y_fun = Dierckx.Spline1D(xgrid, min_y_grid; k=3, bc="extrapolate")
    max_y_fun = Dierckx.Spline1D(xgrid, max_y_grid; k=3, bc="extrapolate")

    ref_xgrid = range(minimum(xgrid), stop=maximum(xgrid), length=N)
    min_y_ref_grid = [minimum([maximum([ymin, min_y_fun(x)]), ymax]) for x in ref_xgrid]
    max_y_ref_grid = [minimum([maximum([ymin, max_y_fun(x)]), ymax]) for x in ref_xgrid]
    
    # return min_y, max_y, ref_xgrid, min_y_ref_grid, max_y_ref_grid
    return ref_xgrid, min_y_ref_grid, max_y_ref_grid
end



function get_eq_contour_mesh_grid(xvals::Array{Float64,1}, yvals::Array{Float64,1},
                                  fun_dict; N::Int64=10^3)

    # eq_bool = (x, y) -> fun_dict[:fi_ind](x, y) + 2 * fun_dict[:sep_ind](x, y) + 3 * fun_dict[:pool_ind](x, y)
    # eq_vals = (x, y) -> (fun_dict[:fi_ind](x,y) * fun_dict[:mbr][:fi](x, y) +
    #                      fun_dict[:sep_ind](x,y) * fun_dict[:mbr][:sep](x, y) +
    #                      fun_dict[:pool_ind](x,y) * fun_dict[:mbr][:pool](x, y))
    
    X, Y, bool_Z = form_mesh_grid(xvals, yvals, fun_dict[:eq_bool], N=N)
    _, _, bool_OTC_EP = form_mesh_grid(xvals, yvals, fun_dict[:bool_otc_ep], N=N)
    _, _, r_MBR = fetch(@spawn form_mesh_grid(xvals, yvals, fun_dict[:r_mbr], N=N))
    _, _, s_FV = fetch(@spawn form_mesh_grid(xvals, yvals, fun_dict[:s_fv], N=N))
    
    return Dict{Symbol, Any}(:X => X, :Y => Y,
                             :bool_Z => bool_Z,
                             :bool_OTC_EP => bool_OTC_EP,
                             :r_MBR => r_MBR,
                             :s_FV => s_FV)
end


# function find_xvar_yvar_values(df::DataFrame, xvar::Symbol, 
#                                yvar::Symbol, zvar::Symbol;
#                                N1::Int64=100, N2::Int64=10^5,
#                                spline_k::Int64=3, 
#                                spline_bc::String="extrapolate")
#     # Separate DataFrames
#     xdf = slice_df(df, xvar)
#     ydf = slice_df(df, yvar)
    
#     # Form Grids
#     xgrid = range(minimum(xdf[xvar]), stop=maximum(xdf[xvar]), length=N2)
#     ygrid = range(minimum(ydf[yvar]), stop=maximum(ydf[yvar]), length=N2)

#     # Interpolate Functions
#     xzfun = Dierckx.Spline1D(xdf[xvar], xdf[zvar], 
#                              k=spline_k, bc=spline_bc)
#     yzfun = Dierckx.Spline1D(ydf[yvar], ydf[zvar], 
#                              k=spline_k, bc=spline_bc)

#     # Form grid of yvals
#     xvals = []
#     yvals = [] 
#     for zval in zgrid
#         push!(xvals, xgrid[argmin(abs.(xzfun(xgrid) .- zval))])
#         push!(yvals, ygrid[argmin(abs.(yzfun(ygrid) .- zval))])
#     end

#     # Interpolate Results
#     tmp = DataFrame(Dict(zip([xvar, yvar, zvar], 
#                              [xvals, yvals, zgrid])))
#     sort!(tmp, [xvar])
#     xyfun = Dierckx.Spline1D(tmp[xvar], tmp[yvar], k=spline_k, bc=spline_bc)
        
#     return xzfun, yzfun, xyfun, tmp
# end



# function store_contour_results(xvar::Symbol, yvar::Symbol,
#                                fidf::DataFrame, misrepdf::DataFrame, 
#                                sepdf::DataFrame, pooldf::DataFrame;
#                                safe_zvar::Symbol=:firm_value,
#                                risky_zvar::Symbol=:MBR,
#                                N1::Int64=100, N2::Int64=10^5)

#     fd = Dict{Symbol, Dict}(:safe => Dict("full_info" => deepcopy(resdict), 
#                                           "pooling" => deepcopy(resdict),
#                                           "separating" => deepcopy(resdict)),
#                             :risky => Dict("full_info" => deepcopy(resdict),
#                                            "misrep" => deepcopy(resdict)))

#     for firm_type in keys(fd)
#         zvar = (firm_type == :safe) ? safe_zvar : risky_zvar
#         for eq_type in keys(fd[firm_type])
#             df = get_eq_type_df(eq_type, fidf, misrepdf, sepdf, pooldf)

#             v_z_prefix = ""
#             if eq_type != "full_info"
#                v_z_prefix = (firm_type == :safe) ? :s_ : :r_
#             end
#             v_xy_prefix = (eq_type != "full_info") ? :r_ : ""
            
#             x_fun, y_fun, xy_fun, tmpdf = find_xvar_yvar_values(df, Symbol(v_xy_prefix, xvar), 
#                                                                 Symbol(v_xy_prefix, yvar), 
#                                                                 Symbol(v_z_prefix, zvar);
#                                                                 N1=N1, N2=N2)

#             fd[firm_type][eq_type][:xvar] = Symbol(v_xy_prefix, xvar)
#             fd[firm_type][eq_type][:yvar] = Symbol(v_xy_prefix, yvar)
#             fd[firm_type][eq_type][:zvar] = Symbol(v_z_prefix, zvar)
#             fd[firm_type][eq_type][:x_fun] = x_fun
#             fd[firm_type][eq_type][:y_fun] = y_fun
#             fd[firm_type][eq_type][:xy_fun] = xy_fun
#             fd[firm_type][eq_type][:df] = tmpdf
#         end
#     end

#     return fd
# end



