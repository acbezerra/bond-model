

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


function interp_z_values(df::DataFrame, xvar::Symbol, 
                         yvar::Symbol, zvars::Array{Symbol, 1};
                         ft_xy::Symbol=Symbol(""),
                         spline_k::Int64=3, 
                         spline_bc::String="extrapolate")
    # Separate DataFrames
    xdf = slice_df(df, Symbol(ft_xy, xvar))
    ydf = slice_df(df, Symbol(ft_xy, yvar))

    # Form Dictionary to store results
    tmpdict = Dict(zip([xvar, yvar], 
                       [Spline1D(1:5, 1:5), Spline1D(1:5, 1:5)]))
    # fd = Dict(zip(zvars, repeat([deepcopy(tmpdict)], 1, size(zvars, 1))))

    fd = Dict()
    fd[:xvals] = xdf[Symbol(ft_xy, xvar)]
    fd[:yvals] = ydf[Symbol(ft_xy, yvar)]
    
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



