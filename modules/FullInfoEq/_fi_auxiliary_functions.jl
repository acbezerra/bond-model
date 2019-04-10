
function full_info_eq_deriv_root_search(svm, df; N::Int64=10^5,
                                        k::Int64=3, bc::String="extrapolate")
    
    # Interpolate Values
    eq_deriv_fun = Dierckx.Spline1D(df[:vb], df[:eq_deriv], k=3, bc="extrapolate")
    
    # Compute optimal VB:
    res = DataFrame()
    vbroots = roots(eq_deriv_fun; maxn=8)
    if !isempty(vbroots)
        eq_min_val_interp = Dierckx.Spline1D(df[:vb], df[:eq_min_val], k=3, bc="extrapolate")
        abs_eq_min_val = abs.(eq_min_val_interp(vbroots))
        res[:vb] = vbroots[argmin(abs_eq_min_val)]
    else
        ref_vbgrid = range(minimum(df[:vb]), stop=maximum(df[:vb]), length=N)
        res[:vb] = ref_vbgrid[argmin(abs.(eq_deriv_fun(ref_vbgrid)))]
    end
    
    # Equity Values
    res[:eq_deriv] = eq_deriv_fun(res[:vb])

    # Interpolate Functions
    interp_cols = vcat([:debt, :equity],
                       [:eq_deriv_min_val, 
                        :eq_min_val, :eq_vb])
    res = interp_values(res, df, :vb, interp_cols; k=k, bc=bc)

    res[:eq_negative] = (res[:eq_min_val] .< -.005)

    # Fixed Values
    return  res #non_interp_values(svm, res)
end

