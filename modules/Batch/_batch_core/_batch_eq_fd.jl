

function df_slicer(df, p)
   return df[abs.((df[:p] .- p)) .< 1e-4, :] 
end


function abs_debt_diff_cond(df, p; tol=.05)
    # Slice DataFrame
    return minimum(df_slicer(df, p)[:abs_debt_diff]) < tol
end


function filter_debt_values(bt, svm, df; N1=50, N2=10)
    pval = unique(df[:p])[1]
    
    debtf = Dierckx.Spline1D(df[:vb], df[:debt], k=3, bc="extrapolate")
    vbgrid = range(minimum(df[:vb]), stop=maximum(df[:vb]), length=N1)
    aggP = get_agg_p(svm, p=pval)
    debt_diff = debtf(vbgrid) .- aggP
    
    # Find VB values for which the absolute debt diff is the smaller:
    pos = partialsortperm(abs.(debt_diff), 1:N2)
    vbvals = vbgrid[pos]
    
    sdf = DataFrame(vb=vbvals, debt=debtf(vbvals))
    sdf[:c] = unique(df[:c])[1]
    sdf[:p] = unique(df[:p])[1]                        
    sdf = debt_at_par_diffs(svm, sdf, :vb, p=pval)
    
    return sdf[vcat([:c, :p, :vb, :debt], bt.dfc.debt_vars)]
end


function filter_batch_I_df(bt, svm, df; tol=.05, N1=50, N2=10)
    LL = fetch(@spawn [filter_debt_values(bt, svm, df_slicer(df, p); N1=N1, N2=N2) 
                       for p in unique(df[:p])])
    sdf=vcat(LL...)
    
    # Get Filtered pgrid
    pgrid =  [p for p in unique(sdf[:p]) if 
              (minimum(abs.(df_slicer(sdf, p)[:debt_diff])) < tol)] 
    while size(pgrid, 1) < 5
        tol = 1.25 * tol
        pgrid =  [p for p in unique(sdf[:p]) if 
                  (minimum(abs.(df_slicer(sdf, p)[:debt_diff])) < tol)] 
    end
                                                        
    # Return Filtered DataFrame
    return df[findall(in(pgrid), df[:p]), :]
end


function eq_fd_method(bt, svm, df; mu_b=NaN, c=NaN)
    if isnan(mu_b)
        mu_b=svm.mu_b 
    end
    
    if isnan(c)
        c=svm.c
    end
    
    res = @time fetch(@spawn [eq_fd(svm, df[i, :vb]; mu_b=mu_b, c=c, p=df[i, :p]) 
                              for i in 1:size(df, 1)])

    # Collect Results
    tmp = vcat(res...)

    # Compute Debt Differences
    eqdf_all = debt_at_par_diffs(svm, tmp, :p; sort_vars=[:p, :vb])

    # Rearrange Columns
    return eqdf_all[bt.dfc.dfcols]
end


function interp_values(res, df, xvar::Symbol, interp_cols::Array{Symbol,1};
                       k::Int64=3, bc::String="extrapolate")
    for col in interp_cols
        ffun = Dierckx.Spline1D(df[xvar], df[col]; k=k, bc=bc)
        res[col] = ffun(res[xvar])
    end 
    return res
end


function non_interp_values(svm, df)
    # Debt Values
    df[:abs_debt_diff] = abs.(df[:debt_diff])
    df[:abs_debt_per_diff] = abs.(df[:debt_per_diff])

    # Equity
    df[:eq_negative] = (df[:eq_min_val] .< -.005)

    # Share Values
    df[:firm_value] = df[:debt] .+ df[:equity]
    df[:leverage] = (df[:debt]./df[:firm_value]) .* 100
    df[:ROE] = (df[:equity]./(svm.pm.V0 .- df[:debt]) .- 1) .* 100

    return df
end


function eq_deriv_root_search(svm, df, p; mu_b=NaN, c=NaN, N=10^5)
    if isnan(mu_b)
        mu_b=svm.mu_b 
    end
    
    if isnan(c)
        c=svm.c
    end
    
    # Create DataFrame
    res = DataFrame()
    res[:p] = p
    
    # Filter DataFrame
    fdf = df[abs.(df[:p] .- p) .< 1e-4, :]
    
    # Interpolate Values
    eq_deriv_fun = Dierckx.Spline1D(fdf[:vb], fdf[:eq_deriv], k=3, bc="extrapolate")
    
    # Compute optimal VB:
    vbroots = roots(eq_deriv_fun; maxn=8)
    if !isempty(vbroots)
        # debt_interp = Dierckx.Spline1D(fdf[:vb], fdf[:debt], k=3, bc="extrapolate")
        # abs_debt_diff = abs.(debt_interp(vbroots) .- get_agg_p(svm, p=p))
        # res[:vb] = vbroots[argmin(abs_debt_diff)]
        eq_min_val_interp = Dierckx.Spline1D(fdf[:vb], fdf[:eq_min_val], k=3, bc="extrapolate")
        abs_eq_min_val = abs.(eq_min_val_interp(vbroots))
        res[:vb] = vbroots[argmin(abs_eq_min_val)]
    else
        ref_vbgrid = range(minimum(fdf[:vb]), stop=maximum(fdf[:vb]), length=N)
        res[:vb] = ref_vbgrid[argmin(abs.(eq_deriv_fun(ref_vbgrid)))]
    end
    
    # Equity Values
    res[:eq_deriv] = eq_deriv_fun(res[:vb])

    # Interpolate Functions
    interp_cols = vcat([:debt, :equity],
                       [:eq_deriv_min_val, 
                        :eq_min_val, :eq_vb])
    res = interp_values(res, fdf, :vb, interp_cols)
    
    # Debt Values
    res = debt_at_par_diffs(svm, res, :p)

    # Fixed Values
    return non_interp_values(svm, res)

    # # Share Values
    # res[:firm_value] = res[:debt] .+ res[:equity]
    # res[:leverage] = (res[:debt]./res[:firm_value]) .* 100
    # res[:ROE] = (res[:equity]./(svm.pm.V0 .- res[:debt]) .- 1) .* 100.
    
    # return res
end


function eq_fd_processing(bt, svm, df; mu_b=NaN, c=NaN, N=10^5)
    if isnan(mu_b)
        mu_b=svm.mu_b 
    end
    
    if isnan(c)
        c=svm.c
    end
    
    params_dict = Dict()
    for par in vcat(bt.dfc.main_params, [:mu_b, :m, :c], bt.dfc.fixed_params)
        params_dict[par] = unique(df[par])[1]
    end

    # Compute the Solutions for each pair (c, p)
    tmp = fetch(@spawn [eq_deriv_root_search(svm, df, p; mu_b=mu_b, c=c, N=N) 
                        for p in unique(df[:p])])

    # Collect results
    eqfinal = hcat(vcat(tmp...), repeat(DataFrame(params_dict), outer=size(tmp, 1)))

    # Rearrange columns
    return eqfinal[bt.dfc.dfcols]
end


function eq_fd_sol(bt, svm, df; N=10^5)
    res=DataFrame()
    
    # Find p for which |D - P| = 0:
    ddiff = Dierckx.Spline1D(df[:p], df[:debt_diff], k=3, bc="extrapolate")
    pgrid = range(minimum(df[:p]), stop=maximum(df[:p]), length=N)
    res[:p] = pgrid[argmin(abs.(ddiff(pgrid)))]
    res[:debt_diff] = ddiff(res[:p])
    
    # Interpolate 
    interp_cols = vcat([:vb, :debt, :equity],
                       [:debt_per_diff],
                       [:eq_deriv, :eq_deriv_min_val, 
                        :eq_min_val, :eq_vb])
    res = interp_values(res, df, :p, interp_cols)
    res = non_interp_values(svm, res)
    
    # Get Non-Variable Values    
    cols = vcat(bt.dfc.main_params, [:mu_b, :m, :c], bt.dfc.fixed_params)
    for col in cols
        res[col] = df[1, col]
    end
    
    return res[bt.dfc.dfcols]
end