
function get_cvm_vb_funs(svm, c, p_grid)
    cvm_vbl = [AnalyticFunctions.get_cvm_vb(svm, svm.pm.sigmal; c=c, p=p) for p in p_grid]
    cvm_vbh = [AnalyticFunctions.get_cvm_vb(svm, svm.pm.sigmah; c=c, p=p) for p in p_grid]
    cvm_vbl_fun = Dierckx.Spline1D(p_grid, cvm_vbl, k=3, bc="extrapolate")
    cvm_vbh_fun = Dierckx.Spline1D(p_grid, cvm_vbh, k=3, bc="extrapolate")

    return cvm_vbl_fun, cvm_vbh_fun
end

# #######################################################
# ##################### DEBT AT PAR #####################
# #######################################################
function debt_at_par_grid_bounds(svm, c, x_min, x_max; p=NaN, vbl=NaN,
                                 f_low=.75, f_high=1.25)
    xgrid = range(f_low * x_min, stop=f_high * x_max, length=50)

    
    vbh = get_cvm_vb(svm, svm.pm.sigmah; c=c, p=p)
    if isnan(vbl)
        vbhl_vec = [vbh/xval for xval in xgrid]
    elseif isnan(p)
        vbhl_vec = [get_cvm_vb(svm, svm.pm.sigmah; 
                               c=c, p=x)/vbl for x in xgrid]
    else
        println("Please enter a value for p or vbl")
        return
    end

    # Interpolate
    vbhlf = Dierckx.Spline1D(xgrid, vbhl_vec, k=3, bc="extrapolate")

    # Refine
    xgrid_ref  = range(f_low * x_min, stop=f_high * x_max, length=10^5)
    vbhl_ref = [vbhlf(xval) for xval in xgrid_ref]
    
    # Bounds
    x0 = xgrid_ref[argmin([abs(vbhl - 1.025 * svm.bi.vbhlmin) for vbhl in vbhl_ref])]
    x1 = xgrid_ref[argmin([abs(vbhl - .975 * svm.bi.vbhlmax) for vbhl in vbhl_ref])]

    
    # Now Make Sure vt bounds are satisfied:
    # xmax = minimum([maximum([x0,x1]),  .975 * svm.pm.V0 * exp(- minimum(svm.bs.vtgrid))])
    
    return minimum([x0,x1, .95 * svm.pm.V0]), minimum([maximum([x0,x1]), svm.pm.V0]) # 1.001 * minimum([x0, x1]), .999 * maximum([x0, x1])
end


function debt_at_par_check_conds(svm, df, xvar, xmin, xmax; disp_option=false)
    # DF should satisfy the following conditions:
    # 1. At least 5 (P, VB) pairs with precision <= 2*1e-2 (2% off-mark)
    # 2. At least 2 out of 5 (P, VB) pairs should be in the interval:
    #                       [.9 * p_min, 1.1 * p_max)]
    # 3. At least 2 out of 5 pairs should have precision <= 1e-2
    cond1 = sum(df[:abs_debt_per_diff] .<= 2e-2) >= 5
    cond2 = sum(.&(df[xvar] .>= xmin * .9, df[xvar] .<= xmax * 1.1)) >= 2
    cond3 = sum(df[:abs_debt_per_diff] .<= 1e-2) >= 2
    conds = cond1 & cond2 & cond3

    if disp_option
        println(string("Condition 1: ", cond1))
        println(string("Condition 2: ", cond2))
        println(string("Condition 3: ", cond3))
        println(string("Conds: ", conds))
    end
    
    return conds
end


function debt_at_par_diffs(svm, df, xvar; p=NaN, sort_vars=[])
    if xvar ==:p
        aggP = [get_agg_p(svm, p=p) for p in df[:p]]
    elseif !isnan(p)
        aggP = get_agg_p(svm, p=p)
    else
        println("Missing p value(s). Exiting...")
    end

    df[:debt_diff] = df[:debt] .- aggP
    df[:abs_debt_diff] = abs.(df[:debt_diff])
    df[:debt_per_diff] = (df[:debt_diff]) ./ aggP
    df[:abs_debt_per_diff] = abs.(df[:debt_per_diff])

    if isempty(sort_vars)
        return sort!(df, xvar)
    else
        return sort!(df, sort_vars)
    end
end


function debt_at_par_new_df(svm, c, p, df, xvar, debt_vars=debt_vars)
    xgrid_new = (Array(df[xvar])[2:end] + Array(df[xvar])[1:end-1])./ 2.
    df_new = DataFrame(hcat(reverse(xgrid_new), 999. * (ones(length(xgrid_new), size(df, 2) - 1))),
                       vcat([xvar, :debt], debt_vars))
    df_new[:debt] = fetch(@spawn [get_svm_debt_price(svm, vbl; c=c, p=p)
                                  for vbl in df_new[xvar]])
    
    return sort!(vcat(df, df_new), :vb)
end


function debt_at_par_main(svm, c, df, xvar, xgrid_ref;
                          vbl=NaN, p=NaN, debt_vars=debt_vars, disp_option=true)
    # Notice here I use the percentage difference, not
    # the absolute percentage difference. I want to see
    # where the function crosses the zero line.
    
    # Interpolate
    if xvar == :p
        agg_p_vec = [get_agg_p(svm, p=p) for p in df[:p]]
        diff_interp_fun = Dierckx.Spline1D(df[:p],
                              (df[:debt] .- agg_p_vec) ./ agg_p_vec,
                                           k=3, bc="extrapolate")
    else
        diff_interp_fun = Dierckx.Spline1D(df[:vb],
                              (df[:debt] .- get_agg_p(svm, p=p)) ./ get_agg_p(svm, p=p),
                                           k=3, bc="extrapolate")
    end

    diffs = diff_interp_fun(xgrid_ref)
    
    # Find Candidate for the Global Minimum:
    xroot = reverse(xgrid_ref)[argmin(abs.(reverse(diffs)))]
    if disp_option
        println(string(xvar, " root: ", xroot))
    end
            
    if sum(abs.(df[xvar] .- xroot) .< 1e-4) == 0
        if xvar == :p
            debt_root = get_svm_debt_price(svm, vbl; p=xroot)
        else
            debt_root = get_svm_debt_price(svm, xroot; c=c, p=p)
        end

        sort!(push!(df, [xroot, debt_root, NaN, NaN, NaN, NaN]),
              vcat([xvar, :debt], debt_vars))
    end
    
    return debt_at_par_diffs(svm, df, xvar; p=p) 
end


# ATTENTION: NEED TO ADJUST P TO ACCOUNT FOR DIFFERENT M!!!!
# ATTENTION: p_min = np.min([cvml_p,cvmh_p])
# ATTENTION: p_max = 1.5 * np.max([cvml_p,cvmh_p])
function svm_debt_at_par(svm, x_min, x_max, c; p=NaN, vbl=NaN,
                         N1=15, N2=10^5, f_low=.85, f_high=1.15,
                         debt_vars=debt_vars, disp_option=true)

    start_tic = time_ns()

    # Initial Check:
    if isnan(p) & isnan(vbl)
        println("Please enter a value for p or vbl")
        return
    end

    # Form Grids:
    xvar = [:p, :vb][isnan.([p, vbl])][1]
    fixed_var = [:p, :vb][isnan.([p, vbl]).==false][1]
    fixed_var_value = [p, vbl][isnan.([p, vbl]).==false][1]
    xmin, xmax = debt_at_par_grid_bounds(svm, c, x_min, x_max; 
                                         p=p, vbl=vbl, f_low=.75, f_high=1.25)
    xgrid = range(xmin, stop=xmax, length=N1)
    xgrid_ref = range(xmin, stop=xmax, length=N2)

    # Compute Debt Values
    df = DataFrame(hcat(reverse(xgrid), 999. * ones(length(xgrid), 5)),
                        vcat([xvar, :debt], debt_vars))

    if xvar == :p
        df[:debt] = fetch(@spawn [get_svm_debt_price(svm, vbl; c=c, p=pval) for pval in df[:p]])
    else
        df[:debt] = fetch(@spawn [get_svm_debt_price(svm, vblval; c=c, p=p) for vblval in df[:vb]])
    end

    # Sort DataFrame
    sort!(df, xvar)

    # Compute Candidate Solution:
    df = debt_at_par_main(svm, c, df, xvar, xgrid_ref,
                          vbl=vbl, p=p, disp_option=disp_option)

    # Check Conditions:
    conds = debt_at_par_check_conds(svm, df, xvar, x_min, x_max,
                                    disp_option=disp_option)
    
    count = 0
    if disp_option
        println(string("Counter: ", count))
        println("Unique ", xvar, " values: ", length(unique(df[xvar])))
    end

    while !conds & (length(unique(df[xvar])) < 20) & (count < 10)
        if disp_option
            println("While loop")
        end

        # Add Data Points:
        df = debt_at_par_new_df(svm, c, p, df, xvar)

        # Compute Candidate Solution:
        df = debt_at_par_main(svm, c, df, xvar,
                              xgrid_ref, vbl=vbl, p=p)
        
        # Check Conditions:
        conds = debt_at_par_check_conds(svm, df, xvar, x_min, x_max,
                                        disp_option=disp_option)

        # Update Counter:
        count += 1

        if disp_option
            println(string("Counter: ", count))
            println("Unique ", xvar, " values: ", length(unique(df[xvar])))
        end
    end

    if disp_option
        println("Debt at par: preparing results...")
    end
    df[:c] = c
    df[fixed_var] = fixed_var_value
    df[:count] = count
    df[:time] = time_ns() - start_tic
                                
    cols = vcat([:c, fixed_var, xvar, :debt],
                debt_vars, [:count, :time])
    df = df[cols]

    if disp_option
        println("returning results...")
        println(string("length: ", size(df, 1)))                            
    end

    return df
end
# #######################################################
