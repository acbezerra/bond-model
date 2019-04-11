

function get_cvm_p_debt_diff(svm, sigma, mu_b, c, pgrid; N=10000)
    # Compute CVM VBs
    cvm_vb = [get_cvm_vb(svm, sigma; mu_b=mu_b, c=c, p=p) for p in pgrid]
    
    # Interpolate CVM VBs
    cvm_vb_fun = Dierckx.Spline1D(pgrid, cvm_vb, k=3, bc="extrapolate")
    
    # Compute Debt Values
    cvm_debt = fetch(@spawn [get_cvm_debt_price(svm, cvm_vb_fun(p), sigma; 
                                            mu_b=mu_b, c=c, p=p) for p in pgrid])
    
    # Compute Debt-P
    aggPgrid = [get_agg_p(svm; mu_b=mu_b, p=p) for p in pgrid]
    cvm_debt_diff = cvm_debt - aggPgrid
    
    # Interpolate and Find Optimal Value
    pgrid_ref = range(pgrid[1], stop=pgrid[end], length=N)
    cvm_dff = Dierckx.Spline1D(pgrid, cvm_debt_diff, k=3, bc="extrapolate")
    # Get highest pj8 value:
    cvm_pOpt = reverse(pgrid_ref)[argmin(abs.(reverse(cvm_dff(pgrid_ref))))]
    cvm_vbOpt = cvm_vb_fun(cvm_pOpt)
    
    return Dict("cvm_vb_fun" => cvm_vb_fun,
                "cvm_debt" => cvm_debt,
                "cvm_dff" => cvm_dff,
                "cvm_pOpt" => cvm_pOpt,
                "cvm_vbOpt" => cvm_vbOpt)
end


function p_candidates(svm; pmin=5, pmax=NaN, mu_b=NaN, c=NaN, N1=20, N2=10000, lower_adj=.85, upper_adj=1.15)
    if isnan(pmax)
        pmax = .9 * svm.pm.V0
    end

    if isnan(mu_b)
        mu_b = svm.mu_b
    end
    
    if isnan(c)
        c = svm.c
    end
        
    # Form p grids
    pgrid = range(pmin, stop=pmax, length=N1)

    # Get VB, Debt, Debt Diff, Debt at Par values
    cvml_dict = get_cvm_p_debt_diff(svm, svm.pm.sigmal, mu_b, c, pgrid; N=N2)
    cvmh_dict = get_cvm_p_debt_diff(svm, svm.pm.sigmah, mu_b, c, pgrid; N=N2)
    
    # Form Adjusted pgrid
    p_lower = lower_adj * minimum([cvml_dict["cvm_pOpt"], cvmh_dict["cvm_pOpt"]])
    p_upper = upper_adj * maximum([cvml_dict["cvm_pOpt"], cvmh_dict["cvm_pOpt"]])
    adj_pgrid = range(p_lower, stop=p_upper, length=N1)
   
    # Ensure VB values are within bounds
    
    return Dict("cvml" => cvml_dict, "cvmh" => cvmh_dict, "adj_pgrid" => adj_pgrid)
end


function bounded_vbs(svm, cont, p)
    
    cvm_vbl = cont["cvml"]["cvm_vb_fun"](p)
    cvm_vbh = cont["cvmh"]["cvm_vb_fun"](p)
    
    # Form VB grid
    vblmin = .5 * minimum([cvm_vbl, cvm_vbh])
    vblmax = 1.5 * maximum([cvm_vbl, cvm_vbh])
    vblgrid = range(vblmin, stop=vblmax, length=20)
    
    # Compute vbh/vbl ratios
    ratio = [cvm_vbh/vbl for vbl in vblgrid]

    # Interpolate ratios
    ratiof = Dierckx.Spline1D(vblgrid, ratio, k=3, bc="extrapolate")

    # Refine vblgrid
    vblgrid_ref = range(vblmin, stop=vblmax, length=10000)
    
    # Ensure vbl values satisfy bhl_min <= cvm_vbh/vbl <= bhl_max
    bounded_vblgrid = [vb for vb in vblgrid_ref if 
                        (ratiof(vb) <= svm.bi.vbhlmax) & (ratiof(vb) >= svm.bi.vbhlmin)]
    vblmin = minimum(bounded_vblgrid)
    vblmax = maximum(bounded_vblgrid)
                
    return vblmin, vblmax
end