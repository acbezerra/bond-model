
function joint_eq_fd_newly_issued_bonds(jf, jks, vbl::Float64,
                                        vgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                            Base.TwicePrecision{Float64}};
                                        vtN::Int64=10^3,
                                        sf_ftype::String="bf",
                                        rf_ftype::String="bf",
                                        spline_k::Int64=3,
                                        spline_bc::String="extrapolate")

    # Common Payoffs ###############################
    rfbond = rfbond_price(jks.m, jks.c, jks.p,
                          jf.sf.pm.r, jf.sf.pm.xi,
                          jf.sf.pm.kappa)
    
    #dpayoff = on_default_payoff(0., vbl, jks.m, jks.mu_b,
    #                            jks.m, jks.c, jks.p,
    #                            jf.sf.pm.r, jf.sf.pm.xi,
    #                            jf.sf.pm.kappa, jf.sf.pm.alpha)
    # ##############################################

    
    _, v_subgrid = grid_creator((1 + 1e-4) * minimum(jf.sf.bs.vtgrid), maximum(vgrid), vtN)

    # ##############################################
    # Whether to compute bond pricing surfaces
    # (fixed maturity) on the fly
    if any([abs.(jf.sf.m - jks.m) > 1e-4, sf_ftype == "bft"])
        sf_ftype = "bft"
        jf.sf = bpr_interp_fixed_ttm(jf.sf; ttm=jks.m)
    end

    if any([abs.(jf.rf.m - jks.m) > 1e-4, rf_ftype == "bft"])
        rf_ftype = "bft"
        jf.rf = bpr_interp_fixed_ttm(jf.rf; ttm=jks.m)
    end
    # ###############################################

    bpr_vec = fetch(@spawn [joint_bond_price(jf, jks.m; vt=v,
                                             jks=jks,
                                             mu_s=jks.mu_s,
                                             sf_ftype=sf_ftype,
                                             rf_ftype=rf_ftype) for v in v_subgrid])
    
    #bpr = Dierckx.Spline1D(vcat(.0, v_subgrid), vcat(dpayoff, bpr_vec); k=3, bc="extrapolate")
    bpr = Dierckx.Spline1D(v_subgrid, bpr_vec; k=spline_k, bc=spline_bc)

    return Array([minimum([bpr(v)[1], rfbond]) for v in vgrid])[2:end-1]
end


function adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)
    # Capital Structure
    for fn in fieldnames(JointKStruct)
        if isnan(getfield(jks, fn))
            setfield!(jks, fn, jf.jks)
        end
    end

    # Measure of Safe Firms
    if isnan(mu_s)
        mu_s = jks.mu_s
    elseif mu_s < 0.
        println("Setting mu_s to zero!")
        mu_s = 0.0
    elseif mu_s > 1.
        println("Setting mu_s to 1.!")
        mu_s = 1.
    end
    

    # Set Asset Value
    if .&(isnan(Vt), isnan(vt))
        Vt = get_param(jf.sf, :V0)
        vt = log(Vt/jks.vbl)
    elseif isnan(vt)
        vt = log(Vt/jks.vbl)
    end
    
    return jks, mu_s, vt
end


function joint_bond_price(jf, ttm::Float64;
                          jks=JointKStruct(fill(NaN, 10)...),
                          mu_s::Float64=NaN,
                          vt::Float64=NaN, Vt::Float64=NaN,
                          sf_ftype::String="bf",
                          rf_ftype::String="bf")

    jks, mu_s, vt = adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)

    # Compute Bond Prices
    if get_obj_model(jf.sf) == "cvm" 
        sf_bpr = get_cvm_bond_price(jf.sf, ttm, jf.sf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt),
                                    vb=jks.sf_vb,
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p)
    else
        sf_bpr = get_svm_bond_price(jf.sf, jks.sf_vb, ttm;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ftype=sf_ftype)
    end

    if get_obj_model(jf.rf) == "cvm" 
        rf_bpr = get_cvm_bond_price(jf.rf, ttm, jf.rf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt),
                                    vb=jks.rf_vb,
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p)
    else   
        rf_bpr = get_svm_bond_price(jf.rf, jks.rf_vb, ttm;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ftype=rf_ftype)
    end
 
    # Joint Price
    return mu_s * sf_bpr + (1. - mu_s) * rf_bpr
end


function joint_debt_price(jf;
                          jks=JointKStruct(fill(NaN, 10)...),
                          mu_s::Float64=NaN,
                          vt::Float64=NaN,
                          Vt::Float64=NaN,
                          ttmN0::Int64=10^2,
                          ttmN::Int64=10^4)

    jks, mu_s, vt = adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)

    
    # Compute Debt Price #################################################
    # either model is CVM or bond surface maturity
    # must coincide with capital structure maturity.
    if get_obj_model(jf.sf) == "cvm" 
        sf_dpr = get_cvm_debt_price(jf.sf, jks.sf_vb,
                                    jf.sf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p,
                                    N1=ttmN0, N2=ttmN)
    elseif abs.(jf.sf.m - jks.m) < 1e-4
        sf_bpr = get_svm_debt_price(jf.sf, jks.sf_vb;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ttmN=ttmN)
    else
        println(string("Capital Structure Maturity of Safe Type ",
                       "does not coincide with ",
                       "m used in bond pricing surfaces computation!"))
        return NaN
    end

    if get_obj_model(jf.rf) == "cvm" 
        rf_dpr = get_cvm_debt_price(jf.rf, jks.rf_vb,
                                    jf.rf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p,
                                    N1=ttmN0, N2=ttmN)
    elseif abs.(jf.rf.m - jks.m) < 1e-4
        rf_dpr = get_svm_debt_price(jf.rf, jks.rf_vb;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ttmN=ttmN)
    else
        println(string("Capital Structure Maturity of Risky Type ",
                       "does not coincide with ",
                       "m used in bond pricing surfaces computation!"))
        return NaN
    end

    # Joint Price
    return mu_s * sf_dpr + (1 - mu_s) * rf_dpr
end
