
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


function find_full_info_vb(svm, jks;
                           lb::Float64=.75, ub::Float64=1.25,
                           vbN::Int64=15, N::Int64=10^5)
    if get_obj_model(svm) == "cvm"
        return get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=jks.mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    else
        vbl = get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=jks.mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    
        vbh = get_cvm_vb(svm, svm.pm.sigmal;
                         mu_b=jks.mu_b, m=jks.m,
                         c=jks.c, p=jks.p)

        vbgrid = range(.75 * minimum([vbl, vbh]),
                       stop=1.25 * maximum([vbl, vbh]), length=vbN)

        res = @time fetch(@spawn [eq_fd(svm; vbl=vbl, mu_b=jks.mu_b,
                                        m=jks.m, c=jks.c, p=jks.p)
                                  for vbl in vbgrid])

        eqdf = full_info_eq_deriv_root_search(svm, vcat(res...); N=N)

        return eqdf[1, :vb]
    end
end


function get_joint_k_struct!(jf;
                             jks=JointKStruct(fill(NaN, 10)...),
                             mu_b::Float64=NaN,
                             m::Float64=NaN,
                             c::Float64=NaN,
                             p::Float64=NaN)

    # if !isnan(mu_s)
    #     jks.mu_s = mu_s
    # else isnan(jks.mu_s)
    #     jks.mu_s = jf.jks.mu_s
    # end
    
    if !isnan(mu_b)
        setfield!(jks, :mu_b, mu_b)
    elseif isnan(jks.mu_b)
        setfield!(jks, :mu_b, jf.jks.mu_b)
    end

    if !isnan(m)
        setfield!(jks, :m, m)
    elseif isnan(jks.m)
        setfield!(jks, :m, jf.jks.m)
    end
    
    if !isnan(c)
        setfield!(jks, :c, c)
    elseif isnan(jks.c)
        setfield!(jks, :c, jf.jks.c)
    end

    if !isnan(p)
        setfield!(jks, :p, p)
    elseif isnan(jks.p)
        setfield!(jks, :p, jf.jks.p)
    end
   
    return jks 
end


function set_full_information_vb!(jf, jks;
                                  rerun_fi_vb::Bool=false,
                                  fi_sf_vb::Float64=NaN,
                                  fi_rf_vb::Float64=NaN,
                                  lb::Float64=.75,
                                  ub::Float64=1.25,
                                  vbN::Int64=20)
    
    if any([.&(isnan(jks.fi_sf_vb), isnan(fi_sf_vb)), rerun_fi_vb])
        fi_sf_vb = find_full_info_vb(jf.sf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_sf_vb, fi_sf_vb)
    elseif !isnan(fi_sf_vb)
         setfield!(jks, :fi_sf_vb, fi_sf_vb)       
    end

    if any([.&(isnan(jks.fi_rf_vb), isnan(fi_rf_vb)), rerun_fi_vb])
        fi_rf_vb = find_full_info_vb(jf.rf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_rf_vb, fi_rf_vb)
    elseif !isnan(fi_rf_vb)
         setfield!(jks, :fi_rf_vb, fi_rf_vb)       
    end

    return jks
end

function joint_eq_set_k_struct!(jf;
                                jks=JointKStruct(fill(NaN, 10)...),
                                mu_s::Float64=NaN,
                                mu_b::Float64=NaN,
                                m::Float64=NaN,
                                c::Float64=NaN,
                                p::Float64=NaN,
                                rerun_fi_vb::Bool=false,
                                fi_sf_vb::Float64=NaN,
                                sf_vb::Float64=NaN,
                                fi_rf_vb::Float64=NaN,
                                rf_vb::Float64=NaN,
                                lb::Float64=.75,
                                ub::Float64=1.25,
                                vbN::Int64=20)
    
    jks = get_joint_k_struct!(jf; jks=jks,
                              mu_b=mu_b,
                              m=m, c=c, p=p)

    # Default Barriers ##############################
    # Full Information: fi_sf_vb, fi_rf_vb
    jks = set_full_information_vb!(jf, jks;
                                   rerun_fi_vb=rerun_fi_vb,
                                   fi_sf_vb=fi_sf_vb,
                                   fi_rf_vb=fi_rf_vb,
                                   lb=lb, ub=ub,
                                   vbN=vbN)
    
    # setfield!(jks, :sf_vb, maximum(x->isnan(x) ? -Inf : x, [sf_vb, jks.fi_sf_vb]))

    # rf_vb = maximum([minimum(x->isnan(x) ? Inf : x, [rf_vb, jks.fi_rf_vb]),
    #                  jks.sf_vb])
    # setfield!(jks, :rf_vb, minimum(x->isnan(x) ? Inf : x, [rf_vb, jks.fi_rf_vb]))

    jks.sf_vb = !isnan(sf_vb) ? sf_vb : jks.fi_sf_vb
    jks.rf_vb = !isnan(rf_vb) ? rf_vb : jks.fi_rf_vb
    # Joint Equilibrium Barrier
    setfield!(jks, :vbl, maximum([jks.sf_vb, jks.rf_vb]))
    # ###############################################

    # Measure of Safe Firms
    if mu_s < 0.
        println("Setting mu_s to zero!")
        setfield!(jks, :mu_s, 0.0)
    elseif mu_s > 1.
        println("Setting mu_s to 1.!")
        setfield!(jks, :mu_s, 1.)
    elseif !isnan(mu_s)
        setfield!(jks, :mu_s, mu_s)
    end
    
    return jks
end


function get_type_contingent_vbs(vbl::Float64, 
                                 fi_sf_vb::Float64, 
                                 fi_rf_vb::Float64;
                                 sf_defaults_first::Bool=true)
    
    sf_vb = fi_sf_vb
    rf_vb = fi_rf_vb
    if vbl > maximum([fi_sf_vb, fi_rf_vb])
        sf_vb = sf_defaults_first ? vbl : fi_sf_vb
        rf_vb = !sf_defaults_first ? vbl : fi_rf_vb
    elseif  vbl < minimum([fi_sf_vb, fi_rf_vb])
        sf_vb = rf_vb = vbl
    elseif fi_sf_vb <= fi_rf_vb
        sf_vb = sf_defaults_first ? vbl : fi_sf_vb
        rf_vb = vbl
    else # fi_sf_vb > vbl > fi_rf_vb
        sf_vb = vbl
        rf_vb = !sf_defaults_first ? vbl : fi_rf_vb
    end
     
    
    return sf_vb, rf_vb
end


function interp_optimal_vbs(jf, jks, df::DataFrame; 
                            sf_defaults_first::Bool=true, 
                            vbN::Int64=10^5,
                            spline_k::Int64=3,
                            spline_bc::String="extrapolate")

    # Form Refined Grid of Unique vbl values
    vb_grid = range(minimum(df[:vbl]), stop=maximum(df[:vbl]), length=vbN)
    
    # Interpolate Equity and Equity Derivative Functions
    tmp = Dict()
    for var in [:rf_eq_deriv, :sf_eq_deriv] 
        tmp[var] = Dierckx.Spline1D(df[:vbl], df[var]; k=spline_k, bc=spline_bc)
    end

    # Find Optimal VBs
    opt_sf_vb = vb_grid[argmin(abs.(tmp[:sf_eq_deriv](vb_grid)))]
    opt_rf_vb = vb_grid[argmin(abs.(tmp[:rf_eq_deriv](vb_grid)))]

    return opt_sf_vb, opt_rf_vb
end


