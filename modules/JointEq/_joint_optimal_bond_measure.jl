
function bond_measure_full_info_vb(svm, jks, mu_b::Float64;
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5)

    jks2 = deepcopy(jks)
    setfield!(jks2, :mu_b, mu_b)

    return find_full_info_vb(svm, jks2;
                             lb=.75, ub=ub,
                             vbN=15, N=vbN2)
end


function find_optimal_bond_measure(svm;
                                   mu_b_min::Float64=.5,
                                   mu_b_max::Float64=15.,
                                   mu_bN::Int64=20,
                                   mu_bN2::Int64=10^5,
                                   jks=JointKStruct(fill(NaN, 10)...),
                                   m::Float64=NaN,
                                   c::Float64=NaN,
                                   p::Float64=NaN,
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5,
                                   spline_k::Int64=3,
                                   spline_bc::String="extrapolate")

    # Set Capital Structure #########################
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
    # ###############################################

    # Form mu_b grid
    mu_b_grid = range(mu_b_min, stop=mu_b_max; length=mu_bN)

    # Get Optimal VB for each mu_b value
    if get_obj_model(svm) == "cvm"
        vbl_list = [get_cvm_vb(svm, svm.pm.sigmal;
                               mu_b=mu_b, m=jks.m, c=jks.c, p=jks.p)
                    for mu_b in mu_b_grid]
        
    else
        vbl_list = fetch(@spawn [bond_measure_full_info_vb(svm, jks, mu_b;
                                                           lb=lb, ub=ub,
                                                           vbN=vbN, vbN2=vbN2)
                                 for mu_b in mu_b_grid])
    end

    # Run Equity Finite Differences Method ###############
    tmp = DataFrame(mu_b = mu_b_grid, vbl = vbl_list)
    res = fetch(@spawn [eq_fd(svm; vbl=tmp[r, :vbl],
                              mu_b=tmp[r, :mu_b],
                              m=jks.m, c=jks.c, p=jks.p)
                        for r in 1:size(tmp, 1)])
    df = vcat(res...)

    # Interpolate VB and Firm Value in mu_b
    funs = Dict()
    for var in [:vb, :firm_value]
        funs[var] = Dierckx.Spline1D(df[:mu_b], df[var];
                                     k=spline_k, bc=spline_bc)
    end
    

    # Form mu_b refined grid
    mu_b_grid_ref = range(minimum(df[:mu_b]),
                          stop=maximum(df[:mu_b]), length=mu_bN2)

    # Get optimal mu_b
    opt_mu_b = mu_b_grid_ref[argmax(funs[:firm_value](mu_b_grid_ref))]

    # Set optimal mu_b
    setfield!(jks, :mu_b, opt_mu_b)

    # Get Optimal VB
    setfield!(jks, :vbl, funs[:vb](opt_mu_b))

    return eq_fd(svm; vbl=jks.vbl,
                 mu_b=jks.mu_b,
                 m=jks.m, c=jks.c, p=jks.p)
end
