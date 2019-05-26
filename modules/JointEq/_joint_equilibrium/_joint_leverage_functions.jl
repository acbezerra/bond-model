
function find_joint_payoffs(jf, jks, mu_b_grid;
                            sf_obj_fun::Symbol=:firm_value,
                            rf_obj_fun::Symbol=:MBR,
                            spline_k::Int64=3,
                            spline_bc::String="extrapolate",
                            N::Int64=10^5)


    dfl = @time fetch(@spawn [find_joint_optimal_vb(jf, jks;
                                                    mu_b=mu_b,
                                                    rerun_fi_vb=true)
                              for mu_b in mu_b_grid])
    
    # Store results in a DataFrame
    df = vcat(dfl...)
    sf_df = df[isnan.(df[:rf_vb]), :]
    rf_df = df[isnan.(df[:sf_vb]), :]
    
    
    # Interpolate Objective Functions in mu_b ################################
    sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    rf_objf = Dierckx.Spline1D(rf_df[:mu_b], rf_df[rf_obj_fun]; k=spline_k, bc=spline_bc)

    # Refine mu_b_grid
    ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=N)

    return sf_objf, rf_objf, ref_mu_b_grid
end




function find_mu_b_intervals(ref_mu_b_grid::StepRangeLen{Float64,
                                                         Base.TwicePrecision{Float64},
                                                         Base.TwicePrecision{Float64}},
                             cond::BitArray{1}; N::Int64=20)
    #ma = ref_mu_b_grid[cond]
    #md = abs.(ma[2:end] - ma[1:end-1])
    mu_b_discarded = ref_mu_b_grid[cond .== false]
    
   
    #maximum(md) - sum(md)/size(md, 1) < 1e-6
    mu_b_dis_min = minimum(mu_b_discarded)
    mu_b_dis_max = maximum(mu_b_discarded)

    # Check if there is only one interval
    if any([abs.(mu_b_dis_min - minimum(ref_mu_b_grid)) < 1e-6,
            abs.(mu_b_dis_max - maximum(ref_mu_b_grid)) < 1e-6])
        filtered_mu_b_grid_1 = range(minimum(ref_mu_b_grid[cond]),
                                     stop=maximum(ref_mu_b_grid[cond]),
                                     length=N)
        filtered_mu_b_grid_2 = range(0, stop=0, length=0)
    else
        # int1_up  = ref_mu_b_grid[argmax(md) + 1]
        # int2_down = reverse(ref_mu_b_grid)[argmax(reverse(md)) + 1]
        # filtered_mu_b_grid_1 = range(minimum(mu_b_grid), stop=int1_up, length=10)
        # filtered_mu_b_grid_2 = range(int2_down, stop=maximum(mu_b_grid), length=10)
        
        # In case there are two intervals:
        filtered_mu_b_grid_1 = range(minimum(ref_mu_b_grid[cond]),
                                     stop=minimum(mu_b_discarded),
                                     length=N)
        filtered_mu_b_grid_2 = range(maximum(mu_b_discarded),
                                     stop=maximum(ref_mu_b_grid[cond]),
                                     length=N)
    end

    return filtered_mu_b_grid_1, filtered_mu_b_grid_2
end


function find_opt_mu_b(sf_objf,
                       filtered_mu_b_grid_1::StepRangeLen{Float64,
                                                          Base.TwicePrecision{Float64},
                                                            Base.TwicePrecision{Float64}},
                       filtered_mu_b_grid_2::StepRangeLen{Float64,
                                                          Base.TwicePrecision{Float64},
                                                          Base.TwicePrecision{Float64}})
    
    # Maximize Safe Firm's Objective Function
    mu_b_opt_1 = filtered_mu_b_grid_1[argmax(sf_objf(filtered_mu_b_grid_1))]
    
    if size(filtered_mu_b_grid_2, 1) == 0
        return mu_b_opt_1
    else
        # Maximize Safe Firm's Objective Function
        mu_b_opt_2 = filtered_mu_b_grid_2[argmax(sf_objf(filtered_mu_b_grid_2))]
        
        return maximum([mu_b_opt_1, mu_b_opt_2])
    end
        
        
end


    # # Now compute payoffs for the safe firm 
    # s_eqdf_list = fetch(@spawn [sep_eq_fd(jf.sf, deepcopy(jks), mu_b) 
    #                             for mu_b in filtered_mu_b_grid])
    # sf_df = vcat(s_eqdf_list...)

    # sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    
    # # Maximize Safe Firm's Objective Function
    # ref_mu_b_grid = range(minimum(filtered_mu_b_grid),
    #                         stop=maximum(filtered_mu_b_grid), length=N2)
    # mu_b_opt = ref_mu_b_grid[argmax(sf_objf(ref_mu_b_grid))]


function sep_eq_fd(sf, jks, mu_b::Float64)
    sf_vb = find_full_info_vb(sf, jks; mu_b=mu_b)
    return eq_fd(sf, vbl=sf_vb, 
                 mu_b=mu_b, m=jks.m, 
                 c=jks.c, p=jks.p)
    
end


function separate_eq_calculator(jf, jks, mu_b_opt::Float64, fi_rf_mu_b::Float64)
    s_eqdf = sep_eq_fd(jf.sf, deepcopy(jks), mu_b_opt)
    r_eqdf = sep_eq_fd(jf.rf, deepcopy(jks), fi_rf_mu_b)

    s_eqdf[:fi_vb] = s_eqdf[1, :vb]
    s_eqdf[:sf_vb] = s_eqdf[1, :vb]
    s_eqdf[:rf_vb] = NaN
    
    r_eqdf[:fi_vb] = r_eqdf[1, :vb]
    r_eqdf[:sf_vb] = NaN 
    r_eqdf[:rf_vb] = r_eqdf[1, :vb]

    eqdf = vcat(s_eqdf, r_eqdf)
    eqdf[:sf_defaults_first] = s_eqdf[1, :vb] > r_eqdf[1, :vb]
    eqdf[:eq_type] = "separating"
    eqdf[:mu_s] = jks.mu_s

    return eqdf
end

