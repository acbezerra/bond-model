
function get_pool_eqdf(jf, jks, mu_b_grid,
                       fi_rf_obj_val::Float64;
                       sf_obj_fun::Symbol=:firm_value,
                       rf_obj_fun::Symbol=:MBR,
                       spline_k::Int64=3,
                       spline_bc::String="extrapolate",
                       mu_bN2::Int64=10^5)
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
    ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=mu_bN2)

    cond = rf_objf(ref_mu_b_grid) .>= fi_rf_obj_val
    filtered_mu_b_grid = ref_mu_b_grid[cond]

    # Maximize Safe Firm's Objective Function
    mu_b_opt = filtered_mu_b_grid[argmax(sf_objf(filtered_mu_b_grid))]
    
    
    # Compute Results
    return find_joint_optimal_vb(jf, jks;
                                 mu_b=mu_b_opt, rerun_fi_vb=true)
    
end


function sep_eq_fd(sf, jks, mu_b)
    sf_vb = find_full_info_vb(sf, jks; mu_b=mu_b)
    return eq_fd(sf, vbl=sf_vb, 
                 mu_b=mu_b, m=jks.m, 
                 c=jks.c, p=jks.p)
    
end


function get_sep_eqdf(jf, jks, mu_b_grid,
                      fi_rf_mu_b::Float64,
                      fi_rf_obj_val::Float64;
                      sf_obj_fun::Symbol=:firm_value,
                      rf_obj_fun::Symbol=:MBR,
                      spline_k::Int64=3,
                      spline_bc::String="extrapolate",
                      mu_bN::Int64=20,
                      mu_bN2::Int64=10^5)

    sep_misrep_jks = deepcopy(jks)
    sep_misrep_jks.mu_s = 1.

    dfl = @time fetch(@spawn [find_joint_optimal_vb(jf, sep_misrep_jks;
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
    ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=mu_bN2)

    # Filter
    cond = rf_objf(ref_mu_b_grid) .<= fi_rf_obj_val
    # filtered_mub_grid = ref_mu_b_grid[cond]
    filtered_mu_b_grid = range(minimum(ref_mu_b_grid[cond]),
                              stop=maximum(ref_mu_b_grid[cond]),
                              length=mu_bN)
    
    

    # Now compute payoffs for the safe firm 
    s_eqdf_list = fetch(@spawn [sep_eq_fd(jf.sf, deepcopy(jks), mu_b) 
                                for mu_b in filtered_mu_b_grid])
    sf_df = vcat(s_eqdf_list...)

    sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    
    # Maximize Safe Firm's Objective Function
    # mu_b_opt = ref_mu_b_grid[argmax(sf_objf(ref_mu_b_grid))]
    ref_mu_b_grid2 = range(minimum(filtered_mu_b_grid),
                            stop=maximum(filtered_mu_b_grid), length=mu_bN2)
    mu_b_opt = ref_mu_b_grid2[argmax(sf_objf(ref_mu_b_grid2))]


    s_eqdf = sep_eq_fd(jf.sf, deepcopy(jks), mu_b_opt)
    s_eqdf[:fi_vb] = s_eqdf[1, :vb]
    s_eqdf[:sf_vb] = s_eqdf[1, :vb]
    s_eqdf[:rf_vb] = NaN
    # s_eqdf[:firm_type] = "safe"
    
    r_eqdf = sep_eq_fd(jf.rf, deepcopy(jks), fi_rf_mu_b)
    r_eqdf[:fi_vb] = r_eqdf[1, :vb]
    r_eqdf[:sf_vb] = NaN 
    r_eqdf[:rf_vb] = r_eqdf[1, :vb]
    # r_eqdf[:firm_type] = "risky"

    
    eqdf = vcat(s_eqdf, r_eqdf)
    eqdf[:sf_defaults_first] = s_eqdf[1, :vb] > r_eqdf[1, :vb]
    eqdf[:eq_type] = "separating"
    eqdf[:mu_s] = jks.mu_s
    
    return eqdf
end


function find_joint_optimal_bond_measure(jf, jks,
                                         fi_sf_mu_b::Float64,
                                         fi_rf_mu_b::Float64,
                                         fi_rf_obj_val::Float64;
                                         equilibrium_type::String="all", 
                                         sf_obj_fun::Symbol=:firm_value,
                                         rf_obj_fun::Symbol=:MBR,
                                         lb::Float64=.75,
                                         ub::Float64=1.25,
                                         mu_bN::Int64=20,
                                         mu_bN2::Int64=10^5,
                                         spline_k::Int64=3,
                                         spline_bc::String="extrapolate")

    # Form Grid of mu_b candidates
    min_mu_b = .75 * minimum([fi_sf_mu_b, fi_rf_mu_b])
    max_mu_b = 1.25 * maximum([fi_sf_mu_b, fi_rf_mu_b])
    mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN)

    # Compute Optimal VB for each mu_b candidate
    jks2 = deepcopy(jks)

    # Misrepresentation case:
    eq_types = ["pooling", "separating"]
    if equilibrium_type == "all"
        println("Computing Pooling and Separating Equilibria")
    elseif equilibrium_type == "pooling"
        println("Computing Pooling Equilibrium")
        eq_types = ["pooling"]
    elseif equilibrium_type == "separating"
        println("Computing Separating Equilibrium")
        eq_types = ["separating"]
        # setfield!(jks2, :mu_s, 1.)
    else
        println("Please set equilibrium_type to pooling or separting. Exiting...")
        return
    end
    
    # dfl = @time fetch(@spawn [find_joint_optimal_vb(jf, jks2;
    #                                                 mu_b=mu_b,
    #                                                 rerun_fi_vb=true)
    #                           for mu_b in mu_b_grid])

    # # Store results in a DataFrame
    # df = vcat(dfl...)
    # sf_df = df[isnan.(df[:rf_vb]), :]
    # rf_df = df[isnan.(df[:sf_vb]), :]


    # # Interpolate Objective Functions in mu_b ################################
    # sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    # rf_objf = Dierckx.Spline1D(rf_df[:mu_b], rf_df[rf_obj_fun]; k=spline_k, bc=spline_bc)

    # # Refine mu_b_grid
    # ref_mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN2)

    optdf = DataFrame()
    for eq_type in eq_types
        # Find optimal mu_b values 
        if eq_type == "pooling"
            # cond = rf_objf(ref_mu_b_grid) .>= fi_rf_obj_val
            # filtered_mub_grid = ref_mu_b_grid[cond]

            # # Maximize Safe Firm's Objective Function
            # mub_opt = filtered_mub_grid[argmax(sf_objf(filtered_mub_grid))]
            
            
            # # Compute Results
            # eqdf = find_joint_optimal_vb(jf, jks2;
            #                              mu_b=mub_opt, rerun_fi_vb=true)
            eqdf = get_pool_eqdf(jf, jks2, mu_b_grid,
                                 fi_rf_obj_val;
                                 sf_obj_fun=sf_obj_fun,
                                 rf_obj_fun=rf_obj_fun,
                                 mu_bN2=mu_bN2)
            
        elseif eq_type == "separating" 
            # cond = rf_objf(ref_mu_b_grid) .<= fi_rf_obj_val
            # filtered_mub_grid = range(minimum(ref_mu_b_grid[cond]),
            #                           stop=maximum(ref_mu_b_grid[cond]),
            #                           length=mu_bN)


            # eqdf = get_sep_eqdf(jf, jks2, filtered_mub_grid,
            #                     fi_rf_mu_b; sf_obj_fun=sf_obj_fun)

            eqdf = get_sep_eqdf(jf, jks2, mu_b_grid,
                                fi_rf_mu_b, fi_rf_obj_val;
                                sf_obj_fun=sf_obj_fun,
                                rf_obj_fun=rf_obj_fun,
                                mu_bN=mu_bN,
                                mu_bN2=mu_bN2)

        end
        # filtered_mub_grid = ref_mu_b_grid[cond]
        
        
        # # Maximize Safe Firm's Objective Function
        # mub_opt = filtered_mub_grid[argmax(sf_objf(filtered_mub_grid))]
        

        # # Compute Results
        # eqdf = find_joint_optimal_vb(jf, jks2;
        #                              mu_b=mub_opt, rerun_fi_vb=true)
        
        # Add Objective Function and Equilibrium Type columns
        eqdf[:obj_fun] = String(sf_obj_fun)
        eqdf[2, :obj_fun] = String(rf_obj_fun)
        eqdf[:eq_type] = eq_type


        optdf = vcat(optdf, reshape_sf_rf_df(eqdf))
    end
    
    return optdf
end

