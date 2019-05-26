


function get_pool_eqdf(jf, jks, mu_b_grid,
                       fi_rf_obj_val::Float64;
                       sf_obj_fun::Symbol=:firm_value,
                       rf_obj_fun::Symbol=:MBR,
                       spline_k::Int64=3,
                       spline_bc::String="extrapolate",
                       N1::Int64=20,
                       N2::Int64=10^5)

    sf_objf, rf_objf, ref_mu_b_grid = find_joint_payoffs(jf, jks, mu_b_grid;
                                                         sf_obj_fun=sf_obj_fun,
                                                         rf_obj_fun=rf_obj_fun,
                                                         spline_k=spline_k,
                                                         spline_bc=spline_bc,
                                                         N=N2)

    # Range of mu_b values in pooling equilibrium
    cond = rf_objf(ref_mu_b_grid) .>= fi_rf_obj_val
    filtered_mu_b_grid_1, filtered_mu_b_grid_2 = find_mu_b_intervals(ref_mu_b_grid,
                                                                     cond; N=N1)

    # Maximize Safe Firm's Objective Function
    mu_b_opt = find_opt_mu_b(sf_objf,
                             filtered_mu_b_grid_1,
                             filtered_mu_b_grid_2)
    
    # Compute Results
    return find_joint_optimal_vb(jf, jks;
                                 mu_b=mu_b_opt, rerun_fi_vb=true)
end



function get_sep_eqdf(jf, jks, mu_b_grid,
                      fi_rf_mu_b::Float64,
                      fi_rf_obj_val::Float64;
                      sf_obj_fun::Symbol=:firm_value,
                      rf_obj_fun::Symbol=:MBR,
                      spline_k::Int64=3,
                      spline_bc::String="extrapolate",
                      N1::Int64=20,
                      N2::Int64=10^5)

    # Compute the Payoffs in case of Misrepresentation
    sep_misrep_jks = deepcopy(jks)
    sep_misrep_jks.mu_s = 1.
    sf_objf, rf_objf, ref_mu_b_grid = find_joint_payoffs(jf, sep_misrep_jks,
                                                         mu_b_grid;
                                                         sf_obj_fun=sf_obj_fun,
                                                         rf_obj_fun=rf_obj_fun,
                                                         spline_k=spline_k,
                                                         spline_bc=spline_bc,
                                                         N=N2)

    # Filter -> leverage values for which misrepresentation is not
    # attractive
    cond = rf_objf(ref_mu_b_grid) .<= fi_rf_obj_val
    filtered_mu_b_grid_1, filtered_mu_b_grid_2 = find_mu_b_intervals(ref_mu_b_grid,
                                                                     cond; N=N1)
    
    # Mu_b that yields maximum Safe Type's Firm Value
    # conditional on misrepresentation not being optimal for risky type.
    # Since mu_s = 1 above, the payoffs for the safe firm coincide with
    # the payoffs under full information, when debt investors can fully
    # observe firms' types.
    mu_b_opt = find_opt_mu_b(sf_objf,
                             filtered_mu_b_grid_1,
                             filtered_mu_b_grid_2)
    

    # Compute Separating Eq. Payoffs
    eqdf = separate_eq_calculator(jf, jks, mu_b_opt, fi_rf_mu_b)

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
    else
        println("Please set equilibrium_type to pooling or separting. Exiting...")
        return
    end
    
    optdf = DataFrame()
    for eq_type in eq_types
        # Find optimal mu_b values 
        if eq_type == "pooling"
            eqdf = get_pool_eqdf(jf, jks2, mu_b_grid,
                                 fi_rf_obj_val;
                                 sf_obj_fun=sf_obj_fun,
                                 rf_obj_fun=rf_obj_fun,
                                 N1=mu_bN,
                                 N2=mu_bN2)
            
        elseif eq_type == "separating" 
            eqdf = get_sep_eqdf(jf, jks2, mu_b_grid,
                                fi_rf_mu_b, fi_rf_obj_val;
                                sf_obj_fun=sf_obj_fun,
                                rf_obj_fun=rf_obj_fun,
                                N1=mu_bN,
                                N2=mu_bN2)

        end
        
        # Add Objective Function and Equilibrium Type columns
        eqdf[:obj_fun] = String(sf_obj_fun)
        eqdf[2, :obj_fun] = String(rf_obj_fun)
        eqdf[:eq_type] = eq_type


        optdf = vcat(optdf, reshape_sf_rf_df(eqdf))
    end
    
    return optdf
end

