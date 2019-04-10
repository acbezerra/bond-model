

function find_joint_optimal_bond_measure(jf, jks,
                                         fi_sf_mu_b::Float64,
                                         fi_rf_mu_b::Float64,
                                         fi_rf_obj_val::Float64;
                                         equilibrium_type::String="pooling", 
                                         sf_obj_fun::Symbol=:firm_value,
                                         rf_obj_fun::Symbol=:ROE,
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
    if equilibrium_type == "pooling"
        println("Computing Pooling Equilibrium")
    elseif equilibrium_type == "separating"
        println("Computing Separating Equilibrium")
        setfield!(jks2, :mu_s, 1.)
    else
        println("Please set equilibrium_type to pooling or separting. Exiting...")
        return
    end
    
    dfl = @time fetch(@spawn [find_joint_optimal_vb(jf, jks2;
                                                    mu_b=mu_b,
                                                    rerun_fi_vb=true)
                              for mu_b in mu_b_grid])

    # Store results in a DataFrame
    df = vcat(dfl...)
    sf_df = df[isnan.(df[:rf_vb]), :]
    rf_df = df[isnan.(df[:sf_vb]), :]


    # Interpolate Objective Functions in mu_b ################################
    sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun], k=spline_k, bc=spline_bc)
    rf_objf = Dierckx.Spline1D(rf_df[:mu_b], rf_df[rf_obj_fun]; k=spline_k, bc=spline_bc)

    # Refine mu_b_grid
    ref_mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN2)

    # Find optimal mu_b values 
    if equilibrium_type == "pooling"
        cond = rf_objf(ref_mu_b_grid) .>= fi_rf_obj_val
    elseif equilibrium_type == "separating" # Misrepresentation
        cond = rf_objf(ref_mu_b_grid) .<= fi_rf_obj_val
    end
    filtered_mub_grid = ref_mu_b_grid[cond]
        

    # Maximize Safe Firm's Objective Function
    mub_opt = filtered_mub_grid[argmax(sf_objf(filtered_mub_grid))]


    # Compute Results
    optdf = find_joint_optimal_vb(jf, jks2;
                                  mu_b=mub_opt, rerun_fi_vb=true)

    # Add Objective Function columns
    optdf[:obj_fun] = String(sf_obj_fun)
    optdf[2, :obj_fun] = String(rf_obj_fun)

    return optdf
end

