
function deriv(f, x, h)
    return (f[x+h] - f[x])/h
end


function eq_fin_diff_core(svm, vbl, v_grid, v_sub_grid, bond_prices)
        V_sub_grid = vbl * exp.(v_sub_grid)
       
        # #################################
        # ######### Equity Values #########
        # #################################
        # Upper Barrier: Value of Equity
        eq_max = get_cvm_eq(vbl * exp(v_grid[1]), svm.pm.sigmah, svm)

        # Note: I can use the CVM value here because both 
        # the cvm and svm equity values will have converged 
        # to the credit-risk-free equity value in the upper
        # barrier.
        
        # ##################################
        # ##### Baseline Model Equity Values #####
        # ##################################
        println("Computing Constant Volatility Equity Values")
        bsm_eqh_sub_vals_Future = @spawn [get_bsm_eq_vals(cvm, V) for V= V_sub_grid]
        bsm_eqh_sub_vals = fetch(bsm_eqh_sub_vals_Future)
    
        # Interpolate
        itp = interpolate(reverse(bsm_eqh_sub_vals), BSpline(Cubic(Line())), OnGrid())
        bsm_interp_eqh = Interpolations.scale(itp, reverse(v_sub_grid))  # Scale

        # Equity Values
        bsm_eqh = bsm_interp_eqh(v_grid[2:end-1])
        println("Finished computing Constant Volatility Equity Values")
                    
        # #################################
        # ######### Coefficients: #########
        # #################################
        deltav = v_grid[1] - v_grid[2] # or v_grid[-2] # v_grid[0] / np.float(len(v_grid))
        nu = get_rgrow(svm) - .5 * get_param(svm, "sigmal")^2

        qu = .5 * (nu / deltav + get_param(svm, "sigmal") ^ 2 / (deltav ^ 2))
        qd = .5 * (-nu / deltav + get_param(svm, "sigmal") ^ 2 / (deltav ^ 2))
        qm = -(get_param(svm, "r") + get_param(svm, "lambda") + 
                                    get_param(svm, "sigmal") ^ 2 / (deltav ^ 2))

        # Present Value of Debt:
        pv_debt = get_pv_rfdebt(svm)
                    
        # Gamma Vector:
        Gamma = get_param(svm, "delta") * vbl * exp(reverse(v_grid)) - \
                (1 - get_param(svm, "pi")) * get_param(svm, "C") + \
                bond_prices - get_param(svm, "p") + \
                get_param(svm, "lambda") * bsm_eqh

        println("Shape of Gamma matrix: ", string(Gamma.shape))
        Gamma[1] += qu * eq_max

        eq_vbl = max(0., get_param(svm, "alpha") * vbl - pv_debt)
        println("eq_vbl: ", string(eq_vbl))
        Gamma[end] += qd * eq_vbl

        # A Matrix:
        A = (qm * Array(Diagonal(ones(length(Gamma)))) +
             qd * Array(LowerTriangular(ones(length(Gamma), length(Gamma))) - 
                            Diagonal(ones(length(Gamma)))) +
             qu * Array(UpperTriangular(ones(length(Gamma), length(Gamma))) - 
                            Diagonal(ones(length(Gamma)))))

    
        # ###### Compute Pre-Volatility Shock Equity Function: ######
        # Form Function and add lim_v->infty E(v)
        eq_vals = vcat(eq_max, - inv(A) * Gamma, eq_vbl)

        # Interpolate to back-out equity value at VB:
        eq_interp_itp = interpolate(reverse(eq_vals), BSpline(Cubic(Line())), OnGrid())
        eq_pchip = Interpolations.scale(itp, vbl * exp(reverse(v_grid)))  # Scale

        # Compute Derivative at Default Barrier:
        h = 1e-5
        eq_pchip_deriv = [deriv(eq_pchip, V, h) for V=vbl * exp(reverse(v_grid))]
        eq_deriv = eq_pchip_deriv[1]

        # Equity Values

        # Equity set as function of V:
###### CHECK ! e0 = float(eq_pchip(get_param(svm, "V0")))

        # Look for negative values
        eq_min_val = minimum(reverse(eq_vals)[2:end])  # interpolate in a neighborhood
        eq_negative = eq_min_val .< -0.1  # & (np.abs(np.min(eq_vals)) > .025)
        eq_deriv_min_val = minimum(eq_pchip_deriv)

                                                        
        # Equity set as function of V:
        e0 = float(eq_pchip[get_V0(svm)])


        # Look for negative values
        eq_min_val = minimum(reverse(eq_vals)[2:end])  # interpolate in a neighborhood
        eq_negative = eq_min_val .< -0.1  # & (np.abs(np.min(eq_vals)) > .025)
        eq_deriv_min_val = minimum(eq_pchip_deriv) 

        eq_dict = Dict("e0" => e0,
                       "eq_max" => eq_max,
                       "eq_vals"=>  eq_vals,
                       "eq_pchip" =>  eq_pchip,
                       "eq_deriv" => eq_deriv,
                       "eq_pchip_deriv" => eq_pchip_deriv,
                       "eq_vb" => eq_vbl,
                       "eq_min_val" => eq_min_val,
                       "eq_negative" => eq_negative,
                       "eq_deriv_min_val" => eq_deriv_min_val)

        println("Total Equity Core Function Computation Time: ", 
                string(tic() - tic_bsm_eq))
                                                                
        return Gamma, A, eq_dict
end



