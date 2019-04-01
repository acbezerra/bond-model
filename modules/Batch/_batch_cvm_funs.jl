

function get_k_struct(cvm, mu_b::Float64, m::Float64, c::Float64, p::Float64)
    if isnan(mu_b)
        mu_b = cvm.mu_b
    end
    if isnan(m)
        m = cvm.m
    end
    if isnan(c)
        c = cvm.c
    end
    if isnan(p)
        p = cvm.p
    end
    return mu_b, m, c, p
end


function get_p_grid(cvm; mu_b::Float64=NaN, m::Float64=NaN, 
                    c::Float64=NaN, step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)
    rdisc = get_rdisc(cvm)
    C = get_agg_c(cvm; mu_b=mu_b, c=c)
    
    # Create grid for P:
    # 1. Set Max p:
    if (m - 1/rdisc) * (1 - exp(rdisc * m)) > 0
        pMax = (C/ m) / rdisc
    else
        pMax = 3.0*cvm.pm.V0 / m
    end

    # 2. Search for minimum acceptable P value. This is the value of P for which vb = 0:
    p_grid = range(.5, stop=cvm.pm.V0, length=N1)
    vb_grid = fetch(@spawn [get_cvm_vb(cvm, cvm.pm.sigmal; 
                                       mu_b=mu_b, c=c, p=p) for p in p_grid])
    vbff = Dierckx.Spline1D(p_grid, vb_grid, k=3, bc="extrapolate")
    p_grid_ref = range(p_grid[1], stop=p_grid[end], length=N2)

    # vb must be positive
    if sum(vbff(p_grid_ref) .< 0) > 0 
        p_lower_bound = maximum(p_grid_ref[(vbff(p_grid_ref) .< 0)])
        p_grid_ref = p_grid_ref[p_grid_ref .> p_lower_bound]
    end
    p0 = p_grid_ref[argmin(abs.(vbff(p_grid_ref) .- 1e-4))]

    # 3. Return Grid (start at 1.1*p0, to avoid potential function errors)
    return range(1.001 * p0, stop = 1.1 * pMax, step = step)
end


function debt_price_diff(cvm, mu_b::Float64, c::Float64, p::Float64)
    vbl = get_cvm_vb(cvm, cvm.pm.sigmal; mu_b=mu_b, c=c, p=p)
    debt = get_cvm_debt_price(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=p)
    
    return debt - get_agg_p(cvm, mu_b=mu_b, p=p)
end


function get_p_debt_at_par(cvm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN,
                           step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)
    
    # Find P such that Debt = P

    # Given C, we search for P such that P = D, at optimal vb.
    # The function below can have more than one root, so initial
    # guess matters!
    # The solution is to the guess the lowest p value possible: 0
    # To avoid this issue, I rely on grid search.

    # Set p grid:
    pGrid = get_p_grid(cvm; mu_b=mu_b, m=m, c=c, step=step, N1=N1, N2=N2)

    # Compute squared difference between Debt and P:
    pDebtDiff = fetch(@spawn [debt_price_diff(cvm, mu_b, c, p) for p in pGrid])
    
    diff = Dierckx.Spline1D(pGrid, pDebtDiff, k=3, bc="extrapolate")
    pvec = range(pGrid[1], stop=pGrid[end], length=N2)
    
    return  pvec[argmin(abs.(diff(pvec)))]
end


function get_cvm_c_debt_at_par(bt, cvm, c::Float64; 
                               mu_b::Float64=NaN, m::Float64=NaN, 
                               step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4)
    # Capital Structure
    mu_b, m, c, _ = get_k_struct(cvm, mu_b, m, c, NaN)
    # P s.t. P = Debt
    pOpt = get_p_debt_at_par(cvm; mu_b=mu_b, m=m, c=c, step=step, N1=N1, N2=N2)

    # Store Results ##############################################
    df = DataFrame()
    
    # Parameters
    varlist = [var for var in vcat(bt.dfc.main_params, bt.dfc.fixed_params) if var != :delta]
    for var in varlist
        df[var] = bt.mi._svm_dict[var]
    end
    df[:delta] = df[:gross_delta] - df[:iota]

    df[:mu_b] = mu_b
    df[:m] = m
    df[:c] = c
    df[:p] = pOpt
    
    # Default Barrier
    vbl = get_cvm_vb(cvm, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
    df[:vb] = vbl
    
    # Debt
    df[:debt] = get_cvm_debt_price(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)  
    df = debt_at_par_diffs(cvm, df, :p)
    
    # Equity
    df[:eq_vb] = get_cvm_eq(cvm, vbl, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
    df[:eq_min_val] = df[:eq_vb][1]
    df[:equity] = get_cvm_eq(cvm, cvm.pm.V0, cvm.pm.sigmal; mu_b=mu_b, c=c, p=pOpt)
                
    # Firm Value, Leverage & ROE
    df = non_interp_values(cvm, df)

    # Rearrange Columns
    return df[bt.dfc.dfcols]
end


function get_cvm_debt_at_par(bt, cvm;
                             mu_b::Float64=NaN, m::Float64=NaN, 
                             step::Float64=.2, N1::Int64=10^2, N2::Int64=10^4,
                             save_soldf::Bool=true,
                             soldf_name::String=soldf_name)

    resdfs = fetch(@spawn [get_cvm_c_debt_at_par(bt, cvm, c;
                                                 mu_b=mu_b, m=m, 
                                                 step=step, N1=N1, N2=N2)
                           for c in bt.coupon_grid])

    soldf = vcat(resdfs...)

    if save_soldf
        println("Saving CVM results...")
        CSV.write(string(bt.mi.comb_res_path, "/", soldf_name, ".csv"), soldf)
    end

    return soldf
end


function optimal_cvm_capital_struct(bt, cvm; df::DataFrame=DataFrame(),
                                    dfname::String=Batch.soldf_name,
                                    interp_polyOrder::Int64=3,
                                    filter_windowSize::Int64=5 * 10^2 + 1, 
                                    filter_polyOrder::Int64=3,
                                    save_results::Bool=true,
                                    opt_k_struct_df_name::String=opt_k_struct_df_name)
    
    # Load Solutions
    if isempty(df)
        df = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"),
                      types=fill(Float64, 28))
    end
    
    # Interpolate and Filter Optimal Capital Structure
    sgF = filter_k_struct(df; interp_polyOrder=interp_polyOrder,
                          filter_windowSize=filter_windowSize,
                          filter_polyOrder=filter_polyOrder,
                          interp_vars=[:p, :vb, :debt, :equity])
    
    # Compute Optimal Capital Structure
    cpos = minimum(findlocalmaxima(sgF[:firm_value]))
    
    optDF = DataFrame()

    # Parameters ######################
    for var in vcat(bt.dfc.main_params, bt.dfc.fixed_params, [:mu_b, :m])
        optDF[var] = df[var][1]
    end
    # #################################
    
    #  Capital Struct #################
    optDF[:c] = sgF[:cgrid][cpos]
    optDF[:p] = sgF[:p][cpos]
    optDF[:vb] = sgF[:vb][cpos]
    # #################################
    
    # Debt ###########################
    optDF[:sg_debt] = sgF[:debt][cpos]
    optDF[:debt] = get_cvm_debt_price(cvm, optDF[:vb][1], 
                                     optDF[:sigmal][1]; 
                                     mu_b=optDF[:mu_b][1], 
                                     c=optDF[:c][1], 
                                     p=optDF[:p][1])
    optDF = debt_at_par_diffs(cvm, optDF, :p)
    # #################################
    
    # Equity #########################
    optDF[:eq_vb] = get_cvm_eq(cvm, 
                              optDF[:vb][1], 
                              optDF[:sigmal][1]; 
                              mu_b=optDF[:mu_b][1], 
                              c=optDF[:c][1], 
                              p=optDF[:p][1])
    optDF[:eq_min_val] = optDF[:eq_vb][1]
    optDF[:sg_equity] = sgF[:equity][cpos]
    optDF[:equity] = get_cvm_eq(cvm, 
                               optDF[:V0][1], 
                               optDF[:sigmal][1]; 
                               mu_b=optDF[:mu_b][1], 
                               c=optDF[:c][1], 
                               p=optDF[:p][1])
    # #################################
    
    optDF = non_interp_values(cvm, optDF)
    
    
    # CVM NaN columns
    optDF[:cvml_vb] = NaN
    optDF[:cvmh_vb] = NaN
    optDF[:eq_deriv] = NaN
    optDF[:eq_min_val] = NaN
     
    
    # Reoder columns
    cols = opt_k_struct_cols(bt)

    if save_results
        CSV.write(string(bt.mi.comb_res_path, "/", opt_k_struct_df_name, ".csv"), optDF)
    end

    return optDF[cols]
end


function cvm_results_loader(comb_num::Int64; 
                            optdf_name::String=opt_k_struct_df_name, 
                            coltypes::Array{DataType,1}=vcat(fill(Float64, 27), 
                                                             [Bool], 
                                                             fill(Float64, 6)),
                            display_msgs::Bool=false)
    
    bt = get_bt(; model="cvm", comb_num=comb_num, display_msgs=display_msgs)
    
    id_cols = [:comb_num, :m, :m_comb_num]
    opt_cols = opt_k_struct_cols(bt)
    
    
    optDF = CSV.read(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), types=coltypes)
    optDF = hcat(get_batch_comb_num(bt; display_msgs=display_msgs)[[:comb_num, :m_comb_num]], optDF)

    cols_order = vcat(id_cols, [x for x in opt_cols if !(x in id_cols)])
    return optDF[cols_order]
end


function compile_cvm_opt_results(; optdf_name::String=opt_k_struct_df_name, 
                                  coltypes::Array{DataType,1}=vcat(fill(Float64, 27), 
                                                                   [Bool], 
                                                                   fill(Float64, 6)),
                                  display_msgs::Bool=false,
                                  save_results::Bool=true)
    
    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)
    
    optDFs = @time fetch(@spawn [cvm_results_loader(cn;
                                                    optdf_name=optdf_name,
                                                    coltypes=coltypes, 
                                                    display_msgs=display_msgs) 
                                 for cn in bt.bp.df[:comb_num]])
    cvmOptDF= sort!(vcat(optDFs...), [:m, :m_comb_num])
    
    if save_results
        CSV.write(string(bt.mi.batch_res_path, "/", optdf_name, ".csv"), cvmOptDF)
    end
    
    return cvmOptDF
end


function load_cvm_opt_results_df(; optdf_name::String=opt_k_struct_df_name, 
                                 coltypes::Array{DataType,1}=cvm_opt_k_struct_df_coltypes,
                                 display_msgs::Bool=false)
    
    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)
    println("Loading optimal results dataframe...")
    return CSV.read(string(bt.mi.batch_res_path, "/", opt_k_struct_df_name, ".csv"); types=coltypes)
end
