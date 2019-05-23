

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
                
    # Firm Value, Leverage & MBR
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


function optimal_cvm_capital_struct(bt, cvm;
                                    firm_obj_fun::Symbol=:firm_value,
                                    df::DataFrame=DataFrame(),
                                    dfname::String=Batch.soldf_name,
                                    interp_polyOrder::Int64=3,
                                    filter_windowSize::Int64=5 * 10^2 + 1, 
                                    filter_polyOrder::Int64=3,
                                    replace_optdf::Bool=false,
                                    save_optdf::Bool=true,
                                    optdf_name::String=optdf_name)
    
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
    sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]
    sgF[:MBR] = get_mbr(get_param(cvm, :V0), sgF[:debt], sgF[:equity])
    
    # Compute Optimal Capital Structure
    cpos = minimum(findlocalmaxima(sgF[firm_obj_fun]))
    
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
    optDF[:eq_deriv_min_val] = NaN

    optDF[:eq_negative] = false
     
    # Add Column with firm_obj_fun
    optDF[:obj_fun] = String(firm_obj_fun)

    # Add combination IDs
    combdf = get_batch_comb_num(bt)
    optDF[:comb_num] =combdf[1, :comb_num]
    optDF[:m_comb_num] =combdf[1, :comb_num]       
    id_cols = [:comb_num, :m, :m_comb_num]
    cols_order = vcat(id_cols, [x for x in opt_k_struct_cols(bt) if !(x in id_cols)])
    optDF = optDF[cols_order]

    if replace_optdf
        CSV.write(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), optDF)
    elseif save_optdf
        idcols = [x for x in vcat([:m, :obj_fun],
                                  bt.dfc.main_params,
                                  bt.dfc.fixed_params) if !(x in [:lambda, :sigmah])]
        save_combination_opt_k_struct(bt, optDF; dfname=optdf_name,
                                      idcols=idcols,
                                      coltypes=cvm_opt_k_struct_df_coltypes)
    end
    
    return optDF
end


function cvm_results_loader(comb_num::Int64; 
                            optdf_name::String=optdf_name, 
                            coltypes::Array{DataType,1}=cvm_opt_k_struct_df_coltypes,
                            display_msgs::Bool=false)
    
    bt = get_bt(; model="cvm", comb_num=comb_num, display_msgs=display_msgs)
    optDF = CSV.read(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), types=coltypes)
    
    return optDF
end


function get_cvm_opt_results(; comb_num::Integer=0,
                             firm_obj_fun::Symbol=:firm_value,
                             m::Float64=0.,
                             m_comb_num::Integer=0,
                             display_msgs::Bool=true,
                             dfname::String=soldf_name,
                             replace_optdf::Bool=false,
                             save_optdf::Bool=true,
                             optdf_name::String=optdf_name)
    
    bt, cvm = get_bt_cvm(;comb_num=comb_num,
                         m=m, m_comb_num=m_comb_num,
                         display_msgs=display_msgs)

    
    return optimal_cvm_capital_struct(bt, cvm;
                                      firm_obj_fun=firm_obj_fun,
                                      dfname=dfname,
                                      replace_optdf=replace_optdf,
                                      save_optdf=save_optdf,
                                      optdf_name=optdf_name)
end


function compile_cvm_opt_results(; m::Float64=NaN,
                                 firm_obj_fun::Symbol=:firm_value,
                                 display_msgs::Bool=false,
                                 dfname::String=soldf_name,
                                 replace_optdf::Bool=false,
                                 save_optdf::Bool=true,
                                 optdf_name::String=optdf_name,
                                 save_results::Bool=true,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)

    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)

    if !isnan(m)
        comb_nums = bt.bp.df[abs.(bt.bp.df[:m] .- m) .< 1e-6, :comb_num]
    else
        comb_nums = bt.bp.df[:comb_num]
    end

    optDFs = @time fetch(@spawn [get_cvm_opt_results(; comb_num=comb,
                                                     firm_obj_fun=firm_obj_fun,
                                                     display_msgs=display_msgs,
                                                     dfname=soldf_name,
                                                     replace_optdf=replace_optdf,
                                                     save_optdf=save_optdf,
                                                     optdf_name=optdf_name)
                                 for comb in comb_nums])

    cvmOptDF = sort!(vcat(optDFs...), [:m, :m_comb_num])
   
    if save_results
        println("Saving compiled results...")
        extension = lowercase(string(firm_obj_fun))
        dfname = string(opt_k_struct_df_name, "_", extension)

        if isnan(m)
            bt = set_par_dict(bt; comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.batch_res_path
        else
            bt = set_par_dict(bt; m=m, m_comb_num=1, display_msgs=display_msgs)
            bt = set_comb_res_paths(bt)
            fpath = bt.mi.maturity_path
        end
        CSV.write(string(fpath, "/", dfname, ".csv"), cvmOptDF)
    end
    
    return cvmOptDF
end


function load_cvm_opt_results_df(; m::Float64=NaN,
                                 optdf_name::String=opt_k_struct_df_name,
                                 firm_obj_fun::Symbol=:firm_value,
                                 coltypes::Array{DataType,1}=cvm_opt_k_struct_df_coltypes,
                                 display_msgs::Bool=false,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)
    
    bt = get_bt(; model="cvm", comb_num=1, display_msgs=display_msgs)

    extension = lowercase(string(firm_obj_fun))
    dfname = string(opt_k_struct_df_name, "_", extension)
    if isnan(m)
        bt = set_par_dict(bt; comb_num=1, display_msgs=display_msgs)
        bt = set_comb_res_paths(bt)
        fpath = bt.mi.batch_res_path
    else
        bt = set_par_dict(bt; m=m, m_comb_num=1, display_msgs=display_msgs)
        bt = set_comb_res_paths(bt)
        fpath = bt.mi.maturity_path
    end

    try 
        println("Loading optimal results dataframe...")
        return CSV.read(string(fpath, "/", dfname, ".csv"), types=cvm_opt_k_struct_df_coltypes)
    catch
        println("Unable to load optimal results dataframe. Recomputing...")
        return compile_cvm_opt_results(; m=m,
                                       firm_obj_fun=firm_obj_fun,
                                       opt_k_struct_df_name=opt_k_struct_df_name)
    end
end




    
#     println("Loading optimal results dataframe...")
#     optDF = CSV.read(string(bt.mi.batch_res_path, "/", opt_k_struct_df_name, ".csv"); types=coltypes)

#     if !(:MBR in names(optDF))
#         optDF[:MBR] = get_mbr(unique(optDF[:V0])[1], optDF[:debt], optDF[:equity])
#     end
#     return optDF
# end
