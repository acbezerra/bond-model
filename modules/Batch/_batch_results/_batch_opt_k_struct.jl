function opt_k_struct_cols(bt; de_cols::Array{Symbol,1}=[:debt_diff, :eq_deriv, :eq_min_val])
    share_cols = vcat([:sg_debt, :debt, :sg_equity], 
                      [x for x in bt.dfc.share_values if x != :debt])
    debt_diff_cols = [x for x in bt.dfc.debt_vars if !(x in de_cols)]
    eqdf_cols = [x for x in bt.dfc.equity_vars if !(x in de_cols)]
    return vcat(:obj_fun, bt.dfc.main_params, bt.dfc.k_struct_params,
                [:cvml_vb], [:cvmh_vb],
                de_cols, share_cols, bt.dfc.fixed_params, 
                debt_diff_cols, eqdf_cols)
end


function optimal_capital_struct(bt, svm;
                                firm_obj_fun::Symbol=:MBR,
                                df::DataFrame=DataFrame(),
                                dfname::String=soldf_name,
                                interp_polyOrder::Integer=3,
                                filter_windowSize::Integer=5 * 10^2 + 1, 
                                filter_polyOrder::Integer=3)
    
    # Load Solutions
    if isempty(df)
        df = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"),
                      types=vcat(fill(Float64, 23), [Bool], fill(Float64, 5)))
    end

    # Discard failed results
    df = df[df[:p_filter_success].==true, :]
    
    # Interpolate and Filter Optimal Capital Structure
    sgF = filter_k_struct(df; interp_polyOrder=interp_polyOrder,
                          filter_windowSize=filter_windowSize,
                          filter_polyOrder=filter_polyOrder)
    
    # Compute Optimal Capital Structure
    sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]
    sgF[:MBR] = get_mbr(get_param(svm, :V0), sgF[:debt], sgF[:equity])
    
    cpos = minimum(findlocalmaxima(sgF[firm_obj_fun]))
    
    # Run Equity Finite Differences Method
    eqDF = eq_fd(svm;
                 vbl=sgF[:opt_vb][cpos],
                 mu_b=svm.mu_b, 
                 c=sgF[:cgrid][cpos], 
                 p=sgF[:p][cpos])
    
    # Add Parameter Values
    eqDF[:cvml_vb] = sgF[:cvml_vb][cpos]
    eqDF[:cvmh_vb] = sgF[:cvmh_vb][cpos]
    eqDF[:sg_debt] = sgF[:debt][cpos]
    eqDF[:sg_equity] = sgF[:equity][cpos]
    
    # Add Debt Functions
    eqDF = debt_at_par_diffs(svm, eqDF, :p)

    # Add Column with firm_obj_fun
    eqDF[:obj_fun] = String(firm_obj_fun)
    
    # Reoder columns
    cols = opt_k_struct_cols(bt)

    return eqDF[cols]
end

function save_combination_opt_k_struct(bt, eqDF::DataFrame;
                                       idcols::Array{Symbol, 1}=vcat([:m, :obj_fun],
                                                                     bt.dfc.main_params,
                                                                     bt.dfc.fixed_params),
                                       dfname::String=optdf_name,
                                       coltypes::Array{DataType,1}=opt_k_struct_df_coltypes)

    if (string(dfname, ".csv") in readdir(bt.mi.comb_res_path))
        try
            seqDF = CSV.read(string(bt.mi.comb_res_path, "/", dfname, ".csv"),
                             types=coltypes)

            cond = .&(vcat([seqDF[:obj_fun] .== eqDF[:obj_fun]], 
                           [(abs.(seqDF[col] .- eqDF[col]) .< 1e-6) for col in idcols if col != :obj_fun])...)
            if sum(cond) == 0
                println("No match found. Appending row...")
                seqDF = vcat(seqDF, eqDF)
            elseif sum(cond) == 1
                println("Match found! Replacing row...")
                seqDF[cond, :] = eqDF
            else 
                println("Multiple matches found. Refine ID columns. Returning...")
                return
            end 
        catch
            println("Unable to load DataFrame.")
            seqDF = eqDF
        end
    else
        seqDF = eqDF
    end
        
    # Save DataFrame 
    CSV.write(string(bt.mi.comb_res_path, "/", dfname, ".csv"), seqDF)
end


function compile_optimal_cap_struct(bt, svm;
                                    firm_obj_fun::Symbol=:MBR,
                                    toldf::DataFrame=toldf,
                                    use_all_eqdf::Bool=true,
                                    drop_fail::Bool=false,
                                    save_soldf::Bool=true,
                                    soldf_name::String=soldf_name,
                                    interp_polyOrder::Integer=3,
                                    filter_windowSize::Integer=5 * 10^2 + 1, 
                                    filter_polyOrder::Integer=3,
                                    replace_optdf::Bool=false,
                                    save_optdf::Bool=true,
                                    optdf_name::String=optdf_name)

    
    # Collect and Process Results
    soldf = process_combination_results(bt, svm;
                                        toldf=toldf,
                                        use_all_eqdf=use_all_eqdf,
                                        drop_fail=drop_fail,
                                        save_df=save_soldf,
                                        dfname=soldf_name)

    # Find Optimal Capital Structure
    eqDF = optimal_capital_struct(bt, svm;
                                  firm_obj_fun=firm_obj_fun,
                                  df=soldf,
                                  dfname=soldf_name,
                                  interp_polyOrder=interp_polyOrder,
                                  filter_windowSize=filter_windowSize,
                                  filter_polyOrder=filter_polyOrder)

    # Add combination IDs
    combdf = get_batch_comb_num(bt)
    eqDF[:comb_num] =combdf[1, :comb_num]
    eqDF[:m_comb_num] =combdf[1, :comb_num]       
    id_cols = [:comb_num, :m, :m_comb_num]
    cols_order = vcat(id_cols, [x for x in names(eqDF) if !(x in id_cols)])
    eqDF = eqDF[cols_order]

    if replace_optdf
        CSV.write(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), eqDF)
    elseif save_optdf
        save_combination_opt_k_struct(bt, eqDF; dfname=optdf_name)
    end
    
    return eqDF
end


function get_opt_results(bt;
                         comb_num::Integer=0,
                         firm_obj_fun::Symbol=:MBR,
                         use_all_eqdf::Bool=true,
                         m::Float64=0.,
                         m_comb_num::Integer=0,
                         display_msgs::Bool=true,
                         toldf::DataFrame=Batch.toldf,
                         drop_fail::Bool=false,
                         save_soldf::Bool=true,
                         soldf_name::String=soldf_name,
                         interp_polyOrder::Integer=3,
                         filter_windowSize::Integer=5 * 10^2 + 1, 
                         filter_polyOrder::Integer=3,
                         replace_optdf::Bool=false,
                         save_optdf::Bool=true,
                         optdf_name::String=optdf_name)
    
    bt, svm = get_bt_svm(;comb_num=comb_num,
                         m=m, m_comb_num=m_comb_num,
                         display_msgs=display_msgs)

    try
        return compile_optimal_cap_struct(bt, svm;
                                          firm_obj_fun=firm_obj_fun,
                                          toldf=toldf,
                                          use_all_eqdf=use_all_eqdf,
                                          drop_fail=drop_fail,
                                          save_soldf=save_soldf,
                                          soldf_name=soldf_name,
                                          interp_polyOrder=interp_polyOrder,
                                          filter_windowSize=filter_windowSize, 
                                          filter_polyOrder=filter_polyOrder,
                                          replace_optdf=replace_optdf,
                                          save_optdf=save_optdf,
                                          optdf_name=optdf_name)
        
        
        # return eqDF[vcat([:comb_num], [x for x in names(eqDF) if x !=:comb_num])]
    catch
        cols = opt_k_struct_cols(bt)
        eqDF = DataFrame()
        for col in cols
            if string(col) in keys(bt.mi._svm_dict)
               eqDF[col] = bt.mi._svm_dict[string(col)]
            elseif col != :eq_negative
               eqDF[col] = NaN
            else
               eqDF[col] = true
            end
        end
        # eqDF = hcat(get_batch_comb_num(bt)[[:comb_num, :m_comb_num]], eqDF)
        # id_cols = [:comb_num, :m, :m_comb_num]
        # cols_order = vcat(id_cols, [x for x in names(eqDF) if !(x in id_cols)])
        # return eqDF[cols_order]
        # Add combination IDs
        eqDF[:comb_num] = comb_num
        combdf = get_batch_comb_num(bt)
        optDF[:comb_num] =combdf[1, :comb_num]
        optDF[:m_comb_num] =combdf[1, :comb_num]       
        id_cols = [:comb_num, :m, :m_comb_num]
        cols_order = vcat(id_cols, [x for x in names(eqDF) if !(x in id_cols)])
        return eqDF[cols_order]
    end
end


function compile_svm_opt_results(bt; m::Float64=NaN,
                                 firm_obj_fun::Symbol=:MBR,
                                 use_all_eqdf::Bool=true,
                                 display_msgs::Bool=false,
                                 toldf::DataFrame=Batch.toldf,
                                 drop_fail::Bool=false,
                                 save_soldf::Bool=true,
                                 soldf_name::String=soldf_name,
                                 interp_polyOrder::Integer=3,
                                 filter_windowSize::Integer=5 * 10^2 + 1, 
                                 filter_polyOrder::Integer=3,
                                 save_optdf::Bool=true,
                                 optdf_name::String=optdf_name,
                                 recompute_comb_opt_res::Bool=true,
                                 replace_optdf::Bool=false,
                                 save_results::Bool=true,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name)

    if !isnan(m)
        comb_nums = bt.bp.df[bt.bp.df[:m] .== m, :comb_num]
    else
        comb_nums = bt.bp.df[:comb_num]
    end
    #    comb_nums = range(1, stop=size(hcat(bt._params_combs...)', 1))
    if recompute_comb_opt_res
        optDF_LL = @time fetch(@spawn [get_opt_results(bt;
                                                       comb_num=comb_num,
                                                       firm_obj_fun=firm_obj_fun,
                                                       use_all_eqdf=use_all_eqdf,
                                                       display_msgs=display_msgs,
                                                       toldf=toldf,
                                                       drop_fail=drop_fail,
                                                       save_soldf=save_soldf,
                                                       soldf_name=soldf_name,
                                                       interp_polyOrder=interp_polyOrder,
                                                       filter_windowSize=filter_windowSize, 
                                                       filter_polyOrder=filter_polyOrder,
                                                       replace_optdf=replace_optdf,
                                                       save_optdf=save_optdf,
                                                       optdf_name=optdf_name)
                                   for comb_num in comb_nums])
    else
        optDF_LL = @time fetch(@spawn [load_svm_opt_results(bt;
                                                            comb_num=comb_num,
                                                            display_msgs=display_msgs,
                                                            toldf=toldf,
                                                            drop_fail=drop_fail,
                                                            save_soldf=save_soldf,
                                                            soldf_name=soldf_name,
                                                            interp_polyOrder=interp_polyOrder,
                                                            filter_windowSize=filter_windowSize, 
                                                            filter_polyOrder=filter_polyOrder,
                                                            save_optdf=save_optdf,
                                                            optdf_name=optdf_name)
                                       for comb_num in comb_nums])
    end
    optDF = vcat(optDF_LL...)

    # Filter
    optDF = optDF[optDF[:obj_fun] .== String(firm_obj_fun), :]

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
        CSV.write(string(fpath, "/", dfname, ".csv"), optDF)
    end
    
    return optDF
end

