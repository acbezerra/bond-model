function load_svm_opt_results(bt;
                          comb_num::Integer=0,
                          m::Float64=NaN,
                          m_comb_num::Integer=0,
                          display_msgs::Bool=true,
                          toldf::DataFrame=Batch.toldf,
                          drop_fail::Bool=false,
                          save_soldf::Bool=true,
                          soldf_name::String=soldf_name,
                          interp_polyOrder::Integer=3,
                          filter_windowSize::Integer=5 * 10^2 + 1, 
                          filter_polyOrder::Integer=3,
                          save_optdf::Bool=true,
                          optdf_name::String=optdf_name)
    
    bt = get_bt(;comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                display_msgs=display_msgs)

    opt_k_struct_found = false
    if string(optdf_name, ".csv") in readdir(bt.mi.comb_res_path)
        cols = opt_k_struct_cols(bt)
        coltypes = map(x -> begin
                              if x == :eq_negative
                                return Bool
                              end
                              return Float64
                            end, cols)
        optDF = CSV.read(string(bt.mi.comb_res_path, "/", optdf_name, ".csv"), types=coltypes)

        if !isnan(optDF[:c][1])
            opt_k_struct_found=true
        end
    end
    

    if opt_k_struct_found
        optDF = hcat(get_batch_comb_num(bt; display_msgs=display_msgs)[[:comb_num, :m_comb_num]], optDF)
        id_cols = [:comb_num, :m, :m_comb_num]
        cols_order = vcat(id_cols, [x for x in names(optDF) if !(x in id_cols)])
        return optDF[cols_order]
    else
        optDF = get_opt_results(bt;
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
    end
    return optDF
end


function load_svm_opt_results_df(bt;
                                 m::Float64=NaN,
                                 opt_k_struct_df_name::String=opt_k_struct_df_name,
                                 coltypes::Array{DataType,1}=opt_k_struct_df_coltypes,
                                 recompute::Bool=false)
    if !isnan(m)
        bt = get_bt(;bt=bt, m=m, m_comb_num=1)
    end
    
    if recompute || !(string(opt_k_struct_df_name, ".csv") in readdir(bt.mi.maturity_path))
        println("Compiling results...")
        return compile_svm_opt_results(bt; m=m, opt_k_struct_df_name=opt_k_struct_df_name)
        
    else
        if isnan(m)
            fpath = bt.mi.batch_res_path
        else
            fpath = bt.mi.maturity_path
        end
        println("Loading optimal results dataframe...")
        return CSV.read(string(fpath, "/", opt_k_struct_df_name, ".csv"); types=coltypes)
    end
end
