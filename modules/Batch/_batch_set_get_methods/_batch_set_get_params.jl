
# ############## Set Parameter Combinations ##############
function set_param_comb_df(bt)
    # Create Parameter Combinations DataFrame
    pcdf = DataFrame(hcat(bt.bp._params_combs...)')
    cols_dict = Dict([names(pcdf)[i] => Symbol(bt.bp._params_order[i]) 
                      for i in 1:size(names(pcdf), 1)])
    rename!(pcdf, cols_dict)
    pcdf[:comb_num] = 1:size(pcdf, 1)

    # Important: sort by :m
    sort!(pcdf, :m)
    # Count rows by :m
    pcdf[:m_comb_num] = by(pcdf, :m, m_comb_num = :m => x -> 1:size(x, 1))[:m_comb_num]

    idcols = [:comb_num, :m, :m_comb_num]
    bt.bp.df = pcdf[vcat(idcols, [x for x in names(pcdf) if !(x in idcols)])]

    return bt
end


function set_params_combs(bt)
    bt.bp._params_combs = [Array([x1, x2, x3, x4, x5, x6, x7, x8, x9]) 
                       for x1 in bt.bp._param_values_dict[bt.bp._params_order[1]], 
                           x2 in bt.bp._param_values_dict[bt.bp._params_order[2]],
                           x3 in bt.bp._param_values_dict[bt.bp._params_order[3]],
                           x4 in bt.bp._param_values_dict[bt.bp._params_order[4]],
                           x5 in bt.bp._param_values_dict[bt.bp._params_order[5]],
                           x6 in bt.bp._param_values_dict[bt.bp._params_order[6]],
                           x7 in bt.bp._param_values_dict[bt.bp._params_order[7]],
                           x8 in bt.bp._param_values_dict[bt.bp._params_order[8]], 
                           x9 in bt.bp._param_values_dict[bt.bp._params_order[9]]]

    return set_param_comb_df(bt)
end
# ########################################################


# ################### Store Parameters ###################
# function set_par_dict(bt, comb_num)
#     bt.mi._svm_dict = Dict()
#     for par in bt._params_order
#         bt.mi._svm_dict[par] = bt._params_combs[comb_num][findall(x-> x==par, bt._params_order)[1]]
#     end
    
#     bt.mi._svm_dict = merge(bt._common_params, bt.mi._svm_dict)
    
#     return bt
# end


function set_par_dict(bt; comb_num::Integer=0,
                      m::Float64=NaN, m_comb_num::Integer=0,
                      display_msgs=true)
    # No combination entered
    cond1 = .&(comb_num == 0, (isnan(m) || m_comb_num == 0))
    # Two combinations entered
    cond2 = .&(comb_num > 0, !isnan(m),
               m_comb_num > 0)
    if cond1
        println("Please enter a combination number. Exiting...")
        return bt
    elseif cond2
        println("Two identifiers entered: (i) unique combination ID, 
                        (ii) (m, m_comb_num) ID pair.")
        println("Function will use the unique combination ID.")
    end
    
    bt.mi._svm_dict = Dict()
    if comb_num > 0
        if display_msgs
            println("Setting parameter dictionary using unique combination ID...")
        end
        for par in bt.bp._params_order
            bt.mi._svm_dict[par] = bt.bp._params_combs[comb_num][
                    findall(x-> x==par, bt.bp._params_order)[1]]
        end
    elseif .&(m > 0., m_comb_num > 0)
        if display_msgs
            println("Setting parameter dictionary using (m, m_comb_num) ID pair...")
        end
        # Slice Parameter Combinations DataFrame
        cond = .&(abs.(bt.bp.df[:m] .- m) .< 1e-6, 
                  abs.(bt.bp.df[:m_comb_num] .- m_comb_num) .< 1e-6)
        for par in bt.bp._params_order
            bt.mi._svm_dict[par] = bt.bp.df[cond, par][1]
        end
    end

    bt.mi._svm_dict = merge(bt.bp._common_params, bt.mi._svm_dict)
    
    return bt
end


function get_par_dict(bt)
    return bt.mi._svm_dict
end


# ########################################################


function get_batch_comb_num(bt; display_msgs::Bool=true, 
                            svm_dict::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
                            iota::Float64=NaN,
                            kappa::Float64=NaN,
                            lambda::Float64=NaN,
                            sigmah::Float64=NaN)
    
    # Check if all parameters are NaN ########
    nancond = .&(isnan(iota), isnan(kappa), isnan(lambda), isnan(sigmah))
    if .&(isempty(svm_dict), isempty(bt.mi._svm_dict))
        println("Parameter dictionary not found. Returning...")
        return 
    # elseif .&(isempty(svm_dict), nancond)
    #     println("Please define at least one parameter value. Returning...")
    #     return
    end
    # ########################################
    
    # Set Parameter Values: ##################
    if isempty(svm_dict)
        if display_msgs
            println("Setting parameter dictionary to batch object's parameter dictionary")
        end
        svm_dict = bt.mi._svm_dict
    end
    
    if !nancond
        if display_msgs
            println("Setting parameter values... ")
        end
        if !isnan(iota)
            svm_dict[:iota] = eval(iota)
        end
        if !isnan(kappa)
            svm_dict[:kappa] = eval(kappa)
        end
        if !isnan(lambda)
            svm_dict[:lambda] = eval(lambda)
        end
        if !isnan(sigmah)
            svm_dict[:sigmah] = eval(sigmah)
        end
    end
    # ########################################
    
    # Find row in the Parameter Combinations DataFrame
    svm_cols = [:lambda, :sigmah]
    cols = [x for x in names(bt.bp.df) if !(x in vcat([:comb_num, :m_comb_num], svm_cols))]
    LL = [abs.(bt.bp.df[col] .- svm_dict[col]) .< 1e-5 for col in cols]
    for col in svm_cols
        if isnan(svm_dict[col])
            LL = vcat(LL, [isnan.(bt.bp.df[col])])
        else
            LL = vcat(LL, [abs.(bt.bp.df[col] .- svm_dict[col]) .< 1e-5])
        end
    end

    return DataFrame(bt.bp.df[.&(LL...), [:comb_num, :m, :m_comb_num]])
end

