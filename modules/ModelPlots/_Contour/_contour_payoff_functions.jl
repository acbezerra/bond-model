


function check_rm_vars(xvar::Symbol, yvar::Symbol)
    if .&(xvar != :iota, yvar != :iota)
        println("Either xvar or yvar must be :iota. Exiting...")
        return false
    end

    return true
end


function get_rm_payoff_funs(fd::Dict, xvar::Symbol, yvar::Symbol, obj_fun::Symbol)   
    rm_vars_present = check_rm_vars(xvar, yvar)
    if !rm_vars_present
        return
    end

    if xvar == :iota
        rm_fun = (x, y) -> fd[obj_fun][xvar](x)
        nrm_fun = (x, y) -> fd[obj_fun][yvar](y)
    else
        nrm_fun = (x, y) -> fd[obj_fun][xvar](x)
        rm_fun = (x, y) -> fd[obj_fun][yvar](y)
    end

    return rm_fun, nrm_fun
end



function fi_payoff_functions(fi_fd::Dict;
                             xvar::Symbol=:iota, yvar::Symbol=:sigmah)

    fi_funs = Dict{Symbol, Any}(:xvar => xvar,
                                :yvar => yvar)
     
    # Firm Value ##################################################
    # Get RMP-Specific Payoffs
    fi_funs[:rm_fv], fi_funs[:nrm_fv] = get_rm_payoff_funs(fi_fd, xvar, yvar, :firm_value)

    
    # RM Payoff > NRM Payoff ?
    fi_funs[:rm_cond] = (x, y) -> fi_funs[:rm_fv](x, y) >= fi_funs[:nrm_fv](x, y)

    
    # Payoff is the Maximum RMP-Conditional Firm Value
    fi_funs[:fv] = (x, y) -> maximum([fi_fd[:firm_value][xvar](x), 
                                      fi_fd[:firm_value][yvar](y)])
    # #############################################################

    for zvar in [z for z in contour_zvars if z != :firm_value]
        zsym = contour_zvars_sym[zvar]
        fi_funs[Symbol(:rm_, zsym)], fi_funs[Symbol(:nrm_, zsym)] = get_rm_payoff_funs(fi_fd, fi_fd[:xvar], fi_fd[:yvar], zvar)
        fi_funs[zsym] = (x, y) -> fi_funs[:rm_cond](x, y) ? fi_funs[Symbol(:rm_, zsym)](x, y) : fi_funs[Symbol(:nrm_, zsym)](x, y)
    end
    
    return fi_funs
end


function misrep_payoff_functions(fi_funs::Dict{Symbol, Any}, mp_fd::Dict{Symbol, Any};
                                 xvar::Symbol=:iota, yvar::Symbol=:sigmah)

    if any([fi_funs[:xvar] != xvar, fi_funs[:yvar] != yvar])
        println("Full Info and Misrep x and y variables do not coincide. Exiting...")
        return
    end

    mp_funs = Dict{Symbol, Any}(:xvar => xvar, :yvar => yvar)

    # Compute Type-Specific MBR under Full Information ################
    mp_funs[:rm_mbr], mp_funs[:nrm_mbr] = get_rm_payoff_funs(mp_fd, mp_fd[:xvar], mp_fd[:yvar], :MBR)
    mp_funs[:rm_cond] = (x, y) -> mp_funs[:rm_mbr](x, y) >= mp_funs[:nrm_mbr](x, y)
    mp_funs[:mbr] = (x, y) -> maximum([mp_fd[:r_MBR][xvar](x), 
                                       mp_fd[:r_MBR][yvar](y)])
    mp_funs[:mp_fi_mbr_diff] = (x, y) -> mp_funs[:mbr](x, y) - fi_funs[:mbr](x, y)
    

    for zvar in [z for z in contour_zvars if z != :MBR]
        zsym = contour_zvars_sym[zvar]
        mp_funs[Symbol(:rm_, zsym)], mp_funs[Symbol(:nrm_, zsym)] = get_rm_payoff_funs(mp_fd, mp_fd[:xvar],
                                                                                       mp_fd[:yvar], zvar)
        mp_funs[zsym] = (x, y) -> mp_funs[:rm_cond](x, y) ? mp_funs[Symbol(:rm_, zsym)](x, y) : mp_funs[Symbol(:nrm_, zsym)](x, y)
        mp_funs[Symbol(:mp_fi_, zsym, :_diff)] = (x, y) -> mp_funs[zsym](x, y) - fi_funs[zsym](x, y)
    end
    
    return mp_funs 
end


function jeq_payoff_functions(fi_funs::Dict{Symbol, Any}, jfd::Dict;
                              eq_type::String="pooling",
                              xvar::Symbol=:iota, yvar::Symbol=:sigmah)


    if any([fi_funs[:xvar] != xvar, fi_funs[:yvar] != yvar])
        println("Full Info and Misrep x and y variables do not coincide. Exiting...")
        return
    end

    jeq_funs = Dict{Symbol, Any}(:xvar => xvar, :yvar => yvar,
                                 :safe => Dict{Symbol, Any}(),
                                 :risky => Dict{Symbol, Any}())

    r_obj_fun = :MBR
    if eq_type == "separating"
        r_obj_fun = :firm_value
    end

    
    # Set objective function symbol
    zsym = contour_zvars_sym[r_obj_fun]

    # Risky Firm's Objective Function Payoff in case of Risk-Management v.s. No Risk-Management
    jeq_funs[:risky][Symbol(:rm_, zsym)], jeq_funs[:risky][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar,
                                                                                                    Symbol(:r_, r_obj_fun))

    # Choose Risk-Management if it maximizes Payoff
    jeq_funs[:risky][:rm_cond] = (x, y) -> jeq_funs[:risky][Symbol(:rm_, zsym)](x, y) >= jeq_funs[:risky][Symbol(:nrm_, zsym)](x, y)


    # Risky Firm's Objective Function Payoff
    jeq_funs[:risky][zsym] = (x, y) -> maximum([jfd[Symbol(:r_, r_obj_fun)][xvar](x), 
                                                jfd[Symbol(:r_, r_obj_fun)][yvar](y)])

    # Difference between Joint and Full Information Equilibrium
    jeq_funs[:risky][Symbol(:jeq_fi_, zsym, :_diff)] = (x, y) -> jeq_funs[:risky][zsym](x, y) - fi_funs[zsym](x, y)

    
    # Safe Firm's Payoff Depends on what Risky Firm chooses 
    jeq_funs[:safe][Symbol(:rm_, zsym)], jeq_funs[:safe][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar,
                                                                                                   Symbol(:s_, r_obj_fun))
    jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[:safe][Symbol(:rm_, zsym)](x, y) : jeq_funs[:safe][Symbol(:nrm_, zsym)](x, y)

    # jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:safe][zsym](x, y) - Need FI payoff #fi_funs[zsym](x, y)


    for zvar in [z for z in contour_zvars if z != r_obj_fun]
        zsym2 = contour_zvars_sym[zvar]

        for ft in keys(contour_firm_types)
            ft_z = contour_firm_types[ft]
            jeq_funs[ft][Symbol(:rm_, zsym2)], jeq_funs[ft][Symbol(:nrm_, zsym2)] = get_rm_payoff_funs(jfd, jfd[:xvar],
                                                                                                     jfd[:yvar],
                                                                                                     Symbol(ft_z, zvar))
            jeq_funs[ft][zsym2] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[ft][Symbol(:rm_, zsym2)](x, y) : jeq_funs[ft][Symbol(:nrm_, zsym2)](x, y)
        end
        
        jeq_funs[:risky][Symbol(:jeq_fi_, zsym2, :_diff)] = (x, y) -> jeq_funs[:risky][zsym2](x, y) - fi_funs[zsym2](x, y)
    end
    
                               
    return deepcopy(jeq_funs)
end


function get_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs,
                                     fi_fv::Float64, fi_fv_fun, k_otc::Float64)
    jeq_ind = (x, y) -> mp_funs[:mbr](x, y) >= fi_funs[:mbr](x, y) 
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:safe][:fv](x, y) >= sep_funs[:safe][:fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false
    
    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :mbr => Dict{Symbol, Any}())
    for zvar in [:fv, :mbr, :lev]
        fun_dict[:mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[zvar](x, y) : .0
        fun_dict[:mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[:risky][zvar](x, y) : .0
        fun_dict[:mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[:risky][zvar](x, y) : .0
    end

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) + 
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) + 
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:risky][:mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:risky][:mbr](x, y))
    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv + #fi_funs[:fv](x, y) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:safe][:fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:safe][:fv](x, y))
    fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : fun_dict[:eq_bool](x, y)

    catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)], 
                    [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
    fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]          

    return fun_dict
end




# function joint_eq_payoff_functions(fi_fd::Dict, jfd::Dict; eq_type::String="pooling",
#                                    xvar::Symbol=:iota, yvar::Symbol=:sigmah)


#     fi_funs = fi_payoff_functions(fi_fd; xvar=xvar, yvar=yvar)
#     r_rm_fv_fun, r_nrm_fv_fun = get_rm_payoff_funs(jfd, xvar, yvar, :s_firm_value)
#     r_rm_mbr_fun, r_nrm_mbr_fun = get_rm_payoff_funs(jfd, xvar, yvar, :r_MBR)


#     r_obj_fun = :MBR
#     if eq_type == "separating"
#         r_obj_fun = :firm_value
#     end

#     if r_obj_fun == :MBR
#         r_rm_cond_fun = (x, y) -> r_rm_mbr_fun(x, y) > r_nrm_mbr_fun(x, y)
#         r_jeq_mbr_fun = (x, y) -> maximum([jfd[:r_MBR][xvar](x), 
#                                            jfd[:r_MBR][yvar](y)])

#         r_jeq_fv_fun = (x, y) -> r_rm_cond_fun(x, y) ? r_rm_fv_fun(x, y) : r_nrm_fv_fun(x, y)
#     elseif r_obj_fun == :firm_value
#         r_rm_cond_fun = (x, y) -> r_rm_fv_fun(x, y) > r_nrm_fv_fun(x, y)
#         r_jeq_fv_fun = (x, y) -> maximum([jfd[:r_firm_value][xvar](x), 
#                                            jfd[:r_firm_value][yvar](y)])

#         r_jeq_mbr_fun = (x, y) -> r_rm_cond_fun(x, y) ? r_rm_mbr_fun(x, y) : r_nrm_mbr_fun(x, y)
#     end


#     r_jeq_fi_mbr_diff_fun = (x, y) -> r_jeq_mbr_fun(x,y) - fi_funs[:mbr](x, y) 
        
   
#     s_jeq_fv_fun = (x, y) -> r_rm_cond_fun(x, y) ? jfd[:s_firm_value][xvar](x) : jfd[:s_firm_value][yvar](y)


#     jeq_funs = Dict{Symbol, Any}(:safe => Dict{Symbol, Any}(:fv => s_jeq_fv_fun),
#                                  :risky => Dict{Symbol, Any}(:rm_cond => r_rm_cond_fun,
#                                                              :rm_mbr => r_rm_mbr_fun,
#                                                              :nrm_mbr => r_nrm_mbr_fun,
#                                                              :mbr => r_jeq_mbr_fun,
#                                                              :rm_fv => r_rm_fv_fun,
#                                                              :nrm_fv => r_nrm_fv_fun,
#                                                              :fv => r_jeq_fv_fun,
#                                                              :jeq_fi_mbr_diff => r_jeq_fi_mbr_diff_fun))
                                 
#     return jeq_funs
# end


# unction get_contour_payoff_functions(fidf::DataFrame; xvar=:iota, yvar=:sigmah,
 #                                      zvars::Array{Symbol, 1}=[:MBR, :firm_value, :leverage],
 #                                      jeq_zvars::Array{Symbol, 1}=[:MBR, :firm_value],
 #                                      )
 #    cond = (fidf[:iota] .!= 2.5)

 #    fi_fd = interp_z_values(fidf, xvar, yvar, zvars)
 #    mp_fd = interp_z_values(misrepdf, xvar, yvar, zvars)
 #    pool_fd = interp_z_values(misrepdf, xvar, yvar, zvars)
 #    sep_fd = interp_z_values(misrepdf, xvar, yvar, zvars)
       


 #    fi_funs, mp_funs = misrep_payoff_functions(fi_fd, mp_fd; xvar=xvar, yvar=yvar)
 #    pool_funs = joint_eq_payoff_functions(fi_fd, pool_fd; eq_type="pooling",
 #                                          xvar=xvar, yvar=yvar)



   


