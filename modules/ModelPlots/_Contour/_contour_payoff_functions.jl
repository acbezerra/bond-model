


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

    # Get RMP-Specific Payoffs
    fi_rm_fv_fun, fi_nrm_fv_fun = get_rm_payoff_funs(fi_fd, xvar, yvar, :firm_value)
    fi_rm_mbr_fun, fi_nrm_mbr_fun = get_rm_payoff_funs(fi_fd, xvar, yvar, :MBR)
    
    # RM Payoff > NRM Payoff ?
    fi_rm_cond_fun = (x, y) -> fi_rm_fv_fun(x, y) >= fi_nrm_fv_fun(x, y)

    
    # Payoff is the Maximum RMP-Conditional Firm Value
    fi_fv_fun = (x, y) -> maximum([fi_fd[:firm_value][xvar](x), 
                                   fi_fd[:firm_value][yvar](y)])

    fi_mbr_fun = (x, y) -> fi_rm_cond_fun(x, y) ? fi_rm_mbr_fun(x, y) : fi_nrm_mbr_fun(x,y)

    return fi_rm_cond_fun, fi_fv_fun, fi_mbr_fun
end


function misrep_payoff_functions(fi_fd::Dict, mp_fd::Dict;
                                 xvar::Symbol=:iota, yvar::Symbol=:sigmah)


    # Compute Type-Specific MBR under Full Information ################
    # If (x, y) type chooses Risk-Management, return MBR associated
    # with Risk-Management; else, return the MBR associated with
    # No-Risk-Management 
    fi_rm_cond_fun, _ = fi_payoff_functions(fi_fd; xvar=xvar, yvar=yvar)
    r_fi_rm_fun,  r_fi_nrm_fun = get_rm_payoff_funs(fi_fd, xvar, yvar, :MBR)
    r_fi_mbr_fun = (x, y) -> fi_rm_cond_fun(x, y) ? r_fi_rm_fun(x, y) : r_fi_nrm_fun(x, y)
    # #################################################################

    
    # Compute MBR of Misrepresenting Risky Types #######################
    r_mp_mbr_fun = (x, y) -> maximum([mp_fd[:r_MBR][xvar](x), 
                                      mp_fd[:r_MBR][yvar](y)])
    # #################################################################

    
    r_mbr_diff_fun = (x, y) -> r_mp_mbr_fun(x,y) - r_fi_mbr_fun(x,y)

    return r_fi_mbr_fun, r_mp_mbr_fun, r_mbr_diff_fun
end


function joint_eq_payoff_functions(jfd::Dict; eq_type::String="pooling",
                                   xvar::Symbol=:iota, yvar::Symbol=:sigmah)

    # rm_cond = check_xy_variables(xvar, yvar)
    # if !rm_cond
    #     return
    # end
    

    r_obj_fun = :MBR
    if eq_type == "separating"
        r_obj_fun = :firm_value
    end

    # if xvar == :iota
    #     r_rm_cond_fun = (x, y) -> jfd[Symbol(:r_, r_obj_f)][xvar](x) > jfd[Symbol(:r_, r_obj_f)][yvar](y)
    # else
    #     r_rm_cond_fun = (x, y) -> jfd[Symbol(:r_, r_obj_f)][xvar](x) < jfd[Symbol(:r_, r_obj_f)][yvar](y)
    # end

    r_rm_payoff_fun, r_nrm_payoff_fun = get_rm_payoff_funs(jfd, xvar, yvar, r_obj_fun)
    r_rm_cond_fun = (x, y) -> r_rm_payoff_fun(x, y) > r_nrm_payoff_fun(x, y)
    

    r_jeq_obj_fun = (x, y) -> maximum([jfd[Symbol(:r_, r_obj_fun)][xvar](x), 
                                       jfd[Symbol(:r_, r_obj_fun)][yvar](y)])
    
    s_jeq_fv_fun = (x, y) -> r_rm_cond_fun(x, y) ? jfd[:s_firm_value][xvar](x) : jfd[:s_firm_value][yvar](y)


    return r_rm_cond_fun, r_jeq_obj_fun, s_jeq_fv_fun
end




# function get_payoff(df::DataFrame, fd::Dict;
#                     ftype::Symbol=:safe,
#                     xvar::Symbol=:iota,
#                     yvar::Symbol=:sigmah,
#                     zvar::Symbol=:firm_value)

    
#     z_fun = (x, y) -> maximum([fd[zvars][xvar](x), 
#                                fd[zvars][yvar](y)])

#     # Choice of Risk-Management Policy
#     eq_type = df[1, :eq_type]
#     if eq_type == "full_info" #, any([xvar == :iota, yvar == :iota]))
#         if xvar == :iota
#             rm_bool = (x, y) -> fd[:firm_value][xvar](x) > fd[:firm_value][yvar](y)
#         elseif yvar == :iota
#             rm_bool = (x, y) -> fd[:firm_value][xvar](x) < fd[:firm_value][yvar](y)
#         else
#             rm_bool = (x, y) -> NaN
#         end
#     elseif eq_type == "pooling"
        

    

#             rm_var = (xvar == :iota) ? xvar : yvar
#         nrm_var = (xvar == :iota) ? yvar : xvar 

# end
