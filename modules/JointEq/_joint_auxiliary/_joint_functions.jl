
function store_joint_eq_parameters(mu_s::Float64,
                                   kep::Float64,
                                   kotc::Float64;
                                   jep=empty_jep,
                                   sfp=FirmSpecificParams(fill(NaN, 3)...),
                                   rfp=FirmSpecificParams(fill(NaN, 3)...),
                                   s_iota::Float64=NaN,
                                   s_lambda::Float64=NaN,
                                   s_sigmah::Float64=NaN,
                                   r_iota::Float64=NaN,
                                   r_lambda::Float64=NaN,
                                   r_sigmah::Float64=NaN)

    # Set Market Parameters
    setfield!(jep, :mu_s, mu_s)
    setfield!(jep, :kep, kep)
    setfield!(jep, :kotc, kotc)

    # Conditional Loading
    if any([!isnan(getfield(sfp, :iota)),
            !isnan(getfield(sfp, :lambda)),
            !isnan(getfield(sfp, :sigmah))])
        setfield!(jep, :sfp, sfp)
    end
    if any([!isnan(getfield(rfp, :iota)),
            !isnan(getfield(rfp, :lambda)),
            !isnan(getfield(rfp, :sigmah))])
        setfield!(jep, :rfp, rfp)
    end
    
    # Set Firm Specific Parameters
    if !isnan(s_iota)
        setfield!(jep.sfp, :iota, s_iota)
    end
    if !isnan(s_lambda)
        setfield!(jep.sfp, :lambda, s_lambda)
    end
    if !isnan(s_sigmah)    
        setfield!(jep.sfp, :sigmah, s_sigmah)
    end
    if !isnan(r_iota)
        setfield!(jep.rfp, :iota, r_iota)
    end
    if !isnan(r_lambda)
        setfield!(jep.rfp, :lambda, r_lambda)
    end
    if !isnan(r_sigmah)    
        setfield!(jep.rfp, :sigmah, r_sigmah)
    end

    # Set Common Parameters
    fcp = FirmCommonParams(common_params[:V0],
                           common_params[:alpha],
                           common_params[:pi],
                           common_params[:r],
                           svm_param_values_dict[:gross_delta][1],
                           svm_param_values_dict[:xi][1],
                           svm_param_values_dict[:sigmal][1])
    setfield!(jep, :fcp, fcp)


    return jep
end


function store_ep_params(mu_s;
                         ep_jks=JointKStruct(fill(NaN, 10)...),
                         ep_m::Float64=NaN,
                         ep_c::Float64=NaN,
                         ep_p::Float64=NaN)

    setfield!(ep_jks, :mu_s, mu_s)
    
    if !isnan(ep_m)
        setfield!(ep_jks, :m, ep_m)
    end
    if !isnan(ep_c)
        setfield!(ep_jks, :c, ep_c)
    end
    if !isnan(ep_p)
        setfield!(ep_jks, :p, ep_p)
    end

    return ep_jks
end        
    

function check_param_consistency(jf)
    # Common Parameters  
    parvec = [par for par in [:r, :xi, :kappa, :alpha] if
              abs.(get_param(jf.sf, par) - get_param(jf.rf, par)) > 1e-6]
            
    if !isempty(parvec)
        [println(string(par, ": values do not match!")) for par in parvec]
        return false
    end
    
    return true
end


function find_risky_combinations(svm::Firm; 
                                 cvmdf::DataFrame=DataFrame(), 
                                 svmdf::DataFrame=DataFrame())
    if isempty(svmdf)
        svmdf = get_bt(;model="svm", comb_num=1).bp.df
    end
    
    # Reorder columns
    cols = vcat([:match, :model], names(svmdf))

    # param_columns = [:m, :xi, :kappa, :gross_delta,
    #                  :lambda, :iota, :sigmal, :sigmah]

    param_columns = [x for x in names(svmdf) if
                     !(x in [:comb_num, :m_comb_num])]
    
    # Match SVM Values
    svm_common_cond = .&([abs.(svmdf[fn] .- get_param(svm, fn)) .< 1e-6
                          for fn in [:mu_b, :m, :kappa]]...)
   
    if get_obj_model(svm) == "cvm"
        if isempty(cvmdf)
            cvmdf = get_bt(;model="cvm", comb_num=1).bp.df
        end
        
        # Find CVM Combinations
        cvm_common_cond = .&([abs.(cvmdf[fn] .- get_param(svm, fn)) .< 1e-6
                              for fn in [:mu_b, :m, :kappa]]...)
        cvm_iota_cond = cvmdf[:iota] .>= get_param(svm, :iota)
        cvm_match = cvmdf[.&([cvm_common_cond, cvm_iota_cond]...), :]
        cvm_match[:model] = "cvm"

        # Find Own Combination 
        cvm_match[:match] = .&([abs.(cvm_match[fn] .- get_param(svm,fn)) .< 1e-6
                                for fn in param_columns
                                if !(fn in [:lambda, :sigmah])]...)
        
        # Find SVM Combinations
        svm_match = svmdf[svm_common_cond, :]
        svm_match[:model] = "svm"
        svm_match[:match] = false
        
        return vcat([cvm_match, svm_match]...)[cols]
    else
        # Find SVM Combinations
        svm_shock_cond = .&([svmdf[fn] .>= get_param(svm, fn) 
                             for fn in [:lambda, :sigmah]]...)

        svm_match = svmdf[.&([svm_common_cond,
                              svm_shock_cond]...), :]
        svm_match[:model] = "svm"

        # Find Own Combination 
        svm_match[:match] = .&([abs.(svm_match[fn] .- get_param(svm,fn)) .< 1e-6
                                for fn in param_columns]...)
        
        return svm_match[cols]
    end
end


# ########################################################################
function check_firm_dict(id)
    for idk in [:comb_num, :m_comb_num]
        if !haskey(id, idk) 
            id[idk] = 0
        end
    end    
    if !(:m in keys(id))
        id[:m] = NaN
    end
    
    # Must have at least one valid identifier:
    id_cond = (id[:comb_num] > 0) || .&(!isnan(id[:m]), id[:m_comb_num] > 0)
    if !id_cond
        return id, false
    end
    return id, true
end

# kappa must match!
function check_risky_firm_dict(bts, idr)
    idr, id_cond1 = check_firm_dict(idr)
    
    # Check if has at at least one valid identifier:
    id_cond2 = sum([x in keys(idr) for x in [:kappa, :lambda, :sigmah]]) > 0 
    if id_cond1
        println("Risky Firm will be constructed from identifiers!")
        return idr, true
    elseif id_cond2
        println("Risky Firm will be constructed by modifying Safe Firm's parameters!")
        for idk in [:kappa, :lambda, :sigmah]
            if !haskey(idr, idk)
                idr[idk] = NaN
            end
        end
        combdf =  get_batch_comb_num(bts;
                                    kappa=idr[:kappa],
                                    lambda=idr[:lambda],
                                    sigmah=idr[:sigmah])
        idr[:comb_num] = combdf[:comb_num][1]

        return idr, true
    end
    return idr, false
end


function round_value(x::Float64)
    xr_vec = [floor(x), (ceil(x) + floor(x))/2, ceil(x)]
    diffs = [x-xr_vec[1], abs(x-xr_vec[2]), xr_vec[3] - x]    
    return xr_vec[argmin(diffs)]
end

