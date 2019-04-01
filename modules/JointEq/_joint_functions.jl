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
