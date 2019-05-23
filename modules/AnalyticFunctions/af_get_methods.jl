function get_k_struct(svm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN, p::Float64=NaN)
    if isnan(mu_b)
        mu_b = svm.mu_b
    end

    if isnan(m)
        m=svm.m
    end
    
    if isnan(c)
        c=svm.c
    end

    if isnan(p)
        p=svm.p
    end
    
    return mu_b, m, c, p
end


function get_rgrow(svm)
   return rgrow(svm.pm.r, svm.pm.gross_delta) 
end


function get_rdisc(svm)
   return rdisc(svm.pm.r, svm.pm.xi, svm.pm.kappa) 
end


function get_cvm_vb(svm, sigma;
                    mu_b::Float64=NaN, m::Float64=NaN,
                    c::Float64=NaN, p::Float64=NaN)
    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)

    return cvm_vb(mu_b, m, c, p, 
                  sigma, svm.pm.r, svm.pm.gross_delta, 
                  svm.pm.iota, svm.pm.xi, svm.pm.kappa,
                  svm.pm.alpha, svm.pm.pi)
end


function get_agg_c(svm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN)
    mu_b, m, c, _ = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=NaN)
   
    return mu_b * c * m
end


function get_agg_p(svm; mu_b::Float64=NaN, m::Float64=NaN, p::Float64=NaN)
    mu_b, m, _, p = get_k_struct(svm; mu_b=mu_b, m=m, c=NaN, p=p)
   
    return mu_b * p * m
end


function get_param(svm, pname::Symbol)
    val = NaN
    # svm_pos = findin([string(x) for x in fieldnames(svm)], [pname])
    # svm_pm_pos = findin([string(x) for x in fieldnames(svm.pm)], [pname])
    
    svm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm))])
    svm_pm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm.pm))])
    
    if pname == :C
        return get_agg_c(svm; mu_b=svm.mu_b, m=svm.m, c=svm.c)
    elseif pname == :P
        return get_agg_p(svm; mu_b=svm.mu_b, m=svm.m, p=svm.p)
    elseif pname == :delta
        return svm.pm.gross_delta - svm.pm.iota
    else
        return extract_param(svm, pname)
    end
    
    # elseif length(svm_pos) > 0
    #     # return getfield(svm, fieldnames(svm)[svm_pos[1]])
    #     return getfield(svm, fieldnames(typeof(svm))[svm_pos[1]])
    # else
    #     # return getfield(svm.pm, fieldnames(svm.pm)[svm_pm_pos[1]])
    #     return getfield(svm.pm, fieldnames(typeof(svm.pm))[svm_pm_pos[1]])
    # end
end


function get_leverage(debt::Union{Float64, Array{Float64, 1}},
                      equity::Union{Float64, Array{Float64, 1}})
    return (debt ./ (debt .+ equity)) .* 100
end

function get_mbr(V0::Float64,
                 debt::Union{Float64, Array{Float64, 1}},
                 equity::Union{Float64, Array{Float64, 1}})
    return (equity ./ (V0 .- debt) .- 1.) .* 100
end
