

function find_full_info_vb(svm, jks;
                           mu_b::Float64=NaN,
                           lb::Float64=.75, ub::Float64=1.25,
                           vbN::Int64=15, N::Int64=10^5)

    if isnan(mu_b)
        mu_b = jks.mu_b
    end
    
    if get_obj_model(svm) == "cvm"
        return get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    else
        vbl = get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    
        vbh = get_cvm_vb(svm, svm.pm.sigmal;
                         mu_b=mu_b, m=jks.m,
                         c=jks.c, p=jks.p)

        vbgrid = range(.75 * minimum([vbl, vbh]),
                       stop=1.25 * maximum([vbl, vbh]), length=vbN)

        res = @time fetch(@spawn [eq_fd(svm; vbl=vbl, mu_b=mu_b,
                                        m=jks.m, c=jks.c, p=jks.p)
                                  for vbl in vbgrid])

        eqdf = full_info_eq_deriv_root_search(svm, vcat(res...); N=N)

        return eqdf[1, :vb]
    end
end


function set_full_information_vb!(jf, jks;
                                  rerun_fi_vb::Bool=false,
                                  fi_sf_vb::Float64=NaN,
                                  fi_rf_vb::Float64=NaN,
                                  lb::Float64=.75,
                                  ub::Float64=1.25,
                                  vbN::Int64=20)
    
    if any([.&(isnan(jks.fi_sf_vb), isnan(fi_sf_vb)), rerun_fi_vb])
        fi_sf_vb = find_full_info_vb(jf.sf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_sf_vb, fi_sf_vb)
    elseif !isnan(fi_sf_vb)
         setfield!(jks, :fi_sf_vb, fi_sf_vb)       
    end

    if any([.&(isnan(jks.fi_rf_vb), isnan(fi_rf_vb)), rerun_fi_vb])
        fi_rf_vb = find_full_info_vb(jf.rf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_rf_vb, fi_rf_vb)
    elseif !isnan(fi_rf_vb)
         setfield!(jks, :fi_rf_vb, fi_rf_vb)       
    end

    return jks
end


