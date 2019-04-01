module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["ModelObj", "AnalyticFunctions", "BondPrInterp", "EqFinDiff"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


module TestFunctions

using Distributed
using AnalyticFunctions: get_cvm_vb
using BondPrInterp: get_cvm_bond_price,
    get_cvm_debt_price,
    get_svm_bond_price,
    get_svm_debt_price
using EqFinDiff: get_cvm_eq, eq_fd


function homogeneity_tests(svm, phi, c, p; model="cvm", fun="bond", Vt=NaN, ttm=NaN, mu_b=NaN)
    if isnan(mu_b)
        mu_b = svm.mu_b
    end
    
    if isnan(Vt)
        Vt = svm.pm.V0
    end
    
    if model== "cvm"
        if fun=="bond"            
            if isnan(ttm)
                ttm=svm.m
            end
            v1 = get_cvm_bond_price(svm, ttm, svm.pm.sigmal, 
                                               Vt=Vt, mu_b=mu_b, c=c, p=p)
            v2 = (1/phi) * get_cvm_bond_price(svm, ttm, svm.pm.sigmal, Vt=Vt, 
                                              mu_b=mu_b/phi, c=phi*c, p=phi*p)    
        elseif fun=="vb"
            v1 = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            v2 = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b/phi,
                                                c=phi * c, p= phi* p)
        elseif fun=="debt"
            vb = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            
            v1 = get_cvm_debt_price(svm, vb, svm.pm.sigmal; 
                                    Vt=Vt, mu_b=mu_b, c=c, p=p)
            v2 = get_cvm_debt_price(svm, vb, svm.pm.sigmal; 
                                    Vt=Vt, mu_b=mu_b/phi,
                                    c=phi * c, p= phi* p)
        elseif fun=="equity"            
            if isnan(ttm)
                ttm=svm.m
            end
            v1 = get_cvm_eq(svm, Vt, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            v2 = get_cvm_eq(svm, Vt, svm.pm.sigmal; mu_b=mu_b/phi,
                                      c=phi * c, p= phi* p)
        end
    elseif model=="svm"
        if fun=="bond"            
            if isnan(ttm)
                ttm=svm.m
            end 
            
            vbl = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            
            v1 = get_svm_bond_price(svm, vbl, ttm; Vt=Vt, mu_b=mu_b, c=c, p=p)
            v2 = (1/phi) * get_svm_bond_price(svm, vbl, ttm; Vt=Vt, 
                                              mu_b=mu_b/phi, c=phi*c, p=phi*p)    
        elseif fun=="debt"
            vbl = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            
            v1 = get_svm_debt_price(svm, vbl; Vt=Vt, mu_b=mu_b, c=c, p=p)
            v2 = get_svm_debt_price(svm, vbl; Vt=Vt, mu_b=mu_b,
                                              c=phi * c, p= phi* p)
        elseif fun=="equity"            
            vbl = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
            
            df1 = eq_fd(svm, vbl; mu_b=mu_b, c=c, p=p)
            df2 = eq_fd(svm, vbl; mu_b=mu_b/phi,
                                  c=phi * c, p= phi* p)
            
            v1 = df1.equity[1]
            v2 = df2.equity[1]
        end
    end

    return v1 - v2
end


function cvm_svm_comparison(svm, c, p; fun="bond", Vt=NaN, ttm=NaN, mu_b=NaN)
    if isnan(mu_b)
        mu_b = svm.mu_b
    end
    
    if isnan(Vt)
        Vt = svm.pm.V0
    end

    vbl = get_cvm_vb(svm, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)

    if fun=="bond"            
        if isnan(ttm)
            ttm=svm.m
        end
        v1 = get_svm_bond_price(svm, vbl, ttm; Vt=Vt, mu_b=mu_b, c=c, p=p)
        v2 = get_cvm_bond_price(svm, ttm, svm.pm.sigmal, Vt=Vt, 
                                mu_b=mu_b, c=c, p=p)
    elseif fun=="debt"
        v1 = get_svm_debt_price(svm, vbl; Vt=Vt, mu_b=mu_b, c=c, p=p)
        v2 = get_cvm_debt_price(svm, vbl, svm.pm.sigmal; 
                                Vt=Vt, mu_b=mu_b, c=c, p=p)
    elseif fun=="equity"            
        df1 = eq_fd(svm, vbl; mu_b=mu_b, c=c, p=p)
        v1 = df1.equity[1]
        v2 = get_cvm_eq(svm, Vt, svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
    end

    return v1 - v2
end


end
