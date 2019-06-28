
LL = [ ]
# for sf_comb_num in cvm_combs[1]

# Preliminary Objects #####################################
sf_bt, sf = Batch.get_bt_mobj(; model=sf_model, comb_num=sf_comb_num)
sf_df = (sf_model == "cvm") ? cvmdf : svmdf
sf = ModelObj.set_opt_k_struct(sf, sf_df)

# Capital Structure -> Fixed
jks = JointEq.JointKStruct(1., 
                           sf.optKS.mu_b,
                           sf.optKS.m, sf.optKS.c, sf.optKS.p, 
                           NaN, NaN, NaN, NaN, NaN)

# rf_bt, rf = Batch.get_bt_svm(; comb_num=svm_combs[1])
rf_bt, rf = Batch.get_bt_svm(; comb_num=svmdf[1, :comb_num])

# Joint Firm Object
jf = JointEq.JointFirms(jks, sf, rf,
                        sf_bt,
                        rf_bt,  
                        cvmdf, svmdf)    
# #########################################################


for rf_comb_num in svmdf[:comb_num] #svm_combs
    rf_bt, rf = Batch.get_bt_svm(; comb_num=rf_comb_num)
    
    # # Joint Firm Object
    # jf = JointEq.JointFirms(jks, sf, rf,
    #                         sf_bt,
    #                         rf_bt,  
    #                         cvmdf, svmdf)
    
    # Joint Equilibrium Parameters
    jep = JointEq.store_joint_eq_parameters(jks.mu_s, sf.pm.kappa, sf.pm.kappa;
                                            s_iota=sf.pm.iota,
                                            s_lambda=sf.pm.lambda,
                                            s_sigmah=sf.pm.sigmah,
                                            r_iota=rf.pm.iota,
                                            r_lambda=rf.pm.lambda,
                                            r_sigmah=rf.pm.sigmah)

    
    # Compute Misrepresentation
    jeq = JointEq.ep_constructor(jep, jf.cvm_bt, jf.svm_bt;
                                 ep_jks=jks,
                                 run_pool_eq=false,
                                 run_sep_eq=false,                       
                                 sf_obj_fun=firm_obj_fun,
                                 rf_obj_fun=firm_obj_fun,
                                 rerun_full_info=false,
                                 run_misrep=true,
                                 rerun_pool=false,
                                 rerun_sep=false)
    
    
    push!(LL, getfield(jeq, :misrep))
end

return vcat(LL...)
# end

# return jf, cvm_combs, svm_combs, LL


