
function joint_firm_constructor(sf::Firm, rf::Firm;
                                sf_comb_num::Int64=0,
                                rf_comb_num::Int64=0,
                                # jks::JointKStruct=JointKStruct(fill(NaN, 10)...),
                                jks=JointKStruct(fill(NaN, 10)...),
                                m::Float64=NaN,
                                firm_obj_fun::Symbol=:firm_value,
                                load_results_dfs::Bool=false,
                                cvmdf::DataFrame=DataFrame(),
                                svmdf::DataFrame=DataFrame(),                               
                                opt_k_struct_df_name::String=opt_k_struct_df_name,
                                recompute_svm::Bool=false)

    # Form Batch Objects ###################################
    sf_bt = BatchObj(; model=sf.model)
    rf_bt = BatchObj(; model=rf.model)
    if sf_comb_num > 0
        sf_bt = set_par_dict(sf_bt; comb_num=sf_comb_num)
    end
    if rf_comb_num > 0
        rf_bt = set_par_dict(rf_bt; comb_num=rf_comb_num)
    end
    # ######################################################
    
    if .&(isnan(m), !isnan(jks.m))
        m = jks.m
    end

    if load_results_dfs
        #cvm_bt = BatchObj(; model="cvm")
        svm_bt = BatchObj(; model="svm")

        cvmdf = load_cvm_opt_results_df(; m=m, firm_obj_fun=firm_obj_fun)
        svmdf = load_svm_opt_results_df(svm_bt; m=m,
                                        firm_obj_fun=firm_obj_fun,
                                        opt_k_struct_df_name=opt_k_struct_df_name,
                                        recompute=recompute_svm)
    end

    # Set Optimal Capital Structure ########################
    df = (sf.model == "cvm") ? cvmdf : svmdf
    if .&(sf_comb_num > 0, !isempty(df))
        sf = set_opt_k_struct(sf, df)
    end
    
    df = (rf.model == "cvm") ? cvmdf : svmdf
    if .&(rf_comb_num > 0, !isempty(df))
        rf = set_opt_k_struct(rf, df)
    end
    # #####################################################
    

    jf = JointFirms(jks, sf, rf, sf_bt, rf_bt, cvmdf, svmdf)
    #jf = JointFirms(jks, sf, rf, cvm_bt, svm_bt, cvmdf, svmdf)
end


function otc_constructor(jep, sf_bt, rf_bt; otc_m::Float64=NaN)
    # Adjust parameter dictionaries
    for var in [:alpha, :pi, :r, :gross_delta, :xi, :sigmal]
        sf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
        rf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
    end
    
    # Form OTC Safe Firm
    otc_sf_comb_num = get_batch_comb_num(sf_bt;
                                         iota=jep.sfp.iota,
                                         kappa = jep.kotc,
                                         lambda = jep.sfp.lambda,
                                         sigmah = jep.sfp.sigmah)[1, :comb_num]
    _, otc_sf = get_bt_mobj(;model=sf_bt.model, comb_num=otc_sf_comb_num)

    # Form OTC Risky Firm
    otc_rf_comb_num = get_batch_comb_num(rf_bt;
                                         iota=jep.rfp.iota,
                                         kappa = jep.kotc,
                                         lambda = jep.rfp.lambda,
                                         sigmah = jep.rfp.sigmah)[1, :comb_num]
    _, otc_rf = get_bt_mobj(;model=rf_bt.model, comb_num=otc_rf_comb_num)


    # Form Joint Equilibrium Firm Object
    otc = joint_firm_constructor(otc_sf, otc_rf;
                                 m=otc_m,
                                 load_results_dfs=true)

    # Extract Results ####################################################
    sfdf = sf_bt.model == "cvm" ? otc.cvmdf : otc.svmdf
    rfdf = rf_bt.model == "cvm" ? otc.cvmdf : otc.svmdf

    # Set Optimal Capital Structure
    otc.sf = set_opt_k_struct(otc.sf, sfdf)
    otc.rf = set_opt_k_struct(otc.rf, rfdf)

    sfdf = sfdf[sfdf[:comb_num] .==otc_sf_comb_num, :]
    rfdf = rfdf[rfdf[:comb_num] .==otc_sf_comb_num, :]

    return OTCStruct(otc_sf, otc_rf, sfdf, rfdf)
end


function ep_otc_constructor(mu_s::Float64,
                            kep::Float64,
                            kotc::Float64;
                            jep=empty_jep,
                            sf_obj_fun::Symbol=:firm_value,
                            rf_obj_fun::Symbol=:firm_value,
                            run_pool_eq::Bool=true,
                            run_sep_eq::Bool=true,                                              
                            sfp=FirmSpecificParams(fill(NaN, 3)...),
                            rfp=FirmSpecificParams(fill(NaN, 3)...),
                            s_iota::Float64=NaN,
                            s_lambda::Float64=NaN,
                            s_sigmah::Float64=NaN,
                            r_iota::Float64=NaN,
                            r_lambda::Float64=NaN,
                            r_sigmah::Float64=NaN,
                            ep_jks=JointKStruct(fill(NaN, 10)...),
                            ep_m::Float64=NaN,
                            ep_c::Float64=NaN,
                            ep_p::Float64=NaN,
                            otc_m::Float64=NaN)

    # #################################################################
    # Parameters ######################################################
    # #################################################################
    jep = store_joint_eq_parameters(mu_s, kep, kotc;
                                    jep=jep,
                                    sfp=sfp, rfp=rfp,
                                    s_iota=s_iota,
                                    s_lambda=s_lambda,
                                    s_sigmah=s_sigmah,
                                    r_iota=r_iota,
                                    r_lambda=r_lambda,
                                    r_sigmah=r_sigmah)

    # Check for missing parameters
    # Lambda and sigmah can be NaN => CVM model
    missing = false
    if any([isnan(getfield(jep, x)) for x in [:mu_s, :kep, :kotc]])
        println("Missing Market Parameters!")
        missing = true
    end
    if any([isnan(getfield(jep.sfp, :iota)) for x in fieldnames(FirmSpecificParams)])
        println("Missing Safe Firm's specific parameters!")
        missing = true
    end
    if any([isnan(getfield(jep.rfp, :iota)) for x in fieldnames(FirmSpecificParams)])
        println("Missing Risky Firm's specific parameters!")
        missing = true
    end

    # Output
    if missing
        return
    end
    # #################################################################

    
    # #################################################################
    # Safe and Risky Firms' Models ####################################
    # #################################################################
    cbt = get_bt(;comb_num=1, model="cvm")
    sbt = get_bt(;comb_num=1, model="svm")
    
    # Identify Firm Model
    sf_model = .&(!isnan(jep.sfp.iota),!isnan(jep.sfp.lambda)) ? "svm" : "cvm"
    rf_model = .&(!isnan(jep.rfp.iota),!isnan(jep.rfp.lambda)) ? "svm" : "cvm"

    # Set Models
    sf_bt = sf_model == "cvm" ? cbt : sbt
    rf_bt = rf_model == "cvm" ? cbt : sbt
    # #################################################################   

    
    # #################################################################
    # Electronic Platform ############################################# 
    # #################################################################
    # Measure of Firms and Standardized Bond
    ep_jks = store_ep_params(jep.mu_s;
                             ep_m=ep_m,
                             ep_c=ep_c,
                             ep_p=ep_p)

    # Check for missing parameters
    if any([isnan(getfield(ep_jks, x)) for x in [:mu_s, :m, :c, :p]])
        println("Missing Electronic Platform parameters")
        return
    end

    ep = ep_constructor(jep, sf_bt, rf_bt;
                        ep_jks=ep_jks,
                        ep_m=ep_m,
                        ep_c=ep_c,
                        ep_p=ep_p,
                        run_pool_eq=run_pool_eq,
                        run_sep_eq=run_sep_eq,                       
                        sf_obj_fun=sf_obj_fun,
                        rf_obj_fun=rf_obj_fun)
    # #################################################################

    
    # #################################################################   
    # Over-the-Counter Markets ########################################   
    # #################################################################
    otc = otc_constructor(jep, sf_bt, rf_bt; otc_m=otc_m)
    # #################################################################

    return JointEquilibrium(jep, ep, otc)
end

