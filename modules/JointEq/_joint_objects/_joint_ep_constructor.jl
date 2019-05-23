


function ep_pool_sep_eq(ep_jf, ep_jks,
                        fi_sf_mu_b::Float64,
                        fi_rf_mu_b::Float64,
                        fi_rf_obj_val::Float64;
                        # ep_ks_path::String;
                        equilibrium_type::String="all",
                        sf_obj_fun::Symbol=:firm_value,
                        rf_obj_fun::Symbol=:firm_value,
                        rerun::Bool=true,
                        lb::Float64=.75,
                        ub::Float64=1.25,
                        mu_bN::Int64=20,
                        mu_bN2::Int64=10^5,
                        spline_k::Int64=3,
                        spline_bc::String="extrapolate")

    
    # if equilibrium_type == "pooling"
    #     dfname =  pooldf_name
    # elseif equilibrium_type == "separating"
    #     dfname = sepdf_name
    # else
    #     println(string("Wrong equilibrium type. ",
    #                    "Please enter either <pooling> or <separating>. ",
    #                    "Exiting..."))
    #     return
    # end                
        
    jks = deepcopy(ep_jks)

    # # Joint Optimal mu_b
    # if .&(exists_ep_df(ep_ks_path, dfname), !rerun)
    #     pmdf = CSV.read(string(ep_ks_path, "/", dfname, ".csv");
    #                         types=mps_col_types)

    #     # Slice DataFrame
    #     ep_pm_eqdf = slice_mps_dataframe(pmdf, jks,
    #                                      ep_jf; rerun=rerun)
    # else
        pmdf = find_joint_optimal_bond_measure(ep_jf, jks,
                                               fi_sf_mu_b,
                                               fi_rf_mu_b,
                                               fi_rf_obj_val;
                                               equilibrium_type=equilibrium_type,
                                               sf_obj_fun=sf_obj_fun,
                                               rf_obj_fun=rf_obj_fun,
                                               lb=lb, ub=ub,
                                               mu_bN=mu_bN, mu_bN2=mu_bN2,
                                               spline_k=spline_k,
                                               spline_bc=spline_bc)

        # Reshape
        #ep_pm_eqdf = reshape_sf_rf_df(pmdf)
 #   end
    return pmdf #ep_pm_eqdf
end


function ep_constructor(jep, sf_bt, rf_bt;
                        ep_jks=JointKStruct(fill(NaN, 10)...),
                        ep_m::Float64=NaN,
                        ep_c::Float64=NaN,
                        ep_p::Float64=NaN,
                        run_misrep::Bool=false,
                        run_pool_eq::Bool=true,
                        run_sep_eq::Bool=true,                       
                        sf_obj_fun::Symbol=:firm_value,
                        rf_obj_fun::Symbol=:MBR,
                        fi_fpath_name::String="",
                        rerun_full_info::Bool=true,
                        rerun_pool::Bool=true,
                        rerun_sep::Bool=true,
                        lb::Float64=.75,
                        ub::Float64=1.25,
                        mu_bN::Int64=20,
                        mu_bN2::Int64=10^5,
                        spline_k::Int64=3,
                        spline_bc::String="extrapolate")


    # Measure of Firms and Standardized Bond
    ep_jks = store_ep_params(jep.mu_s;
                             ep_jks=ep_jks,
                             ep_m=ep_m,
                             ep_c=ep_c,
                             ep_p=ep_p)

    # Check for missing parameters
    if any([isnan(getfield(ep_jks, x)) for x in [:mu_s, :m, :c, :p]])
        println("Missing Electronic Platform parameters")
        return
    end

    # Adjust parameter dictionaries
    for var in [:alpha, :pi, :r, :gross_delta, :xi, :sigmal]
        sf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
        rf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
    end

    # Form EP Safe Firm
    ep_sf_comb_num = get_batch_comb_num(sf_bt;
                                        iota=jep.sfp.iota,
                                        kappa=jep.kep,
                                        lambda=jep.sfp.lambda,
                                        sigmah=jep.sfp.sigmah)[1, :comb_num]
    _, ep_sf_svm = get_bt_mobj(;model=sf_bt.model, comb_num=ep_sf_comb_num)

    # Form EP Risky Firm
    ep_rf_comb_num = get_batch_comb_num(rf_bt;
                                        iota=jep.rfp.iota,
                                        kappa=jep.kep,
                                        lambda=jep.rfp.lambda,
                                        sigmah=jep.rfp.sigmah)[1, :comb_num]
    _, ep_rf_svm = get_bt_mobj(;model=rf_bt.model, comb_num=ep_rf_comb_num)
    

    # Joint Firm Constructor ##########################################
    ep_jf = joint_firm_constructor(ep_sf_svm, ep_rf_svm;
                                   jks=ep_jks,
                                   load_results_dfs=false)
    # #################################################################


    # Path to Results Directory
    #jks_fpath = make_jeq_jks_fpath(ep_jf)
    #fi_fname = get_jeq_mus_fname(ep_jf; eq_type="full_info")
    #fi_fpath_name = string(jks_fpath, "/", fi_fname)


    # Electronic Platform Full-Information Equilibrium ################
    if .&(isfile(fi_fpath_name), !rerun_full_info)
        fidf = CSV.read(fi_fpath_name)
                        #types=fidf_col_types)
        
        # Slice DataFrame
        ep_sf_eqdf = DataFrame(fidf[1, :]) # fidf[fidf[:firm_type] .== "safe", :]
        ep_rf_eqdf = DataFrame(fidf[2, :]) # fidf[fidf[:firm_type] .== "risky", :]

        # Capture Full Information Optimal mu_b and MBR ###################
        fi_sf_mu_b = ep_sf_eqdf[1, :mu_b]
        fi_rf_mu_b = ep_rf_eqdf[1, :mu_b]
        fi_sf_obj_val = ep_sf_eqdf[1, sf_obj_fun]
        fi_rf_obj_val = ep_rf_eqdf[1, rf_obj_fun]
        # #################################################################
    elseif rerun_full_info
        ep_sf_eqdf = find_optimal_bond_measure(ep_sf_svm; jks=ep_jks)
        ep_rf_eqdf = find_optimal_bond_measure(ep_rf_svm; jks=ep_jks)

        # Capture Full Information Optimal mu_b and MBR ###################
        fi_sf_mu_b = ep_sf_eqdf[1, :mu_b]
        fi_rf_mu_b = ep_rf_eqdf[1, :mu_b]
        fi_sf_obj_val = ep_sf_eqdf[1, sf_obj_fun]
        fi_rf_obj_val = ep_rf_eqdf[1, rf_obj_fun]
        # #################################################################
    else
        ep_sf_eqdf = DataFrame()
        ep_rf_eqdf = DataFrame()

        fi_sf_mu_b = NaN
        fi_rf_mu_b =  NaN
        fi_sf_obj_val =  NaN
        fi_rf_obj_val =  NaN

        run_pool_eq = false
        run_sep_eq = false
    end
    # #################################################################

            
    # Electronic Platform Misrepresentation ###########################
    # Do risky firms have an incentive to copy the capital structure
    # of the safe firms?
    #ep_ks_path = get_ep_results_path(ep_jf; ep_dir=ep_dir)
    # misrep_fname = get_jeq_mus_fname(ep_jf; eq_type="misrep")
    # misrep_fpath_name = string(jks_fpath, "/", misrep_fname)
    
    misrep_jks = deepcopy(ep_jks)
    setfield!(misrep_jks, :mu_s, 1.)
    setfield!(misrep_jks, :mu_b, fi_sf_mu_b)
    # if .&(exists_ep_df(ep_ks_path, misrepdf_name), !rerun_misrep)
    # if .&(isfile(misrep_fpath_name), !rerun_misrep)
    #     misrepdf = CSV.read(string(jks_fpath, "/", misrep_fname))#;
    #                       #  types=mps_col_types)

    #     # Slice DataFrame
    #     ep_misrep_eqdf = slice_mps_dataframe(misrepdf, misrep_jks,
    #                                          ep_jf; rerun=rerun_misrep)
    if run_misrep
        ep_misrep_eqdf = find_joint_optimal_vb(ep_jf, misrep_jks;
                                               mu_b=fi_sf_mu_b,
                                               rerun_fi_vb=true)

        # Add Objective Function columns
        ep_misrep_eqdf[:obj_fun] = String(sf_obj_fun)
        ep_misrep_eqdf[2, :obj_fun] = "misrep"
        ep_misrep_eqdf[:eq_type] = "misrep"
        
        # Reshape
        ep_misrep_eqdf = reshape_sf_rf_df(ep_misrep_eqdf)
    else
        ep_misrep_eqdf = DataFrame()
    end
    # ##################################################################

    eq_type = "all"
    run_pool_sep_eq = true
    if .&(run_pool_eq, !run_sep_eq)
        eq_type = "pooling"
    elseif .&(!run_pool_eq, run_sep_eq)
        eq_type = "separating"
    elseif .&(!run_pool_eq, !run_sep_eq)
        run_pool_sep_eq = false
    end
    ep_pool_eqdf = DataFrame()
    ep_sep_eqdf = DataFrame()
    # Electronic Platform Pooling and Separating Equilibria ############
    ep_eqdf = DataFrame()
    if run_pool_sep_eq
        ep_eqdf = ep_pool_sep_eq(ep_jf, deepcopy(ep_jks),
                                 fi_sf_mu_b,
                                 fi_rf_mu_b,
                                 fi_rf_obj_val;
                                 # ep_ks_path;
                                 equilibrium_type=eq_type,
                                 sf_obj_fun=sf_obj_fun,
                                 rf_obj_fun=rf_obj_fun,
                                 rerun=true,
                                 lb=lb, ub=ub,
                                 mu_bN=mu_bN, mu_bN2=mu_bN2,
                                 spline_k=spline_k,
                                 spline_bc=spline_bc)

        if ("pooling" in ep_eqdf[:eq_type])
            ep_pool_eqdf = ep_eqdf[ep_eqdf[:eq_type] .== "pooling", :]
        end
        if ("separating" in ep_eqdf[:eq_type])
            ep_sep_eqdf = ep_eqdf[ep_eqdf[:eq_type] .== "separating", :]
        end
    end
    # #################################################################


   # Form Electronic Platform Struct #################################
   return EPStruct(ep_jf, ep_sf_eqdf, ep_rf_eqdf,
                   ep_misrep_eqdf, ep_pool_eqdf, ep_sep_eqdf)
end






# TRASH ##############################################################
    # if .&(exists_ep_df(ep_ks_path, fidf_name), !rerun_full_info)
        
        # fidf = CSV.read(string(ep_ks_path, "/", fidf_name, ".csv");
        #                 types=fidf_col_types)


        # ep_sf_eqdf = slice_full_info_dataframe(fidf, ep_jks,
        #                                        ep_sf_svm;
        #                                        rerun=rerun_full_info)
        # ep_rf_eqdf = slice_full_info_dataframe(fidf, ep_jks,
        #                                        ep_rf_svm;
        #                                        rerun=rerun_full_info)

    # # Electronic Platform Separating Equilibrium ######################
    # # Separting Optimal mu_b
    # if run_sep_eq
    #     ep_sep_eqdf = ep_pool_sep_eq(ep_jf, deepcopy(ep_jks),
    #                                  fi_sf_mu_b,
    #                                  fi_rf_mu_b,
    #                                  fi_rf_obj_val,
    #                                  ep_ks_path;
    #                                  equilibrium_type="separating",
    #                                  sf_obj_fun=sf_obj_fun,
    #                                  rf_obj_fun=rf_obj_fun,
    #                                  rerun=rerun_misrep,
    #                                  lb=lb, ub=ub,
    #                                  mu_bN=mu_bN, mu_bN2=mu_bN2,
    #                                  spline_k=spline_k,
    #                                  spline_bc=spline_bc)
    # else
    #     ep_sep_eqdf = DataFrame()
    # end

    # sep_jks = deepcopy(ep_jks)
    # if .&(exists_ep_df(ep_ks_path, sepdf_name), !rerun_sep)
    #     sepdf = CSV.read(string(ep_ks_path, "/", sepdf_name, ".csv");
    #                      types=mps_col_types)

    #     # Slice DataFrame
    #     ep_sep_eqdf = slice_mps_dataframe(sepdf, sep_jks,
    #                                       ep_jf; rerun=rerun_pool)
    # else
    #     sepdf = find_joint_optimal_bond_measure(ep_jf, sep_jks,
    #                                             fi_sf_mu_b,
    #                                             fi_rf_mu_b,
    #                                             fi_rf_obj_val;
    #                                             equilibrium_type="separating",
    #                                             sf_obj_fun=sf_obj_fun,
    #                                             rf_obj_fun=rf_obj_fun,
    #                                             lb=lb, ub=ub,
    #                                             mu_bN=mu_bN, mu_bN2=mu_bN2,
    #                                             spline_k=spline_k,
    #                                             spline_bc=spline_bc)
        
    #     # Reshape
    #     ep_sep_eqdf = reshape_sf_rf_df(sepdf)
    # end
    # #################################################################


# ####################################################################
