
function joint_vb_extract_results(jf, jks,
                                   vbl::Float64;
                                   fi_sf_vb::Float64=NaN,
                                   fi_rf_vb::Float64=NaN,
                                   sf_defaults_first::Bool=true)

    if isnan(fi_sf_vb)
        fi_sf_vb = find_full_info_vb(jf.sf, jks)
    end
    if isnan(fi_rf_vb)
        fi_rf_vb = find_full_info_vb(jf.rf, jks)
    end

    sf_vb, rf_vb = get_type_contingent_vbs(vbl, fi_sf_vb, fi_rf_vb;
                                           sf_defaults_first=sf_defaults_first)
    
    # Run Equity Finite Differences Method
    df = joint_eq_fd(jf; jks=jks,
                     fi_sf_vb=fi_sf_vb,
                     fi_rf_vb=fi_rf_vb,
                     sf_vb=sf_vb,
                     rf_vb=rf_vb)
    

    # Store Results
    res = Dict{Symbol, Float64}(:fi_sf_vb => df[1, :fi_vb],
                                :fi_rf_vb => df[2, :fi_vb],
                                :sf_vb => df[1, :sf_vb],
                                :rf_vb => df[2, :rf_vb],
                                :vbl => df[1, :vb])
    for x in [:eq_deriv, :eq_deriv_min_val, :eq_min_val]
        res[Symbol(:sf_, x)] = df[1, x]
        res[Symbol(:rf_, x)] = df[2, x]
    end

    return res
end

function compute_joint_eq_vb_results(jf, jks;
                                     rerun_fi_vb::Bool=false,
                                     fi_sf_vb::Float64=NaN,
                                     fi_rf_vb::Float64=NaN,
                                     lb1::Float64=.75,
                                     ub1::Float64=1.25,
                                     vbN1::Int64=20,
                                     vbl_min::Float64=NaN,
                                     vbl_max::Float64=NaN,
                                     vbN2::Int64=20,
                                     lb2::Float64=.9, ub2::Float64=1.1)

    # Full Information: fi_sf_vb, fi_rf_vb
    jks = set_full_information_vb!(jf, jks;
                                   rerun_fi_vb=rerun_fi_vb,
                                   fi_sf_vb=fi_sf_vb,
                                   fi_rf_vb=fi_rf_vb,
                                   lb=lb1, ub=ub1,
                                   vbN=vbN1)

    # Form Grid of VB candidates
    vbl_min = minimum(x -> isnan(x) ? Inf : x, [jks.fi_sf_vb, jks.fi_rf_vb, vbl_min])
    vbl_max = maximum(x -> isnan(x) ? -Inf : x, [jks.fi_sf_vb, jks.fi_rf_vb, vbl_max])
    vbl_grid = range(lb2 * vbl_min, stop= ub2 * vbl_max, length=vbN2)



    # ####################################################################
    # Compute Joint Equity Finite Differences Method
    sf_res = fetch(@spawn [joint_vb_extract_results(jf, jks, vbl;
                                                    fi_sf_vb=jks.fi_sf_vb,
                                                    fi_rf_vb=jks.fi_rf_vb,
                                                    sf_defaults_first=true)
                           for vbl in vbl_grid])
    # Collect Results
    sf_resdf = vcat([DataFrame(x) for x in sf_res]...)

    # Risky Firm
    rf_res = fetch(@spawn [joint_vb_extract_results(jf, jks, vbl;
                                                    fi_sf_vb=jks.fi_sf_vb,
                                                    fi_rf_vb=jks.fi_rf_vb,
                                                    sf_defaults_first=false)
                           for vbl in vbl_grid])
    # Collect Results
    rf_resdf = vcat([DataFrame(x) for x in rf_res]...)


    return sf_resdf, rf_resdf
    # ####################################################################

    
    # # Compute Joint Equity Finite Differences Method
    # res = fetch(@spawn [joint_vb_extract_results(jf, jks, vbl;
    #                                              fi_sf_vb=jks.fi_sf_vb,
    #                                              fi_rf_vb=jks.fi_rf_vb)
    #                     for vbl in vbl_grid])
    # # Collect Results
    # resdf = vcat([DataFrame(x) for x in res]...)

    # # Extract Unique VB Values
    # resdf[:uniq_vb] = Array(resdf[:vbl])
    # cond = resdf[:vbl] .== resdf[:fi_sf_vb]
    # resdf[cond, :uniq_vb] = resdf[cond, :rf_vb]
    
    # # Sort Columns
    # cols = [:fi_rf_vb, :fi_sf_vb,
    #         :sf_vb, :rf_vb,
    #         :vbl, :uniq_vb,
    #         :rf_eq_deriv, :sf_eq_deriv,     
    #         :rf_eq_deriv_min_val, :sf_eq_deriv_min_val,
    #         :rf_eq_min_val, :sf_eq_min_val]
    
    # return resdf[cols]
end


function compile_opt_vb_results(jf, jks, sfdf::DataFrame, rfdf::DataFrame)

    # Case 1: Safe Firm defaults first! ############################
    opt_sf_vb, opt_rf_vb = interp_optimal_vbs(jf, jks, sfdf)

    # vbl that sets E'_s(vbl) to zero
    sf_vb, rf_vb = get_type_contingent_vbs(opt_sf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=true)
    s1 = joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)


    # vbl that sets E'_r(vbl) to zero
    sf_vb, rf_vb = get_type_contingent_vbs(opt_rf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=true)
    s2 = JointEq.joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
    
    # sfdf = vcat([s1, s2]...)
    # ###############################################################

    # Case 2: Risky Firm defaults first! ############################
    opt_sf_vb, opt_rf_vb = interp_optimal_vbs(jf, jks, rfdf)

    # vbl that sets E'_s(vbl) to zero
    sf_vb, rf_vb = get_type_contingent_vbs(opt_sf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=false)   
    r1 = JointEq.joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)

    # vbl that sets E'_r(vbl) to zero   
    sf_vb, rf_vb = get_type_contingent_vbs(opt_rf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=false)
    r2 = joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
    if abs.(r2[isnan.(r2[:sf_vb]), :eq_deriv][1]) > 1e-2
        r2 = refine_contingent_vbs(jf, jks, sf_vb, rf_vb)
        println(r2)
    end
    # ###############################################################
    println("=====================================================")
    println("=====================================================")
    println(r2)
    println("=====================================================")
    println("=====================================================")
    # Compile Results
    println("compiling results")
    df =  vcat([s1, s2, r1, r2]...)
    df[:sf_defaults_first] = vcat([fill(true, 4), fill(false, 4)]...)

    cols1 = [:sf_defaults_first,
             :fi_vb, :sf_vb, :rf_vb, :vb,
             :eq_deriv, :eq_min_val]
    return df[vcat(cols1, [x for x in names(df) if !(x in cols1)]...)]
end


function filter_joint_vb_results(df::DataFrame;
                                 tol1::Float64=-1e-2, tol2::Float64=1e-2)

    # Limited Liability Conditions ############################
    cond1 = df -> .&([df[x] .>= tol1 for x in [:eq_deriv, :eq_min_val]]...)

    # When the Safe Firm defaults first, 
    # its equity and equity derivative should be zero
    # at the joint default barrier
    sf_cond = df -> .&(df[:sf_defaults_first], 
                       isnan.(df[:sf_vb]) .| (df[:eq_deriv] .<= tol2))

    # When the Risky Firm defaults first, 
    # its equity and equity derivative should be zero
    # at the joint default barrier
    rf_cond = df -> .&(df[:sf_defaults_first] .==false, 
                       isnan.(df[:rf_vb]) .| (df[:eq_deriv] .<= tol2))

    # Limited Liability Conditions must be satisfied by both firms:
    llcond = df -> sum(.&(cond1(df), (sf_cond(df) .| rf_cond(df)))) == 2
    # #########################################################

    ## Find vb candidates
    vbdf = by(df, :vb, llcond = names(df) => x -> llcond(x))
    vblist = vbdf[vbdf[:llcond], :vb]

    # Filter DataFrame
    if !isempty(vblist)
        return df[in.(df[:vb], vblist), :]
    else
        println(string("No results for mu_b in ", unique(df[:mu_b])))
    end
end    


function find_joint_optimal_vb(jf, jks;
                               mu_s::Float64=NaN,
                               mu_b::Float64=NaN,
                               rerun_fi_vb::Bool=true,
                               fi_sf_vb::Float64=NaN,
                               fi_rf_vb::Float64=NaN,
                               lb1::Float64=.75,
                               ub1::Float64=1.25,
                               vbN1::Int64=20,
                               vbl_min::Float64=NaN,
                               vbl_max::Float64=NaN,
                               vbN2::Int64=20,
                               lb2::Float64=.9, ub2::Float64=1.1,
                               vbN3::Int64=10^5,
                               k_spline::Int64=3, bc_spline="extrapolate")

    # Measure of Safe Firms
    if !isnan(mu_s)
        setfield!(jks, :mu_s, mu_s)
    end

    # Measure of Bonds
    if !isnan(mu_b)
        setfield!(jks, :mu_b, mu_b)
    end

    # Set Full Information VBs
    if any([isnan(jks.fi_sf_vb),
            isnan(jks.fi_rf_vb),
            rerun_fi_vb])

        # Full Information: fi_sf_vb, fi_rf_vb
        jks = set_full_information_vb!(jf, jks;
                                       rerun_fi_vb=rerun_fi_vb,
                                       fi_sf_vb=fi_sf_vb,
                                       fi_rf_vb=fi_rf_vb,
                                       lb=lb1, ub=ub1,
                                       vbN=vbN1)
    end

    # Compute Equity Finite Differences Method
    sfdf, rfdf = compute_joint_eq_vb_results(jf, jks;
                                             rerun_fi_vb=rerun_fi_vb,
                                             fi_sf_vb=jks.fi_sf_vb,
                                             fi_rf_vb=jks.fi_rf_vb,
                                             lb1=lb1, ub1=ub1,
                                             vbN1=vbN1,
                                             vbl_min=vbl_min,
                                             vbl_max=vbl_max,
                                             vbN2=vbN2,
                                             lb2=lb2, ub2=ub2)

    # Get Candidates
    df = compile_opt_vb_results(jf, jks, sfdf, rfdf)

    # Filter to Back out optimal vb
    return filter_joint_vb_results(df)
end



