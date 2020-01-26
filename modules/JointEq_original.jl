module_path = "/home/artur/BondPricing/Julia/modules/"
push!(LOAD_PATH, module_path)
# modnames = ["ModelObj", "AnalyticFunctions",
#             "BondPrInterp", "EqFinDiff",
#             "Batch", "FullInfoEq"]
# for modl in modnames
#     if !(joinpath(module_path, modl) in LOAD_PATH)
#         push!(LOAD_PATH, joinpath(module_path, modl))
#     end
# end

include("FullInfoEq.jl")

module JointEq 

using Distributed
using Dierckx
using Parameters
using Printf
using DataFrames
using CSV
using Dates

using ModelObj: set_opt_k_struct,
                Firm, grid_creator,
                get_obj_model,
                set_opt_k_struct

using AnalyticFunctions: get_cvm_vb,
                         get_param,
                         rfbond_price,
                         on_default_payoff

using BondPrInterp: get_pv_rfdebt,
                    get_cvm_bond_price,
                    get_cvm_debt_price,
                    get_svm_bond_price,
                    get_svm_debt_price
    
using EqFinDiff: get_cvm_eq,
                 get_eq_Vmax,
                 eq_set_k_struct,
                 eq_fd_core,
                 eq_fd_export_results,
                 eq_fd

using Batch: BatchObj,
             get_bt,
             get_bt_svm,
             get_bt_mobj,
             get_batch_comb_num,
             load_cvm_opt_results_df,
             load_svm_opt_results_df,
             opt_k_struct_df_name,
             opt_k_struct_df_coltypes,
             BatchStruct,
             interp_values,
             svm_param_values_dict,
             common_params,
             _params_order,
             comb_folder_dict,
             str_format_fun,
             set_par_dict,
             comb_folder_dict,
             form_main_dir_path,
             main_dir, res_dir


using FullInfoEq: set_full_information_vb!,
                  find_optimal_bond_measure,
                  find_full_info_vb


# * Structs, Inputs Objects and Constructor Methods ############
# ** Joint Structs
# include("_joint_objects/_joint_structs.jl")
mutable struct JointKStruct
    mu_s::Float64
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64

    # Safe Firm
    fi_sf_vb::Float64
    sf_vb::Float64

    # Risk Firm
    fi_rf_vb::Float64
    rf_vb::Float64
end


@with_kw mutable struct FirmObj
    bt
    svm
end

mutable struct JointFirms
    jks #::JointKStruct
    sf::Firm
    rf::Firm
    cvm_bt
    svm_bt
    cvmdf::DataFrame
    svmdf::DataFrame
end


mutable struct JointFDParams
    # Safe Firm
    sf_eq_vbl::Float64
    sf_eq_max::Float64

    # Risk Firm
    rf_eq_vbl::Float64
    rf_eq_max::Float64
end



# #############################################################
mutable struct FirmSpecificParams
    iota::Float64
    lambda::Float64
    sigmal::Float64
    sigmah::Float64
end


mutable struct FirmCommonParams
    V0::Float64
    alpha::Float64
    pi::Float64
    r::Float64
    gross_delta::Float64
    xi::Float64
    # sigmal::Float64
end


mutable struct JointEqParams
    # Measure of Bonds
    mu_s::Float64
    
    # Transaction Costs
    kep::Float64
    kotc::Float64

    # Safe Firm Params
    sfp::FirmSpecificParams

    # Risky Firm Params
    rfp::FirmSpecificParams

    # Firm Common Params
    fcp::FirmCommonParams
end


mutable struct EPStruct
    # ep_ks::JointKStruct
    # sf::Firm
    # rf::Firm
    jf::JointFirms
    sfdf::DataFrame
    rfdf::DataFrame
    misrep::DataFrame
    pool::DataFrame
    sep::DataFrame
end


mutable struct OTCStruct
    sf::Firm
    rf::Firm
    sfdf::DataFrame
    rfdf::DataFrame
end


mutable struct JointEquilibrium
    # Joint Equilibrium Params
    jep::JointEqParams

    # Electronic Platform Data
    ep::EPStruct

    # OTC Markets Data
    otc::OTCStruct
end


# ** Joint Inputs
# Inputs -> in this order (need to define structs first)
# include("_joint_inputs.jl")
empty_jep = JointEqParams(vcat(fill(NaN,3), 
                               FirmSpecificParams(fill(NaN,4)...),
                               FirmSpecificParams(fill(NaN,4)...),
                               FirmCommonParams(fill(NaN, 6)...))...)

# Directories and File Names
ep_dir = "EP"
fi_dir = "FI"
jeq_dir = "JEQ"
fidf_name = "fidf"
misrepdf_name = "misrepdf"
pooldf_name = "pooldf"
sepdf_name = "sepdf"


# DataFrames Column Types
fidf_col_types = vcat(fill(Float64, 9), [Bool], fill(Float64, 17))
# Misrepresentation, Pool and Separating DFs Column Types
mps_col_types = vcat(String, Bool, fill(Float64, 11), Bool,
                     fill(Float64, 10), String, fill(Float64, 5),
                     Bool, fill(Float64, 10), String,
                     fill(Float64, 7))

# Electronic Market Parameters
epmcols = [:kappa, :mu_s, :m, :c, :p]
epm_eq_cols = [:obj_fun, :sf_defaults_first]

# Common Parameters
commoncols = vcat(fieldnames(FirmCommonParams)...)

# Firm Specific Parameters
fspcols = vcat(fieldnames(FirmSpecificParams)...)

# Default Boundaries
vbcols = [:fi_vb, :sf_vb, :rf_vb, :vb]


# Results Directories and File Names ################################ 
jeq_comb_folder_dict = deepcopy(comb_folder_dict)
jeq_comb_folder_dict[:kappa] = vcat(["kappa_", comb_folder_dict[:kappa][2:3]...])
jeq_comb_folder_dict[:alpha] = ["__alpha_", "%.2f"]
jeq_comb_folder_dict[:pi] = ["__pi_", "%.2f"]
jeq_comb_folder_dict[:r] = ["__r_", "%.3f"]

# Capital Structure
jeq_comb_folder_dict[:c] = ["__c_", comb_folder_dict[:c][2]]
jeq_comb_folder_dict[:p] = ["__p_", comb_folder_dict[:c][2]]
jeq_comb_folder_dict[:mu_s] = ["_mus_", "%.4f"]

# P/C ratio:
jeq_comb_folder_dict[:pcr] = ["_pcr_", "%.2f"]

common_dir_par_list = [:kappa, :sigmal] #, :gross_delta, :r, :alpha, :pi, :xi]
file_name_par_list = [:iota, :lambda, :sigmah]

eq_type_dict = Dict{String, Dict{Symbol, String}}(
    "full_info" => Dict{Symbol, String}(:dfn => fidf_name,
                                        :fn_prefix => "fi"),
    "misrep" => Dict{Symbol, String}(:dfn => misrepdf_name,
                                       :fn_prefix => "misrep"),
    "pooling" => Dict{Symbol, String}(:dfn => pooldf_name,
                                     :fn_prefix => "pool"),
    "separating" => Dict{Symbol, String}(:dfn => sepdf_name,
                                        :fn_prefix => "sep"))

dup_rows_params = [:eq_type, :kappa, :mu_s, 
                   :m, :c, :p,
                   :s_iota, :s_lambda, :s_sigmah,         
                   :s_delta, :s_obj_fun, #:s_eq_type,
                   :r_iota, :r_lambda, :r_sigmah,         
                   :r_delta, :r_obj_fun, #:r_eq_type,
                   :V0, :alpha, :pi, :r,
                   :gross_delta, :xi, :sigmal]
# ###################################################################


# ** Joint Ep Constructor
# include("_joint_objects/_joint_ep_constructor.jl")
function ep_pool_sep_eq(ep_jf, ep_jks,
                        fi_sf_mu_b::Float64,
                        fi_rf_mu_b::Float64,
                        fi_rf_obj_val::Float64;
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

    if !(equilibrium_type in ["all", "pooling", "separating"])
        println("Please set equilibrium_type to pooling or separting. Exiting...")
        return
    end
    
    pmdf = find_joint_optimal_bond_measure(ep_jf, deepcopy(ep_jks),
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
    
    return pmdf 
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
    for var in [:alpha, :pi, :r, :gross_delta, :xi] #, :sigmal]
        sf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
        rf_bt.mi._svm_dict[var] = getfield(jep.fcp, var)
    end
    sf_bt.mi._svm_dict[:sigmal] = jep.sfp.sigmal
    rf_bt.mi._svm_dict[:sigmal] = jep.rfp.sigmal   

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


    # Electronic Platform Full-Information Equilibrium ################
    if .&(isfile(fi_fpath_name), !rerun_full_info)
        fidf = CSV.read(fi_fpath_name)
                        #types=fidf_col_types)
        
        # Slice DataFrame
        ep_sf_eqdf = DataFrame(fidf[1, :]) 
        ep_rf_eqdf = DataFrame(fidf[2, :]) 

        # Capture Full Information Optimal mu_b and MBR ###################
        fi_sf_mu_b = ep_sf_eqdf[1, :mu_b]
        fi_rf_mu_b = ep_rf_eqdf[1, :mu_b]
        fi_sf_obj_val = ep_sf_eqdf[1, sf_obj_fun]
        fi_rf_obj_val = ep_rf_eqdf[1, rf_obj_fun]
        # #################################################################
    elseif rerun_full_info
        println(ep_sf_svm)
        println(ep_jks)

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
    if run_misrep
        misrep_jks = deepcopy(ep_jks)
        setfield!(misrep_jks, :mu_s, 1.)
        if !isnan(fi_sf_mu_b)
            setfield!(misrep_jks, :mu_b, fi_sf_mu_b)
        end
        
        ep_misrep_eqdf = find_joint_optimal_vb(ep_jf, misrep_jks;
                                               mu_b=fi_sf_mu_b,
                                               rerun_fi_vb=true)

        # Add Objective Function columns
        ep_misrep_eqdf[!, :obj_fun] .= String(sf_obj_fun)
        ep_misrep_eqdf[2, :obj_fun] = "misrep"
        ep_misrep_eqdf[!, :eq_type] .= "misrep"
        
        # Reshape
        ep_misrep_eqdf = reshape_sf_rf_df(ep_misrep_eqdf)
    else
        ep_misrep_eqdf = DataFrame()
    end
    # ##################################################################


    # Electronic Platform Pooling and Separating Equilibria ############
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
    if run_pool_sep_eq
        ep_eqdf = ep_pool_sep_eq(ep_jf, deepcopy(ep_jks),
                                 fi_sf_mu_b,
                                 fi_rf_mu_b,
                                 fi_rf_obj_val;
                                 equilibrium_type=eq_type,
                                 sf_obj_fun=sf_obj_fun,
                                 rf_obj_fun=rf_obj_fun,
                                 rerun=true,
                                 lb=lb, ub=ub,
                                 mu_bN=mu_bN, mu_bN2=mu_bN2,
                                 spline_k=spline_k,
                                 spline_bc=spline_bc)

        if ("pooling" in ep_eqdf[:, :eq_type])
            ep_pool_eqdf = ep_eqdf[ep_eqdf[:, :eq_type] .== "pooling", :]
        end
        if ("separating" in ep_eqdf[:, :eq_type])
            ep_sep_eqdf = ep_eqdf[ep_eqdf[:, :eq_type] .== "separating", :]
        end
    end
    # #################################################################


   # Form Electronic Platform Struct #################################
   return EPStruct(ep_jf, ep_sf_eqdf, ep_rf_eqdf,
                   ep_misrep_eqdf, ep_pool_eqdf, ep_sep_eqdf)
end


# ** Joint Constructor
# include("_joint_objects/_joint_constructors.jl")
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


# ** Joint Capital Structure Functions
# include("_joint_objects/_joint_k_struct_funs.jl")
function get_joint_k_struct!(jf;
                             jks=JointKStruct(fill(NaN, 10)...),
                             mu_b::Float64=NaN,
                             m::Float64=NaN,
                             c::Float64=NaN,
                             p::Float64=NaN)

    # if !isnan(mu_s)
    #     jks.mu_s = mu_s
    # else isnan(jks.mu_s)
    #     jks.mu_s = jf.jks.mu_s
    # end
    
    if !isnan(mu_b)
        setfield!(jks, :mu_b, mu_b)
    elseif isnan(jks.mu_b)
        setfield!(jks, :mu_b, jf.jks.mu_b)
    end

    if !isnan(m)
        setfield!(jks, :m, m)
    elseif isnan(jks.m)
        setfield!(jks, :m, jf.jks.m)
    end
    
    if !isnan(c)
        setfield!(jks, :c, c)
    elseif isnan(jks.c)
        setfield!(jks, :c, jf.jks.c)
    end

    if !isnan(p)
        setfield!(jks, :p, p)
    elseif isnan(jks.p)
        setfield!(jks, :p, jf.jks.p)
    end
   
    return jks 
end


function joint_eq_set_k_struct!(jf;
                                jks=JointKStruct(fill(NaN, 10)...),
                                mu_s::Float64=NaN,
                                mu_b::Float64=NaN,
                                m::Float64=NaN,
                                c::Float64=NaN,
                                p::Float64=NaN,
                                rerun_fi_vb::Bool=false,
                                fi_sf_vb::Float64=NaN,
                                sf_vb::Float64=NaN,
                                fi_rf_vb::Float64=NaN,
                                rf_vb::Float64=NaN,
                                lb::Float64=.75,
                                ub::Float64=1.25,
                                vbN::Int64=20)
    
    jks = get_joint_k_struct!(jf; jks=jks,
                              mu_b=mu_b,
                              m=m, c=c, p=p)

    # Default Barriers ##############################
    # Full Information: fi_sf_vb, fi_rf_vb
    jks = set_full_information_vb!(jf, jks;
                                   rerun_fi_vb=rerun_fi_vb,
                                   fi_sf_vb=fi_sf_vb,
                                   fi_rf_vb=fi_rf_vb,
                                   lb=lb, ub=ub,
                                   vbN=vbN)
    
    # setfield!(jks, :sf_vb, maximum(x->isnan(x) ? -Inf : x, [sf_vb, jks.fi_sf_vb]))

    # rf_vb = maximum([minimum(x->isnan(x) ? Inf : x, [rf_vb, jks.fi_rf_vb]),
    #                  jks.sf_vb])
    # setfield!(jks, :rf_vb, minimum(x->isnan(x) ? Inf : x, [rf_vb, jks.fi_rf_vb]))

    jks.sf_vb = !isnan(sf_vb) ? sf_vb : jks.fi_sf_vb
    jks.rf_vb = !isnan(rf_vb) ? rf_vb : jks.fi_rf_vb
    # Joint Equilibrium Barrier
    setfield!(jks, :vbl, maximum([jks.sf_vb, jks.rf_vb]))
    # ###############################################

    # Measure of Safe Firms
    if mu_s < 0.
        println("Setting mu_s to zero!")
        setfield!(jks, :mu_s, 0.0)
    elseif mu_s > 1.
        println("Setting mu_s to 1.!")
        setfield!(jks, :mu_s, 1.)
    elseif !isnan(mu_s)
        setfield!(jks, :mu_s, mu_s)
    end
    
    return jks
end


function get_type_contingent_vbs(vbl::Float64, 
                                 fi_sf_vb::Float64, 
                                 fi_rf_vb::Float64;
                                 sf_defaults_first::Bool=true)
    
    sf_vb = fi_sf_vb
    rf_vb = fi_rf_vb
    if vbl > maximum([fi_sf_vb, fi_rf_vb])
        sf_vb = sf_defaults_first ? vbl : fi_sf_vb
        rf_vb = !sf_defaults_first ? vbl : fi_rf_vb
    elseif  vbl < minimum([fi_sf_vb, fi_rf_vb])
        sf_vb = rf_vb = vbl
    elseif fi_sf_vb <= fi_rf_vb
        sf_vb = sf_defaults_first ? vbl : fi_sf_vb
        rf_vb = vbl
    else # fi_sf_vb > vbl > fi_rf_vb
        sf_vb = vbl
        rf_vb = !sf_defaults_first ? vbl : fi_rf_vb
    end
     
    
    return sf_vb, rf_vb
end


function interp_optimal_vbs(jf, jks, df::DataFrame; 
                            sf_defaults_first::Bool=true, 
                            vbN::Int64=10^5,
                            spline_k::Int64=3,
                            spline_bc::String="extrapolate")

    # Form Refined Grid of Unique vbl values
    vb_grid = range(minimum(df[:, :vbl]), stop=maximum(df[:, :vbl]), length=vbN)
    
    # Interpolate Equity and Equity Derivative Functions
    tmp = Dict()
    for var in [:rf_eq_deriv, :sf_eq_deriv] 
        tmp[var] = Dierckx.Spline1D(df[:, :vbl], df[:, var]; k=spline_k, bc=spline_bc)
    end

    # Find Optimal VBs
    opt_sf_vb = vb_grid[argmin(abs.(tmp[:sf_eq_deriv](vb_grid)))]
    opt_rf_vb = vb_grid[argmin(abs.(tmp[:rf_eq_deriv](vb_grid)))]

    return opt_sf_vb, opt_rf_vb
end


function refine_contingent_vbs(jf, jks, 
                               sf_vb::Float64, rf_vb::Float64; 
                               sf_defaults_first::Bool=false,
                               N::Int64=7,
                               spline_k::Int64=3, 
                               spline_bc::String="extrapolate")
    
    tvb = !(sf_defaults_first) ? rf_vb : sf_vb
    tvar = !(sf_defaults_first) ? :rf_vb : :sf_vb
    cond_var = !(sf_defaults_first) ? :sf_vb : :rf_vb
    
    dfL = []
    for vb in range(.95 * tvb, stop = 1.05 * tvb, length=N)
        sf_vb, rf_vb = JointEq.get_type_contingent_vbs(vb,
                                                       jks.fi_sf_vb,
                                                       jks.fi_rf_vb; 
                                                       sf_defaults_first=sf_defaults_first)
        tmp = JointEq.joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
        push!(dfL, tmp)
    end
    df = vcat(dfL...)
    cond = isnan.(df[cond_var])
    
    tvb_grid = range(minimum(df[cond, tvar]), stop=maximum(df[cond, tvar]), length=10^5)
    tvb_interp = Dierckx.Spline1D(df[cond, tvar], df[cond, :eq_deriv]; 
                                  k=spline_k, bc=spline_bc)
    
    
    opt_tvb = tvb_grid[argmin(abs.(tvb_interp(tvb_grid)))]
    sf_vb, rf_vb = get_type_contingent_vbs(opt_tvb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=sf_defaults_first)
    
    return joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
end
# ##########################################################


# * Auxiliary File and DataFrame Methods ########################
# ** Joint Auxiliary Functions
# include("_joint_auxiliary/_joint_functions.jl")
function store_joint_eq_parameters(mu_s::Float64,
                                   kep::Float64,
                                   kotc::Float64;
                                   jep=empty_jep,
                                   sfp=FirmSpecificParams(fill(NaN, 4)...),
                                   rfp=FirmSpecificParams(fill(NaN, 4)...),
                                   s_iota::Float64=NaN,
                                   s_lambda::Float64=NaN,
                                   s_sigmal::Float64=NaN,
                                   s_sigmah::Float64=NaN,
                                   r_iota::Float64=NaN,
                                   r_lambda::Float64=NaN,
                                   r_sigmal::Float64=NaN,
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
    if !isnan(s_sigmal)    
        setfield!(jep.sfp, :sigmal, s_sigmal)
    else
        setfield!(jep.sfp, :sigmal, svm_param_values_dict[:sigmal][1])
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
    if !isnan(r_sigmal)    
        setfield!(jep.rfp, :sigmal, r_sigmal)
    else
        setfield!(jep.rfp, :sigmal, svm_param_values_dict[:sigmal][1])
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
                           svm_param_values_dict[:xi][1]) #,
                           # svm_param_values_dict[:sigmal][1])
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


# ** Joint File Methods
# include("_joint_auxiliary/_joint_file_methods.jl")
function get_jeq_mus_fname(jks_fpath::String; eq_type::String="pooling", mu_s::Float64=NaN)
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end

    musn="_tmp"
    if !isnan(mu_s)
        musn = string(jeq_comb_folder_dict[:mu_s][1], 
                      str_format_fun(jeq_comb_folder_dict[:mu_s][2], mu_s))
    end
    return string(jks_fpath, "/", eq_type_dict[eq_type][:fn_prefix], musn, ".csv")
end


function get_jeq_contour_fname(jks_fpath::String, comb_num::Int64; eq_type::String="full_info")
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end

    return string(jks_fpath, "/", eq_type_dict[eq_type][:fn_prefix], "_", comb_num, ".csv")
end


function collect_joint_eq_files(jks_fpath::String; eq_type::String="full_info")
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end
    
    return [CSV.read(string(jks_fpath, "/", x)) for x in readdir(jks_fpath) 
                if occursin(string(eq_type_dict[eq_type][:fn_prefix], "_"), x)]
end


function process_results(jks_fpath::String, eq_type::String)
    if !(eq_type in keys(eq_type_dict))
        println("Please set equilibrium_type to 'full_info', 'misrep', 'pooling' or 'separating'. Exiting...")
        return
    end
    
    df_l = [CSV.read(string(jks_fpath, "/", x)) for x in readdir(jks_fpath) 
            if occursin(string(eq_type_dict[eq_type][:fn_prefix], "_"), x)]
    df = vcat(df_l...)
    if eq_type == "full_info"
        df = df[nonunique(df[:, [:iota, :sigmah]]) .== false, :]
        df = vcat(sort(df[isnan.(df[:, :sigmah]), :], [:iota]),
                  sort(df[df[:, :iota] .== .0, :], [:sigmah]))   
    else
        df = vcat(sort(df[isnan.(df[:, :r_sigmah]), :], [:r_iota]),
                  sort(df[df[:, :r_iota] .== .0, :], [:r_sigmah]))
    end            
end


# ** Joint Dataframe Methods
# include("_joint_auxiliary/_joint_dataframe_methods.jl")
# Match DataFrames' Rows ###########################################################
function identify_match(df::DataFrame, dfrow::DataFrame)
    compf = var -> abs.(df[:, var] .- dfrow[1, var]) .< 1e-6
    fspf = var -> isnan.(dfrow[1, var]) ? isnan.(df[:, var]) : compf(var)
    
    fcp_cond  = [compf(var) for var in commoncols]
    fsp_cond = [fspf(var) for var in fspcols]
    
    return .&([.&(x...) for x in [fcp_cond, fsp_cond]]...)
end

    
function identify_matching_row(df::DataFrame, svm)
    compf = var -> abs.(df[:, var] .- get_param(svm, var)) .< 1e-6
    fspf = var -> isnan.(get_param(svm, var)) ? isnan.(df[:, var]) : compf(var)
    
    fcp_cond  = [compf(var) for var in commoncols]
    fsp_cond = [fspf(var) for var in fspcols]
    
    return .&([.&(x...) for x in [fcp_cond, fsp_cond]]...)
end


function identify_match2(df::DataFrame, dfrow::DataFrame)
    compf = var -> abs.(df[:, var] .- dfrow[1, var]) .< 1e-6
    scompf = (var, prefix) -> (abs.(df[:, Symbol(prefix, var)] .- 
                               dfrow[1, Symbol(prefix, var)]) .< 1e-6)
    fspf = (var, prefix) -> (isnan.(dfrow[1, Symbol(prefix, var)]) ? 
                             isnan.(df[:, Symbol(prefix, var)]) : scompf(var, prefix))

    # Check Parameters
    fcp_cond  = [compf(var) for var in vcat(:kappa, commoncols)]
    sfsp_cond = [fspf(var, :s_) for var in fspcols]
    rfsp_cond = [fspf(var, :r_) for var in fspcols]
    epm_cond = [compf(var) for var in epmcols if var != :kappa]
    
    return .&([.&(x...) for x in [fcp_cond, sfsp_cond, rfsp_cond, epm_cond]]...)
end


function identify_matching_row2(df::DataFrame, jks, jf)
    compf = (svm, var) -> abs.(df[:, var] .- get_param(svm, var)) .< 1e-6
    scompf = (svm, var, prefix) -> abs.(df[:, Symbol(prefix, var)] .- get_param(svm, var)) .< 1e-6
    fspf = (svm, var, prefix) -> isnan.(get_param(svm, var)) ? isnan.(df[:, Symbol(prefix, var)]) : scompf(svm, var, prefix)

    # Parameter Comparison Function
    pcf = var -> abs.(df[:, var] .- getfield(jks, var)) .< 1e-6

    # Check Parameters
    fcp_cond  = [compf(jf.sf, var) for var in vcat(:kappa, commoncols)]
    sfsp_cond = [fspf(jf.sf, var, :s_) for var in fspcols]
    rfsp_cond = [fspf(jf.rf, var, :r_) for var in fspcols]
    epm_cond = [pcf(var) for var in epmcols if var != :kappa]
    
    return .&([.&(x...) for x in [fcp_cond, sfsp_cond, rfsp_cond, epm_cond]]...)
end
# ##################################################################################


# Reshape Misrep, Pool and Separating Eq DFs #######################################
function reshape_sf_rf_df(df::DataFrame)
    jeq_commoncols = commoncols
    
    # Case in which both firms are CVM
    if .&(isnan.(df[:, :sigmah])...)
        jeq_commoncols = [x for x in jeq_commoncols if x != :sigmal] 
    end
    
    # List of firm-specific variables
    non_specific_cols  = vcat(:sf_defaults_first,
                              :eq_type,
                              epmcols,
                              jeq_commoncols, vbcols)

    specific_cols = [col for col in names(df) if 
                     !(col in non_specific_cols)]
        

    # Capital Structure
    if df[1, :eq_type] == "full_info"
        ksdf = DataFrame(df[1, [:eq_type]])
        sfdf = DataFrame(df[1, specific_cols])
        rfdf = DataFrame(df[2, specific_cols])
    else
        ksdf = DataFrame(df[1, [:eq_type, :sf_defaults_first, epmcols..., :vb]])
        sfdf = DataFrame(df[1, vcat([:sf_vb], specific_cols)])
        rfdf = DataFrame(df[2, vcat([:rf_vb], specific_cols)])
    end

    # Common Parameters
    commondf =  DataFrame(df[1, jeq_commoncols])

    # Safe Firm Results
    names!(sfdf, vcat([:sf_vb], [Symbol(:s_, col) for col in specific_cols
                                 if (col != :sf_defaults_first)]))

    # Risky Firm Results
    names!(rfdf, vcat([:rf_vb], [Symbol(:r_, col) for col in specific_cols
                                 if (col != :sf_defaults_first)]))

    return hcat(ksdf, sfdf, rfdf, commondf)
end
# ###############################################################################


# Slice DataFrames ##############################################################
function slice_df_cond(df::DataFrame, svm, rerun::Bool)
    row_cond = false
    if !isempty(df)
        row_cond = identify_matching_row(df, svm)
    end

    if (sum(row_cond) == 0)
        println("No matches found!")
        rerun = true
    elseif .&(sum(row_cond) == 1, !rerun)
        println("Match found!")
    end
        
    return row_cond, rerun
end


function slice_full_info_dataframe(df::DataFrame, jks, svm; rerun::Bool=true)
    row_cond, rerun = slice_df_cond(df, svm, rerun)

    if rerun
        println("Generating results...")
        dfrow = find_optimal_bond_measure(svm; jks=jks) 
    else
        println("Extracting row...")
        dfrow = DataFrame(df[row_cond, :])
    end
        
    return dfrow
end


function slice_mps_dataframe(df::DataFrame, jks, jf; rerun::Bool=true)
    row_cond = false
    if !isempty(df)
        row_cond = identify_matching_row2(df, jks, jf)
    end

    if (sum(row_cond) == 0)
        println("No matches found!")
        rerun = true
    elseif .&(sum(row_cond) == 1, !rerun)
        println("Match found!")
    end

    if rerun
        println("Generating results...")
        dfrow = find_joint_optimal_vb(jf, jks;
                                      rerun_fi_vb=true)
    else
        println("Extracting row...")
        dfrow = DataFrame(df[row_cond, :])
    end
       
    return dfrow
end
# ###############################################################################


function joint_eq_form_dataframes(;pool_list::Array{DataFrame,1}=[DataFrame()],
                                  sep_list::Array{DataFrame,1}=[DataFrame()])
   
    pooldf = DataFrame()
    sepdf = DataFrame()


    dtime = Dates.now()
    if !.&(size(pool_list, 1) == 1, isempty(pool_list[1]))
        pooldf = sort!(DataFrame(vcat(pool_list...)), :mu_s)
        pooldf[:datetime] = dtime
    end
    if !.&(size(sep_list, 1) == 1, isempty(sep_list[1]))
        sepdf = sort!(DataFrame(vcat(sep_list...)), :mu_s)
        sepdf[:datetime] = dtime
    end

    return pooldf, sepdf
end
                                    

# Save DataFrames ###############################################################
function update_file(fpath::String, fname::String, df::DataFrame)
    if !JointEq.exists_ep_df(fpath, fname)
        CSV.write(string(fpath, "/", fname, ".csv"), df)
    else
        if fname == fidf_name
            resdf = CSV.read(string(fpath, "/", fname, ".csv"); types=fidf_col_types)
        else
            resdf = CSV.read(string(fpath, "/", fname, ".csv"); types=mps_col_types)
        end
        
        for i in 1:nrow(df)
            # check if there is a match
            if fname == fidf_name
                row_cond = identify_match(resdf, DataFrame(df[i, :]))
            else
                row_cond = identify_match2(resdf, DataFrame(df[i, :]))
            end
            
            if sum(row_cond) == 1
                println("Match found! Replacing row...")
                resdf[row_cond, :] = DataFrame(df[i, :])
            elseif sum(row_cond) == 0
                println("No match found. Appending row...")
                resdf = vcat(resdf, DataFrame(df[i, :]))
            else
                println(string("Multiple matches for row ", i, 
                               "! Please refine ID columns. Exiting..."))
                return
            end
        end
        
        # Save File
        println(string("Saving File ", fname," ..."))
        CSV.write(string(fpath, "/", fname, ".csv"), resdf)
    end
end


function remove_dup_save_df(df::DataFrame, eq_type::String, jks_fpath::String)
    if eq_type in keys(eq_type_dict)
        df_name = eq_type_dict[eq_type][:dfn]
    else
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end
    
    df_fpath_name = string(jks_fpath, "/", df_name, ".csv")
    cols = (eq_type == "full_info") ? [x for x in names(df) if x != :datetime] : dup_rows_params
    
    if (string(df_name, ".csv") in readdir(jks_fpath))
        df_all = CSV.read(df_fpath_name)

        # Add new row
        df_all = vcat(df_all, df)

        # Move new row to the top
        sort!(df_all, :datetime, rev=true)

        # Remove old duplicate row
        df_all = unique(df_all, cols)
    else
        df_all = df
    end
    # Save file
    CSV.write(df_fpath_name, df_all)
    
    return df_all
end    
# ###############################################################################


# ** Joint Set and Get Dir Methods
# include("_joint_auxiliary/_joint_set_get_dirs.jl")
# Get Directories ##############################################
function get_jeq_common_params_dir(jf)
    common_param_dir = ""
    for par in common_dir_par_list
        if par in [:kappa, :gross_delta]
            par_string =  string(str_format_fun(jeq_comb_folder_dict[par][2],
                                                1e4 * get_param(jf.sf, par)),
                                 jeq_comb_folder_dict[par][3])
        else
            par_string =  string(str_format_fun(jeq_comb_folder_dict[par][2],
                                                get_param(jf.sf, par)))
        end
        common_param_dir = string(common_param_dir, jeq_comb_folder_dict[par][1], par_string)
    end

    return common_param_dir
end


function get_jks_dir(jf)
    pcr = jf.jks.p/jf.jks.c
    
    return string(jeq_comb_folder_dict[:m][1], str_format_fun(jeq_comb_folder_dict[:m][2], jf.jks.m), 
                  jeq_comb_folder_dict[:pcr][1], str_format_fun(jeq_comb_folder_dict[:pcr][2], pcr))
end


function get_jeq_jks_fpath(jf)
    jeq_fpath = make_jeq_res_fpath()
    common_params_dir = get_jeq_common_params_dir(jf)
    
    common_params_fpath = string(jeq_fpath, "/", common_params_dir)
    jks_dir = get_jks_dir(jf)

    return string(common_params_fpath, "/", jks_dir)
end
# ##############################################################


# Get Paths ####################################################
function get_res_fpath(jf)
    return string(jf.svm_bt.mi.main_dir_path, "/", res_dir)
end


function get_jeq_res_fpath(jf)
    res_dir_path = make_res_fpath(jf)

    return string(res_dir_path, "/", jeq_dir)
end


function get_jeq_common_params_fpath(jf)
    jeq_fpath = make_jeq_res_fpath(jf)
    common_params_dir = get_jeq_common_params_dir(jf)
    
    return string(jeq_fpath, "/", common_params_dir)
end


function get_jeq_jks_fpath(jf)
    common_params_fpath =  get_jeq_common_params_fpath(jf)
    jks_dir = get_jks_dir(jf)

    return string(common_params_fpath, "/", jks_dir)
end
# ##############################################################


# Make Directories #############################################
function make_res_fpath(jf)
    res_dir_path = get_res_fpath(jf)
    if !isdir(res_dir_path)
        mkdir(res_dir_path)
    end

    return res_dir_path
end


function make_jeq_res_fpath(jf)
    jeq_fpath = get_jeq_res_fpath(jf)
    if !isdir(jeq_fpath)
        mkdir(jeq_fpath)
    end

    return jeq_fpath
end


function make_jeq_common_params_fpath(jf)
    common_params_fpath = get_jeq_common_params_fpath(jf)
    if !isdir(common_params_fpath)
        mkdir(common_params_fpath)
    end

    return common_params_fpath
    
end


function make_jeq_jks_fpath(jf)
    common_params_fpath = make_jeq_common_params_fpath(jf)
    jks_fpath = get_jeq_jks_fpath(jf)
    if !isdir(jks_fpath)
        mkdir(jks_fpath)
    end

    return jks_fpath 
end
# ##############################################################


# * Joint Pricing Functions #####################################
# include("_joint_pricing.jl")
function adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)
    # Capital Structure
    for fn in fieldnames(JointKStruct)
        if isnan(getfield(jks, fn))
            setfield!(jks, fn, getfield(jf.jks, fn))
        end
    end

    # Measure of Safe Firms
    if isnan(mu_s)
        mu_s = jks.mu_s
    elseif mu_s < 0.
        println("Setting mu_s to zero!")
        mu_s = 0.0
    elseif mu_s > 1.
        println("Setting mu_s to 1.!")
        mu_s = 1.
    end
    

    # Set Asset Value
    if .&(isnan(Vt), isnan(vt))
        Vt = get_param(jf.sf, :V0)
        vt = log(Vt/jks.vbl)
    elseif isnan(vt)
        vt = log(Vt/jks.vbl)
    end
    
    return jks, mu_s, vt
end


function joint_bond_price(jf, ttm::Float64;
                          jks=JointKStruct(fill(NaN, 10)...),
                          mu_s::Float64=NaN,
                          vt::Float64=NaN, Vt::Float64=NaN,
                          sf_ftype::String="bf",
                          rf_ftype::String="bf")

    jks, mu_s, vt = adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)

    # Compute Bond Prices
    if get_obj_model(jf.sf) == "cvm" 
        sf_bpr = get_cvm_bond_price(jf.sf, ttm, jf.sf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt),
                                    vb=jks.sf_vb,
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p)
    else
        sf_bpr = get_svm_bond_price(jf.sf, jks.sf_vb, ttm;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ftype=sf_ftype)
    end

    if get_obj_model(jf.rf) == "cvm" 
        rf_bpr = get_cvm_bond_price(jf.rf, ttm, jf.rf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt),
                                    vb=jks.rf_vb,
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p)
    else   
        rf_bpr = get_svm_bond_price(jf.rf, jks.rf_vb, ttm;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ftype=rf_ftype)
    end
 
    # Joint Price
    return mu_s * sf_bpr + (1. - mu_s) * rf_bpr
end


function joint_debt_price(jf;
                          jks=JointKStruct(fill(NaN, 10)...),
                          mu_s::Float64=NaN,
                          vt::Float64=NaN,
                          Vt::Float64=NaN,
                          ttmN0::Int64=10^2,
                          ttmN::Int64=10^4)

    jks, mu_s, vt = adjust_pricing_parameters(jf, jks, mu_s, Vt, vt)

    
    # Compute Debt Price #################################################
    # either model is CVM or bond surface maturity
    # must coincide with capital structure maturity.
    if get_obj_model(jf.sf) == "cvm" 
        sf_dpr = get_cvm_debt_price(jf.sf, jks.sf_vb,
                                    jf.sf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p,
                                    N1=ttmN0, N2=ttmN)
    elseif abs.(jf.sf.m - jks.m) < 1e-4
        sf_dpr = get_svm_debt_price(jf.sf, jks.sf_vb;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ttmN=ttmN)
    else
        println(string("Capital Structure Maturity of Safe Type ",
                       "does not coincide with ",
                       "m used in bond pricing surfaces computation!"))
        return NaN
    end

    if get_obj_model(jf.rf) == "cvm" 
        rf_dpr = get_cvm_debt_price(jf.rf, jks.rf_vb,
                                    jf.rf.pm.sigmal;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, m=jks.m,
                                    c=jks.c, p=jks.p,
                                    N1=ttmN0, N2=ttmN)
    elseif abs.(jf.rf.m - jks.m) < 1e-4
        rf_dpr = get_svm_debt_price(jf.rf, jks.rf_vb;
                                    Vt=jks.vbl * exp(vt), 
                                    mu_b=jks.mu_b, # m=jks.m,
                                    c=jks.c, p=jks.p, ttmN=ttmN)
    else
        println(string("Capital Structure Maturity of Risky Type ",
                       "does not coincide with ",
                       "m used in bond pricing surfaces computation!"))
        return NaN
    end

    # Joint Price
    return mu_s * sf_dpr + (1. - mu_s) * rf_dpr
end


function joint_eq_fd_newly_issued_bonds(jf, jks, vbl::Float64,
                                        vgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
                                                            Base.TwicePrecision{Float64}};
                                        vtN::Int64=10^3,
                                        sf_ftype::String="bf",
                                        rf_ftype::String="bf",
                                        spline_k::Int64=3,
                                        spline_bc::String="extrapolate")

    # Common Payoffs ###############################
    rfbond = rfbond_price(jks.m, jks.c, jks.p,
                          jf.sf.pm.r, jf.sf.pm.xi,
                          jf.sf.pm.kappa)
    
    #dpayoff = on_default_payoff(0., vbl, jks.m, jks.mu_b,
    #                            jks.m, jks.c, jks.p,
    #                            jf.sf.pm.r, jf.sf.pm.xi,
    #                            jf.sf.pm.kappa, jf.sf.pm.alpha)
    # ##############################################

    
    _, v_subgrid = grid_creator((1 + 1e-4) * minimum(jf.sf.bs.vtgrid), maximum(vgrid), vtN)

    # ##############################################
    # Whether to compute bond pricing surfaces
    # (fixed maturity) on the fly
    if any([abs.(jf.sf.m - jks.m) > 1e-4, sf_ftype == "bft"])
        sf_ftype = "bft"
        jf.sf = bpr_interp_fixed_ttm(jf.sf; ttm=jks.m)
    end

    if any([abs.(jf.rf.m - jks.m) > 1e-4, rf_ftype == "bft"])
        rf_ftype = "bft"
        jf.rf = bpr_interp_fixed_ttm(jf.rf; ttm=jks.m)
    end
    # ###############################################

    bpr_vec = fetch(@spawn [joint_bond_price(jf, jks.m; vt=v,
                                             jks=jks,
                                             mu_s=jks.mu_s,
                                             sf_ftype=sf_ftype,
                                             rf_ftype=rf_ftype) for v in v_subgrid])
    
    #bpr = Dierckx.Spline1D(vcat(.0, v_subgrid), vcat(dpayoff, bpr_vec); k=3, bc="extrapolate")
    bpr = Dierckx.Spline1D(v_subgrid, bpr_vec; k=spline_k, bc=spline_bc)

    return Array([minimum([bpr(v)[1], rfbond]) for v in vgrid])[2:end-1]
end


# * Equilibrium Methods -> Equity, Vb, mu_b Functions ###########
# ** Joint Equity Finite Differences
# include("_joint_equilibrium/_joint_eq_fin_diff.jl")
function joint_eq_get_Vmax_vgrid(jf, jks, vbl::Float64; vtN::Int64=1500)
    # V MAX ######################################
    println("Computing Equity Vmax")
    # Safe Firm
    sf_eq_Vmax = get_eq_Vmax(jf.sf; mu_b=jks.mu_b, m=jks.m, c=jks.c, p=jks.p)

    # Risky Firm
    rf_eq_Vmax = get_eq_Vmax(jf.rf; mu_b=jks.mu_b, m=jks.m, c=jks.c, p=jks.p)

    # Take the Maximum Value:
    eq_Vmax = maximum([sf_eq_Vmax, rf_eq_Vmax])
    println(string("Equity Vmax: ", eq_Vmax))
    println(" ")
    # ############################################
    
    # vtgrid
    vtgrid = reverse(range(0.0, stop=log(eq_Vmax/float(vbl)), length=vtN))
    
    # #################################
    # ######## Boundary Values ########
    # #################################
        
    println(string("eq_max: ", vtgrid[1]))
    
    # Lower Barrier:
    
    return eq_Vmax, vtgrid
end 

                   
function joint_eq_get_boundary_values(tj, jks,
                                      vbj::Float64, vbl::Float64, eq_Vmax::Float64)
    # Lower Barriers
    if abs(vbj - vbl) < 1e-4 #|| get_obj_model(tj) == "cvm"   
        tj_eq_vbl = maximum([0., get_param(tj, :alpha) * vbl -
                             get_pv_rfdebt(tj; mu_b=jks.mu_b,
                                           m=jks.m, c=jks.c, p=jks.p)])
    else
        eqvals = eq_fd(tj; vbl=vbj, mu_b=jks.mu_b, m=jks.m, c=jks.c, p=jks.p, V0=vbl)

        if .&(eqvals[1, :eq_min_val] > -1e-2,
              eqvals[1, :eq_deriv_min_val] > -1e-2)
            tj_eq_vbl = maximum([0.0, eqvals[1, :equity]])
        else
            tj_eq_vbl = 0.0
        end
    end
    println(string("eq_vbl: ", tj_eq_vbl))

    # Upper Barriers: Value of Equity
    tj_eq_max = get_cvm_eq(tj, eq_Vmax, tj.pm.sigmal;
                           mu_b=jks.mu_b, m=jks.m, c=jks.c, p=jks.p)

    return tj_eq_vbl, tj_eq_max 
end


# DATAFRAME COLUMNS
# K Structure
jks_cols = [:mu_s, :m, :mu_b, :c, :p]

# Default Barrier
vb_cols = [:fi_vb, :sf_vb, :rf_vb, :vb]

# EFD
share_cols = [:eq_deriv, :eq_deriv_min_val, 
              :eq_min_val, :eq_negative, 
              :eq_vb, :MBR, :debt, :equity, 
              :firm_value, :leverage]

# Parameters
param_cols = [:iota, :lambda, :sigmah, 
              :gross_delta, :delta, :kappa, 
              :sigmal, :V0, :xi, :r, :alpha, :pi]

jks_eq_fd_cols = vcat(jks_cols, vb_cols, share_cols, param_cols)


# Joint Equilibrium Equity Finite Differences
function joint_eq_fd(jf; V0::Float64=NaN,
                     jks=JointKStruct(fill(NaN, 10)...),
                     mu_s::Float64=NaN,
                     mu_b::Float64=NaN,
                     m::Float64=NaN,
                     c::Float64=NaN,
                     p::Float64=NaN,
                     vbl::Float64=NaN,
                     sf_vb::Float64=NaN,
                     rf_vb::Float64=NaN,
                     fi_sf_vb::Float64=NaN,
                     fi_rf_vb::Float64=NaN,
                     cov_gross_delta::Float64=NaN,
                     debt::Float64=NaN,
                     sf_ftype::String="bf",
                     rf_ftype::String="bf",
                     lb::Float64=.75, ub::Float64=1.25, vbN::Int64=15,
                     vtN::Int64=1500)

    tic = time()

    # Common Parameters #############################
    if !check_param_consistency(jf)
        println("Exiting...")
        return
    end
    # ###############################################

    # Set Capital Structure #########################
    jks = joint_eq_set_k_struct!(jf; jks=jks,
                                 mu_s=mu_s,
                                 mu_b=mu_b,
                                 m=m, c=c, p=p,
                                 fi_sf_vb=fi_sf_vb,
                                 sf_vb=sf_vb,
                                 fi_rf_vb=fi_rf_vb,
                                 rf_vb=rf_vb)
    # ###############################################

    println("================================================================")
    println("  ")
    if !isnan(vbl)
        setfield!(jks, :vbl, vbl)
    end
    println(string("vbl: ", getfield(jks, :vbl)))
    println("  ")
    println("================================================================")
    
    # Equity Boundary Values ########################
    eq_Vmax, vtgrid = joint_eq_get_Vmax_vgrid(jf, jks, jks.vbl; vtN=vtN)

    sf_eq_vbl, sf_eq_max = joint_eq_get_boundary_values(jf.sf, jks,
                                                        jks.sf_vb,
                                                        jks.vbl, eq_Vmax)
    rf_eq_vbl, rf_eq_max = joint_eq_get_boundary_values(jf.rf, jks,
                                                        jks.rf_vb,
                                                        jks.vbl, eq_Vmax)
    # ###############################################

    # Store Parameters ############################## 
    fdp = JointFDParams(sf_eq_vbl,
                        sf_eq_max,
                        rf_eq_vbl,
                        rf_eq_max)
     # ###############################################   

    # Newly-Issued Bond Prices ######################
    bond_prices = joint_eq_fd_newly_issued_bonds(jf, jks, jks.vbl,
                                                 vtgrid;
                                                 vtN=vtN,
                                                 sf_ftype=sf_ftype,
                                                 rf_ftype=rf_ftype)

    # ###############################################

    # Debt Price ###################################
    debt_pr = joint_debt_price(jf; jks=jks)

    # ##############################################
    
    # Compute Equity Values #########################
    # No adjustments to vtgrid, because
    # bond_prices and sf_eq_vbl already adjust for 
    # the differences in the default barrier.
    sf_eq_dict = eq_fd_core(jf.sf, jks, jks.vbl,
                            sf_eq_vbl, sf_eq_max,
                            vtgrid, bond_prices) # no adjustments
                            # vtgrid .-log(jks.vbl/jks.sf_vb), bond_prices)
    _, sf_df = eq_fd_export_results(jf.sf, jks, jks.vbl, sf_eq_dict; debt=debt_pr)

    sf_df[!, :mu_s] .= jks.mu_s
    sf_df[!, :fi_vb] .= jks.fi_sf_vb 
    sf_df[!, :sf_vb] .= jks.sf_vb
    sf_df[!, :rf_vb] .= NaN   
    
    rf_eq_dict = eq_fd_core(jf.rf, jks, jks.vbl,
                            rf_eq_vbl, rf_eq_max,
                            vtgrid, bond_prices;
                            cov_gross_delta=cov_gross_delta)
                            # vtgrid .-log(jks.vbl/jks.rf_vb), bond_prices)
    _, rf_df = eq_fd_export_results(jf.rf, jks, jks.vbl, rf_eq_dict; debt=debt_pr)
    rf_df[!, :mu_s] .= jks.mu_s
    rf_df[!, :fi_vb] .= jks.fi_rf_vb 
    rf_df[!, :sf_vb] .= NaN
    rf_df[!, :rf_vb] .= jks.rf_vb
    # ##############################################

    println(string("Total computation time: ", time() - tic))

    return vcat([sf_df, rf_df]...)[:, jks_eq_fd_cols]
end


# ** Joint Optimal VB
# include("_joint_equilibrium/_joint_optimal_vb.jl")
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

    
    # Set vbl bounds and form vbl grid ######################################
    s_vbh = NaN
    if get_obj_model(jf.sf) == "svm"
        s_vbh = get_cvm_vb(jf.sf, jf.sf.pm.sigmah; mu_b=jks.mu_b, c=jks.c, p=jks.p)
    end
    
    r_vbh = NaN
    if get_obj_model(jf.rf) == "svm"
        r_vbh = get_cvm_vb(jf.rf, jf.rf.pm.sigmah; mu_b=jks.mu_b, c=jks.c, p=jks.p)
    end

    if .&(isnan(s_vbh), isnan(r_vbh))
        s_vbh = get_cvm_vb(jf.sf, jf.sf.pm.sigmal; mu_b=jks.mu_b, c=jks.c, p=jks.p)
        r_vbh = get_cvm_vb(jf.rf, jf.rf.pm.sigmal; mu_b=jks.mu_b, c=jks.c, p=jks.p)
    end
    vbh_max = maximum([x for x in [s_vbh, r_vbh] if !isnan(x)])
    vbh_min = minimum([x for x in [s_vbh, r_vbh] if !isnan(x)])
    
    # Form Grid of VB candidates
    vbl_min = minimum(x -> isnan(x) ? Inf : x, [jks.fi_sf_vb, jks.fi_rf_vb, vbl_min])
    vbl_max = maximum(x -> isnan(x) ? -Inf : x, [jks.fi_sf_vb, jks.fi_rf_vb, vbl_max])

    # vbh_max / vbl_min < jf.rf.bi.vbhlmax
    vbl_min = maximum([lb2 * vbl_min, vbh_max/ jf.rf.bi.vbhlmax])

    # vbh_min / vbl_max > jf.rf.bi.vbhlmin
    vbl_max = minimum([ub2 * vbl_max, vbh_min/ jf.rf.bi.vbhlmin])
    
    # Form vbl grid
    vbl_grid = range(vbl_min, stop= vbl_max, length=vbN2)
    # #######################################################################


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
    s2 = joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
    
    # sfdf = vcat([s1, s2]...)
    # ###############################################################

    # Case 2: Risky Firm defaults first! ############################
    opt_sf_vb, opt_rf_vb = interp_optimal_vbs(jf, jks, rfdf)

    # vbl that sets E'_s(vbl) to zero
    sf_vb, rf_vb = get_type_contingent_vbs(opt_sf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=false)   
    r1 = joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)

    # vbl that sets E'_r(vbl) to zero   
    sf_vb, rf_vb = get_type_contingent_vbs(opt_rf_vb,
                                           jks.fi_sf_vb,
                                           jks.fi_rf_vb; 
                                           sf_defaults_first=false)
    r2 = joint_eq_fd(jf; jks=jks, sf_vb=sf_vb, rf_vb=rf_vb)
    if abs.(r2[isnan.(r2[:, :sf_vb]), :eq_deriv][1]) > 1e-2
        r2 = refine_contingent_vbs(jf, jks, sf_vb, rf_vb)
    end
    # ###############################################################
    # println("=====================================================")
    # println("=====================================================")
    # println(r2)
    println("=====================================================")
    println("=====================================================")
    # Compile Results
    println("compiling results")
    df =  vcat([s1, s2, r1, r2]...)
    df[!, :sf_defaults_first] .= vcat([fill(true, 4), fill(false, 4)]...)

    cols1 = [:sf_defaults_first,
             :fi_vb, :sf_vb, :rf_vb, :vb,
             :eq_deriv, :eq_min_val]
    return df[:, vcat(cols1, [x for x in names(df) if !(x in cols1)]...)]
end


function filter_joint_vb_results(df::DataFrame;
                                 tol1::Float64=-1e-2, tol2::Float64=1e-2)

    # Limited Liability Conditions ############################
    cond1 = df -> .&([df[:, x] .>= tol1 for x in [:eq_deriv, :eq_min_val]]...)

    # When the Safe Firm defaults first, 
    # its equity and equity derivative should be zero
    # at the joint default barrier
    sf_cond = df -> .&(df[:, :sf_defaults_first], 
                       isnan.(df[:, :sf_vb]) .| (df[:, :eq_deriv] .<= tol2))

    # When the Risky Firm defaults first, 
    # its equity and equity derivative should be zero
    # at the joint default barrier
    rf_cond = df -> .&(df[:, :sf_defaults_first] .==false, 
                       isnan.(df[:, :rf_vb]) .| (df[:, :eq_deriv] .<= tol2))

    # Limited Liability Conditions must be satisfied by both firms:
    llcond = df -> sum(.&(cond1(df), (sf_cond(df) .| rf_cond(df)))) == 2
    # #########################################################

    ## Find vb candidates
    # vbdf = by(df, :vb, llcond = names(df) => x -> llcond(x))
    # vblist = vbdf[vbdf[:, :llcond], :vb]
    # Filter DataFrame
    # if !isempty(vblist)
    #     return df[in.(df[:, :vb], vblist), :]
    
    LL = []
    for x in groupby(df, [:vb, :sf_defaults_first])
        if llcond(x)
            push!(LL, DataFrame(x))
        end
    end
    vbdf = vcat(LL...)
    
    if isempty(vbdf)
        println(string("No results for mu_b in ", unique(df[:, :mu_b])))
        return
    end
    return vbdf
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


# ** Joint Leverage Functions
# include("_joint_equilibrium/_joint_leverage_functions.jl")
function find_joint_payoffs(jf, jks, mu_b_grid;
                            sf_obj_fun::Symbol=:firm_value,
                            rf_obj_fun::Symbol=:MBR,
                            spline_k::Int64=3,
                            spline_bc::String="extrapolate",
                            N::Int64=10^5)


    dfl = @time fetch(@spawn [find_joint_optimal_vb(jf, jks;
                                                    mu_b=mu_b,
                                                    rerun_fi_vb=true)
                              for mu_b in mu_b_grid])
    
    # Store results in a DataFrame
    df = vcat(dfl...)
    sf_df = df[isnan.(df[:rf_vb]), :]
    rf_df = df[isnan.(df[:sf_vb]), :]
    
    
    # Interpolate Objective Functions in mu_b ################################
    sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    rf_objf = Dierckx.Spline1D(rf_df[:mu_b], rf_df[rf_obj_fun]; k=spline_k, bc=spline_bc)

    # Refine mu_b_grid
    ref_mu_b_grid = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=N)

    return sf_objf, rf_objf, ref_mu_b_grid
end




function find_mu_b_intervals(ref_mu_b_grid::StepRangeLen{Float64,
                                                         Base.TwicePrecision{Float64},
                                                         Base.TwicePrecision{Float64}},
                             cond::BitArray{1}; N::Int64=20)
    #ma = ref_mu_b_grid[cond]
    #md = abs.(ma[2:end] - ma[1:end-1])
    mu_b_discarded = ref_mu_b_grid[cond .== false]
    
   
    #maximum(md) - sum(md)/size(md, 1) < 1e-6
    mu_b_dis_min = minimum(mu_b_discarded)
    mu_b_dis_max = maximum(mu_b_discarded)

    # Check if there is only one interval
    if any([abs.(mu_b_dis_min - minimum(ref_mu_b_grid)) < 1e-6,
            abs.(mu_b_dis_max - maximum(ref_mu_b_grid)) < 1e-6])
        filtered_mu_b_grid_1 = range(minimum(ref_mu_b_grid[cond]),
                                     stop=maximum(ref_mu_b_grid[cond]),
                                     length=N)
        filtered_mu_b_grid_2 = range(0, stop=0, length=0)
    else
        # int1_up  = ref_mu_b_grid[argmax(md) + 1]
        # int2_down = reverse(ref_mu_b_grid)[argmax(reverse(md)) + 1]
        # filtered_mu_b_grid_1 = range(minimum(mu_b_grid), stop=int1_up, length=10)
        # filtered_mu_b_grid_2 = range(int2_down, stop=maximum(mu_b_grid), length=10)
        
        # In case there are two intervals:
        filtered_mu_b_grid_1 = range(minimum(ref_mu_b_grid[cond]),
                                     stop=minimum(mu_b_discarded),
                                     length=N)
        filtered_mu_b_grid_2 = range(maximum(mu_b_discarded),
                                     stop=maximum(ref_mu_b_grid[cond]),
                                     length=N)
    end

    return filtered_mu_b_grid_1, filtered_mu_b_grid_2
end


function find_opt_mu_b(sf_objf,
                       filtered_mu_b_grid_1::StepRangeLen{Float64,
                                                          Base.TwicePrecision{Float64},
                                                            Base.TwicePrecision{Float64}},
                       filtered_mu_b_grid_2::StepRangeLen{Float64,
                                                          Base.TwicePrecision{Float64},
                                                          Base.TwicePrecision{Float64}})
    
    # Maximize Safe Firm's Objective Function
    mu_b_opt_1 = filtered_mu_b_grid_1[argmax(sf_objf(filtered_mu_b_grid_1))]
    
    if size(filtered_mu_b_grid_2, 1) == 0
        return mu_b_opt_1
    else
        # Maximize Safe Firm's Objective Function
        mu_b_opt_2 = filtered_mu_b_grid_2[argmax(sf_objf(filtered_mu_b_grid_2))]
        
        return maximum([mu_b_opt_1, mu_b_opt_2])
    end
        
        
end


    # # Now compute payoffs for the safe firm 
    # s_eqdf_list = fetch(@spawn [sep_eq_fd(jf.sf, deepcopy(jks), mu_b) 
    #                             for mu_b in filtered_mu_b_grid])
    # sf_df = vcat(s_eqdf_list...)

    # sf_objf = Dierckx.Spline1D(sf_df[:mu_b], sf_df[sf_obj_fun]; k=spline_k, bc=spline_bc)
    
    # # Maximize Safe Firm's Objective Function
    # ref_mu_b_grid = range(minimum(filtered_mu_b_grid),
    #                         stop=maximum(filtered_mu_b_grid), length=N2)
    # mu_b_opt = ref_mu_b_grid[argmax(sf_objf(ref_mu_b_grid))]


function sep_eq_fd(sf, jks, mu_b::Float64)
    sf_vb = find_full_info_vb(sf, jks; mu_b=mu_b)
    return eq_fd(sf, vbl=sf_vb, 
                 mu_b=mu_b, m=jks.m, 
                 c=jks.c, p=jks.p)
    
end


function separate_eq_calculator(jf, jks, mu_b_opt::Float64, fi_rf_mu_b::Float64)
    s_eqdf = sep_eq_fd(jf.sf, deepcopy(jks), mu_b_opt)
    r_eqdf = sep_eq_fd(jf.rf, deepcopy(jks), fi_rf_mu_b)

    s_eqdf[:fi_vb] = s_eqdf[1, :vb]
    s_eqdf[:sf_vb] = s_eqdf[1, :vb]
    s_eqdf[:rf_vb] = NaN
    
    r_eqdf[:fi_vb] = r_eqdf[1, :vb]
    r_eqdf[:sf_vb] = NaN 
    r_eqdf[:rf_vb] = r_eqdf[1, :vb]

    eqdf = vcat(s_eqdf, r_eqdf)
    eqdf[:sf_defaults_first] = s_eqdf[1, :vb] > r_eqdf[1, :vb]
    eqdf[:eq_type] = "separating"
    eqdf[:mu_s] = jks.mu_s

    return eqdf
end


# ** Joint Optimal Bond Measure
# include("_joint_equilibrium/_joint_optimal_bond_measure.jl")
function get_pool_eqdf(jf, jks, mu_b_grid,
                       fi_rf_obj_val::Float64;
                       sf_obj_fun::Symbol=:firm_value,
                       rf_obj_fun::Symbol=:MBR,
                       spline_k::Int64=3,
                       spline_bc::String="extrapolate",
                       N1::Int64=20,
                       N2::Int64=10^5)

    sf_objf, rf_objf, ref_mu_b_grid = find_joint_payoffs(jf, jks, mu_b_grid;
                                                         sf_obj_fun=sf_obj_fun,
                                                         rf_obj_fun=rf_obj_fun,
                                                         spline_k=spline_k,
                                                         spline_bc=spline_bc,
                                                         N=N2)

    # Range of mu_b values in pooling equilibrium
    cond = rf_objf(ref_mu_b_grid) .>= fi_rf_obj_val
    filtered_mu_b_grid_1, filtered_mu_b_grid_2 = find_mu_b_intervals(ref_mu_b_grid,
                                                                     cond; N=N1)

    # Maximize Safe Firm's Objective Function
    mu_b_opt = find_opt_mu_b(sf_objf,
                             filtered_mu_b_grid_1,
                             filtered_mu_b_grid_2)
    
    # Compute Results
    return find_joint_optimal_vb(jf, jks;
                                 mu_b=mu_b_opt, rerun_fi_vb=true)
end



function get_sep_eqdf(jf, jks, mu_b_grid,
                      fi_rf_mu_b::Float64,
                      fi_rf_obj_val::Float64;
                      sf_obj_fun::Symbol=:firm_value,
                      rf_obj_fun::Symbol=:MBR,
                      spline_k::Int64=3,
                      spline_bc::String="extrapolate",
                      N1::Int64=20,
                      N2::Int64=10^5)

    # Compute the Payoffs in case of Misrepresentation
    sep_misrep_jks = deepcopy(jks)
    sep_misrep_jks.mu_s = 1.
    sf_objf, rf_objf, ref_mu_b_grid = find_joint_payoffs(jf, sep_misrep_jks,
                                                         mu_b_grid;
                                                         sf_obj_fun=sf_obj_fun,
                                                         rf_obj_fun=rf_obj_fun,
                                                         spline_k=spline_k,
                                                         spline_bc=spline_bc,
                                                         N=N2)

    # Filter -> leverage values for which misrepresentation is not
    # attractive
    cond = rf_objf(ref_mu_b_grid) .<= fi_rf_obj_val
    filtered_mu_b_grid_1, filtered_mu_b_grid_2 = find_mu_b_intervals(ref_mu_b_grid,
                                                                     cond; N=N1)
    
    # Mu_b that yields maximum Safe Type's Firm Value
    # conditional on misrepresentation not being optimal for risky type.
    # Since mu_s = 1 above, the payoffs for the safe firm coincide with
    # the payoffs under full information, when debt investors can fully
    # observe firms' types.
    mu_b_opt = find_opt_mu_b(sf_objf,
                             filtered_mu_b_grid_1,
                             filtered_mu_b_grid_2)
    

    # Compute Separating Eq. Payoffs
    eqdf = separate_eq_calculator(jf, jks, mu_b_opt, fi_rf_mu_b)

    return eqdf
end


function find_joint_optimal_bond_measure(jf, jks,
                                         fi_sf_mu_b::Float64,
                                         fi_rf_mu_b::Float64,
                                         fi_rf_obj_val::Float64;
                                         equilibrium_type::String="all", 
                                         sf_obj_fun::Symbol=:firm_value,
                                         rf_obj_fun::Symbol=:MBR,
                                         lb::Float64=.75,
                                         ub::Float64=1.25,
                                         mu_bN::Int64=20,
                                         mu_bN2::Int64=10^5,
                                         spline_k::Int64=3,
                                         spline_bc::String="extrapolate")

    # Form Grid of mu_b candidates
    min_mu_b = .75 * minimum([fi_sf_mu_b, fi_rf_mu_b])
    max_mu_b = 1.25 * maximum([fi_sf_mu_b, fi_rf_mu_b])
    mu_b_grid = range(min_mu_b, stop=max_mu_b, length=mu_bN)

    # Compute Optimal VB for each mu_b candidate
    jks2 = deepcopy(jks)

    # Misrepresentation case:
    eq_types = ["pooling", "separating"]
    if equilibrium_type == "all"
        println("Computing Pooling and Separating Equilibria")
    elseif equilibrium_type == "pooling"
        println("Computing Pooling Equilibrium")
        eq_types = ["pooling"]
    elseif equilibrium_type == "separating"
        println("Computing Separating Equilibrium")
        eq_types = ["separating"]
    else
        println("Please set equilibrium_type to pooling or separting. Exiting...")
        return
    end
    
    optdf = DataFrame()
    for eq_type in eq_types
        # Find optimal mu_b values 
        if eq_type == "pooling"
            eqdf = get_pool_eqdf(jf, jks2, mu_b_grid,
                                 fi_rf_obj_val;
                                 sf_obj_fun=sf_obj_fun,
                                 rf_obj_fun=rf_obj_fun,
                                 N1=mu_bN,
                                 N2=mu_bN2)
            
        elseif eq_type == "separating" 
            eqdf = get_sep_eqdf(jf, jks2, mu_b_grid,
                                fi_rf_mu_b, fi_rf_obj_val;
                                sf_obj_fun=sf_obj_fun,
                                rf_obj_fun=rf_obj_fun,
                                N1=mu_bN,
                                N2=mu_bN2)

        end
        
        # Add Objective Function and Equilibrium Type columns
        eqdf[:obj_fun] = String(sf_obj_fun)
        eqdf[2, :obj_fun] = String(rf_obj_fun)
        eqdf[:eq_type] = eq_type


        optdf = vcat(optdf, reshape_sf_rf_df(eqdf))
    end
    
    return optdf
end
# ##########################################################

# * END MODULE
end
