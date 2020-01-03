module_path = "/home/artur/BondPricing/Julia/modules/"
push!(LOAD_PATH, module_path)


module JointEqStructs

using DataFrames
using Dates

using ModelObj: Firm, set_opt_k_struct

using BondPrInterp: get_cvm_debt_price
    
using EqFinDiff: eq_fd

using Batch: BatchObj,
             get_bt,
             get_bt_mobj,
             load_cvm_opt_results_df,
             load_svm_opt_results_df,
             opt_k_struct_df_name,
             opt_k_struct_df_coltypes,
             get_batch_comb_numbers,
             comb_folder_dict


# * To be removed Struct
mutable struct FirmSpecificParams
    iota::Float64
    lambda::Float64
    sigmal::Float64
    sigmah::Float64
end


mutable struct MarketTypeDist
    mu_s::Float64
    kappa::Float64
end


# * Structs
mutable struct RMPCFirmObj
    fr  # Firm
    bt  # BatchStruct
end


mutable struct TypesDist
    mu_s::Float64
    st_iota::Float64
    st_lambda::Float64
    st_sigmah::Float64
    rt_iota::Float64
    rt_lambda::Float64
    rt_sigmah::Float64
end


mutable struct FirmType
    rm::RMPCFirmObj
    nrm::RMPCFirmObj
end


mutable struct BondContract
    m::Float64
    c::Float64
    p::Float64
end


mutable struct TypesCommonParams
    V0::Float64
    alpha::Float64
    pi::Float64
    r::Float64
    gross_delta::Float64
    xi::Float64
    sigmal::Float64
end


mutable struct JointFirms
    kappa::Float64
    td::TypesDist
    cp::TypesCommonParams
    st::FirmType
    rt::FirmType
    bc::BondContract
end


mutable struct JointKStruct
    mu_s::Float64
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64

    # Safe Firm
    fi_st_vb::Float64
    st_vb::Float64

    # Risk Firm
    fi_rt_vb::Float64
    rt_vb::Float64
end


# * Joint Inputs
# Directories and File Names
ep_dir = "EP"
fi_dir = "FI"
jeq_dir = "JEQ"
fidf_name = "fidf"
misrepdf_name = "misrepdf"
pooldf_name = "pooldf"
sepdf_name = "sepdf"


# DataFrames Column Types
fidf_col_types = vcat(fill(Float64, 9), [Bool], fill(Float64, 17)) # remove
# Misrepresentation, Pool and Separating DFs Column Types
mps_col_types = vcat(String, Bool, fill(Float64, 11), Bool,
                     fill(Float64, 10), String, fill(Float64, 5),
                     Bool, fill(Float64, 10), String,
                     fill(Float64, 7)) # remove

# Electronic Market Parameters
epmcols = [:kappa, :mu_s, :m, :c, :p]
epm_eq_cols = [:obj_fun, :sf_defaults_first]

# Common Parameters
commoncols = vcat(fieldnames(TypesCommonParams)...)

# Firm Specific Parameters
fspcols = vcat(fieldnames(FirmSpecificParams)...)

# Default Boundaries
vbcols = [:fi_vb, :sf_vb, :rf_vb, :vb]


# DATAFRAME COLUMNS
# K Structure
jks_cols = [:mu_s, :m, :mu_b, :c, :p]

# Default Barrier
vb_cols = [:fi_vb, :st_vb, :rt_vb, :vb]

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



# * Auxiliary Functions
function round_value(x::Float64)
    xr_vec = [floor(x), (ceil(x) + floor(x))/2, ceil(x)]
    diffs = [x-xr_vec[1], abs(x-xr_vec[2]), xr_vec[3] - x]    
    return xr_vec[argmin(diffs)]
end


function type_fun(x::Symbol)
    if any([occursin("neg", String(x)),
            occursin("defaults_first", String(x))])
        return Bool
    elseif occursin("time", String(x))
        return DateTime
    elseif any([occursin("type", String(x)),
                String(x) == "rmp"])
        return Symbol
    end
    return Float64
end


# * Initialization Functions
function get_joint_types_dfs(pardict::Dict{Symbol,Array{Float64,1}},
                             iota_vec::Array{Float64,1},
                             lambda_vec::Array{Float64,1},
                             sigmah_vec::Array{Float64,1};
                             firm_obj_fun::Symbol=:firm_value)
    
    cvmdict = deepcopy(pardict)
    cvmdict[:iota] = [x for x in iota_vec if !isnan(x)]
    
    svmdict = deepcopy(pardict)
    svmdict[:iota] = [.0]
    svmdict[:lambda] = [x for x in lambda_vec if !isnan(x)]
    svmdict[:sigmah] = [x for x in sigmah_vec if !isnan(x)]


    # Load CVM and SVM Data
    cbt = get_bt(;model="cvm", m=cvmdict[:m][1], m_comb_num=1)
    cvm_data = load_cvm_opt_results_df(; m=cvmdict[:m][1],
                                       firm_obj_fun=firm_obj_fun)
    
    sbt = get_bt(;model="svm", m=svmdict[:m][1], m_comb_num=1)
    svm_data = load_svm_opt_results_df(sbt;
                                       m=svmdict[:m][1],
                                       firm_obj_fun=firm_obj_fun)

    # Get Combination Numbers
    cvm_combs = get_batch_comb_numbers(cbt, cvmdict)[:, :comb_num]
    svm_combs = get_batch_comb_numbers(sbt, svmdict)[:, :comb_num]
    # #########################################################
    
    
    # CVM and SVM Firms' DFs ###############################
    cvmdf = cvm_data[[(x in cvm_combs) for x in cvm_data[:, :comb_num]], :]
    svmdf = svm_data[[(x in svm_combs) for x in svm_data[:, :comb_num]], :]
    # #########################################################


    return cvmdf, svmdf
end


function get_rm_nrm_comb_nums(iota::Float64, 
                              lambda::Float64, sigmah::Float64, 
                              cvmdf::DataFrame, svmdf::DataFrame)
    
    rm_comb_num = cvmdf[abs.(cvmdf[:, :iota] .- iota) .< 1e-5, 
                        :comb_num][1]

    nrm_comb_num = 0
    if .&(!isnan(lambda), !isnan(sigmah))
        svmdf_slice = svmdf[.&(abs.(svmdf[:, :lambda] .- lambda) .< 1e-5,
                                   abs.(svmdf[:, :sigmah] .- sigmah) .< 1e-5), 
                            :comb_num]
        nrm_comb_num = !isempty(svmdf_slice) ? svmdf_slice[1] : 0
    end
        
                            
    return rm_comb_num, nrm_comb_num
end


function firm_type_constructor(; cvm_comb_num::Int64=0,
                               svm_comb_num::Int64=0,
                               firm_obj_fun::Symbol=:firm_value,
                               set_rmpc_opt_k_struct::Bool=true,
                               cvmdf::DataFrame=DataFrame(),
                               svmdf::DataFrame=DataFrame(),                               
                               opt_k_struct_df_name::String=opt_k_struct_df_name,
                               recompute_svm::Bool=false)

    # Risk-Management: Constant Volatility Model
    if cvm_comb_num > 0
        rm_bt, rmf = get_bt_mobj(;model="cvm", comb_num=cvm_comb_num)
    end
    

    # No Risk-Management: Stochastic Volatility Model
    if svm_comb_num > 0
        nrm_bt, nrmf = get_bt_mobj(;model="svm", comb_num=svm_comb_num)
    end


    # Verify that the transaction costs are the same
    if .&(cvm_comb_num > 0, svm_comb_num > 0)
        if abs(rmf.pm.kappa - nrmf.pm.kappa) > 1e-5
            println("Error! Inconsistent transaction costs (kappa). Exiting...")
            return
        end
    end


    # Set RMP-Contingent Optimal Capital Structures ########
    if .&(set_rmpc_opt_k_struct, cvm_comb_num > 0)
        if isempty(cvmdf)
            cvmdf = load_cvm_opt_results_df(; m=rmf.m, firm_obj_fun=firm_obj_fun)
        end
        
        rmf = set_opt_k_struct(rmf, cvmdf)
    end
    
    if .&(set_rmpc_opt_k_struct, svm_comb_num > 0)
        if isempty(svmdf)
            svmdf = load_svm_opt_results_df(nrm_bt; m=nrmf.m,
                                            firm_obj_fun=firm_obj_fun,
                                            opt_k_struct_df_name=opt_k_struct_df_name,
                                            recompute=recompute_svm)
        end
        
        nrmf = set_opt_k_struct(nrmf, svmdf)
    end
    # #####################################################


    rm = cvm_comb_num > 0 ? RMPCFirmObj(rmf, rm_bt) : RMPCFirmObj(nothing, nothing)
    nrm = svm_comb_num > 0 ? RMPCFirmObj(nrmf, nrm_bt) : RMPCFirmObj(nothing, nothing)
    return FirmType(rm, nrm)
end


function get_firm_type(iota::Float64, 
                       lambda::Float64, sigmah::Float64, 
                       cvmdf::DataFrame, svmdf::DataFrame;
                       set_rmpc_opt_k_struct::Bool=true)

    # Get Combination Numbers
    rm_comb_num, nrm_comb_num = get_rm_nrm_comb_nums(iota, 
                                                     lambda, sigmah,
                                                     cvmdf, svmdf)

    # Construct Type
    return firm_type_constructor(; cvm_comb_num=rm_comb_num,
                                 svm_comb_num=nrm_comb_num,
                                 set_rmpc_opt_k_struct=set_rmpc_opt_k_struct,
                                 cvmdf=cvmdf, svmdf=svmdf)                            
end


function get_ep_bond_contract(ft::FirmType)
    ft_opt_rmp = get_otc_fi_opt_rmp(ft)
    ft_optKS = ft_opt_rmp == :rm ? ft.rm.fr.optKS : ft.nrm.fr.optKS
    ep_c = round_value(ft_optKS.c)
    ep_ratio = round_value(ft_optKS.p/ft_optKS.c)
    ep_p = ep_ratio * ep_c

   return BondContract(ft_optKS.m, ep_c, ep_p)
end


function form_joint_types(kappa::Float64, mu_s::Float64,
                          pardict::Dict{Symbol,Array{Float64,1}},
                          typesdict::Dict{Symbol,Float64};
                          cvmdf::DataFrame=DataFrame(),
                          svmdf::DataFrame=DataFrame(),
                          firm_obj_fun::Symbol=:firm_value,
                          bc=nothing)

    # Type Distribution
    td = TypesDist(mu_s,
                   typesdict[:st_iota],
    	           typesdict[:st_lambda],
    		   typesdict[:st_sigmah],
    		   typesdict[:rt_iota],
    		   typesdict[:rt_lambda],
    		   typesdict[:rt_sigmah])

    # Get Safe and Risky Firms' Full Info Optimal Results #####
    if any([isempty(cvmdf), isempty(svmdf)])
            cvmdf, svmdf = get_joint_types_dfs(pardict,
                                               [td.st_iota, td.rt_iota],
                                               [td.st_lambda, td.rt_lambda],
                                               [td.st_sigmah, td.rt_sigmah])
    end
    # #########################################################

    # Form Types
    st = get_firm_type(td.st_iota, td.st_lambda, td.st_sigmah,
                       cvmdf, svmdf;
                       set_rmpc_opt_k_struct=true)

    rt = get_firm_type(td.rt_iota, td.rt_lambda, td.rt_sigmah,
                       cvmdf, svmdf;
                       set_rmpc_opt_k_struct=true)

    # Common Parameters
    cp = TypesCommonParams(st.rm.fr.pm.V0,
                           st.rm.fr.pm.alpha,
                           st.rm.fr.pm.pi,
                           st.rm.fr.pm.r,
                           st.rm.fr.pm.gross_delta,
                           st.rm.fr.pm.xi,
                           st.rm.fr.pm.sigmal)

    # Form Bond Contract
    if isnothing(bc)
        bc = get_ep_bond_contract(st)
    end

    return JointFirms(kappa, td, cp, st, rt, bc)
end


# ** Set and Get Methods
function get_otc_fi_opt_rmp(ft::FirmType)
    if .&(!isnothing(ft.rm.fr), isnothing(ft.nrm.fr))
        opt_rmp = :rm
    elseif .&(!isnothing(ft.rm.fr), !isnothing(ft.nrm.fr))
        rm_debt = get_cvm_debt_price(ft.rm.fr, ft.rm.fr.optKS.vbl, 
                                     ft.rm.fr.pm.sigmal;
                                     mu_b=ft.rm.fr.optKS.mu_b, 
                                     m=ft.rm.fr.optKS.m, 
                                     c=ft.rm.fr.optKS.c,    
                                     p=ft.rm.fr.optKS.p)

        rm_eq = get_cvm_eq(ft.rm.fr, ft.rm.fr.pm.V0,
                           ft.rm.fr.pm.sigmal;
                           mu_b=ft.rm.fr.optKS.mu_b, 
                           m=ft.rm.fr.optKS.m, 
                           c=ft.rm.fr.optKS.c,    
                           p=ft.rm.fr.optKS.p)
        
        rm_fv = rm_debt + rm_eq
        
        nrm_res = eq_fd(ft.nrm.fr; vbl=ft.nrm.fr.optKS.vbl, 
                        mu_b=ft.nrm.fr.optKS.mu_b, 
                        m=ft.nrm.fr.optKS.m, 
                        c=ft.nrm.fr.optKS.c,    
                        p=ft.nrm.fr.optKS.p)
        
        nrm_fv = nrm_res[:firm_value]
        
        opt_rmp = rm_fv > nrm_fv ? :rm : :nrm
    else
        println("Error! Missing Risk-Management results. Exiting...")
        return
    end
    
    return opt_rmp
end


function get_types_common_params(ft::FirmType)
        return TypesCommonParams(ft.rm.fr.pm.V0,
                                 ft.rm.fr.pm.alpha,
                                 ft.rm.fr.pm.pi,
                                 ft.rm.fr.pm.r,
                                 ft.rm.fr.pm.gross_delta,
                                 ft.rm.fr.pm.xi,
                                 ft.rm.fr.pm.sigmal)
end


# * Joint Capital Struct
# df is the dataframe with the Full Information Equilibrium results
function form_firms_jks(jf, df::DataFrame;
                        jks::JointKStruct=JointKStruct(fill(NaN, 10)...),
                        mu_s::Float64=NaN,
                        set_fi_vb::Bool=true)

    # Transaction Costs
    cond = .&(size(unique(df[:, :kappa]), 1) == 1,  
              abs.(unique(df[:, :kappa])[1] - jf.kappa < 1e-5))
    if !cond
        println(string("Error! Transaction cost parameter kappa in", 
                       " Joint Firm Struct and dataframe do not coincide!"))
        println("Exiting...")
        return
    end

    # Common Parameters
    for x in fieldnames(TypesCommonParams)
        cond = .&(cond, size(unique(df[:, x]), 1) == 1, 
                  abs.(unique(df[:, x])[1] - getfield(jf.cp, x)) < 1e-5)
        if !cond
            println(string("Error! Common parameters in Joint Firm Struct ", 
                           " and dataframe do not coincide!"))
            println("Exiting...")
            return
        end
    end

    # Bond Contract
    for x in fieldnames(BondContract)
        cond = .&(cond, size(unique(df[:, x]), 1) == 1, 
                  abs.(unique(df[:, x])[1] - getfield(jf.bc, x)) < 1e-5)
        if cond
            setfield!(jks, x, getfield(jf.bc, x))
        else
            println(string("Error! Bond Contract parameters in Joint Firm Struct",
                           " and dataframe do not coincide!"))
            println("Exiting...")
            return
        end
    end

    firms=Dict{Symbol, Firm}()
    # Safe and Risky Types
    for ftype in [:st, :rt]
        # Form Firm Type
        ftdf = df[df[:, :type] .== ftype, :]
        cond = .&(cond, size(ftdf, 1) == 1)
        if !cond
            println(string("Error! Type", ft, "sub-dataframe has more than one row."))
            println("Exiting...")
        end

        ft_rmp = ftdf[1, :rmp]
        ft = getfield(getfield(jf, ftype), ft_rmp).fr
        for x in [:iota, :lambda, :sigmah]    
            xval = ftdf[1, x]
            if isnan(getfield(ft.pm, x))
                cond = .&(cond, isnan(xval))
            elseif !isnan(xval)
                cond = .&(cond, abs.(getfield(ft.pm, x) - xval) < 1e-5)
            else
                cond = false
            end

            if !cond
                println(string("Error! Type ", ftype, " parameter ", x, " in Joint Firm",
                           "Struct does not coincide with value in dataframe."))
                println("Exiting...")
                return
            end

        end
        
        # Set Full Information VB
        if set_fi_vb
            setfield!(jks, Symbol(:fi_, ftype, :_vb), ftdf[1, :vb])
        end

        # Store firm
        firms[ftype] = ft
    end
    
    # Set measure of safe types
    mu_s = isnan(mu_s) ? jf.td.mu_s : mu_s
    setfield!(jks, :mu_s, mu_s)
                
    return firms[:st], firms[:rt], jks
end


function update_jks(jks;
                    mu_s::Float64=NaN,
                    mu_b::Float64=NaN,
                    m::Float64=NaN,
                    c::Float64=NaN,
                    p::Float64=NaN,
                    fi_st_vb::Float64=NaN,
                    fi_rt_vb::Float64=NaN,
                    st_vb::Float64=NaN,
                    rt_vb::Float64=NaN,
                    vbl::Float64=NaN)
    
    if !isnan(mu_s)
        setfield!(jks, :mu_s, mu_s)
    end
    
    if !isnan(mu_b)
        setfield!(jks, :mu_b, mu_b)
    end
    
    if !isnan(m)
        setfield!(jks, :m, m)
    end

    if !isnan(c)
        setfield!(jks, :c, c)
    end

    if !isnan(p)
        setfield!(jks, :p, p)
    end

    if !isnan(fi_st_vb)
        setfield!(jks, :fi_st_vb, fi_st_vb)
    end

    if !isnan(fi_rt_vb)
        setfield!(jks, :fi_rt_vb, fi_rt_vb)
    end
    
    if !isnan(st_vb)
        setfield!(jks, :st_vb, st_vb)
    end
    
    if !isnan(rt_vb)
        setfield!(jks, :rt_vb, rt_vb)
    end

    if !isnan(vbl)
        setfield!(jks, :vbl, vbl)
    end

   return jks
end
    
# * END
end
