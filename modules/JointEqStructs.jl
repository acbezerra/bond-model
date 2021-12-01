module_path = "/home/artur/BondPricing/bond-model/modules/"
push!(LOAD_PATH, module_path)


module JointEqStructs

using DataFrames
using Dates
using CSV

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
             comb_folder_dict,
             dfcols


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
# vbcols = [:fi_vb, :sf_vb, :rf_vb, :vb]



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
    "full_info" => Dict{Symbol, String}(:dfn => string("z", fidf_name),
                                        :fn_prefix => "fi"),
    "misrep" => Dict{Symbol, String}(:dfn => string("z", misrepdf_name),
                                       :fn_prefix =>  "misrep"),
    "pooling" => Dict{Symbol, String}(:dfn => string("z", pooldf_name),
                                     :fn_prefix => "pool"),
    "separating" => Dict{Symbol, String}(:dfn => string("z", sepdf_name),
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


function is_svm(jf)
        cvm_cond = false
        for ft in [:st, :rt]
            fr = getfield(getfield(jf, ft), :nrm).fr
            if !isnothing(fr)
                 cvm_cond = cvm_cond | .&(fr.pm.sigmah > fr.pm.sigmal, isnan(fr.pm.lambda))
            end
        end
        return !cvm_cond
end


# * DataFrame Columns and Column-Types
# ** DataFrame Columns
# DATAFRAME COLUMNS
# K Structure
jks_cols = [:mu_s, :m, :mu_b, :c, :p]

# Default Barrier
vb_cols = [:fi_vb, :st_vb, :rt_vb, :vb]

# EFD
share_cols = [:eq_deriv, :eq_deriv_min_val,
              :eq_min_val, :eq_negative,
              :eq_vb, :MBR, :debt, :equity,
              :firm_value, :leverage, :yield, :yield_spd]

# Parameters
param_cols = [:iota, :lambda, :sigmah,
              :gross_delta, :delta, :kappa,
              :sigmal, :V0, :xi, :r, :alpha, :pi]

jks_eq_fd_cols = vcat(jks_cols, vb_cols, share_cols, param_cols)

extra_cols = [:eq_type, :datetime, :type, :rmp]

rmp_param_cols = [:rm_iota, :nrm_lambda, :nrm_sigmah]
jeq_rmp_param_cols = vcat([Symbol(x, :_, y) for y in rmp_param_cols, x in [:s, :r]]...)


# fi_df_names =  vcat([x for x in dfcols if !occursin("diff", String(x))],
#                     [:eq_type, :datetime, :type, :rmp])
# misrep_df_names = vcat([unique(vcat([fi_df_names, jks_eq_fd_cols]...)),
#                         [:sf_defaults_first, :misrep_type]]...)


function get_type_identifiers(jf, ft::Symbol)

    if !(ft in [:st, :rt])
        println(string("Wrong type of argument. Firm type (ft) must be :st or :rt",
                       "\n Exiting..."))
        return
    end

    rm_iota = 0.0
    fr = getfield(getfield(jf, ft), :rm).fr
    if !isnothing(fr)
        rm_iota = fr.pm.iota
    end

    # No Risk Management
    nrm_lambda = NaN
    nrm_sigmah = NaN

    fr = getfield(getfield(jf, ft), :nrm).fr
    if !isnothing(fr)
        nrm_lambda = fr.pm.lambda
        nrm_sigmah = fr.pm.sigmah
    end

    return Dict{Symbol, Float64}(:rm_iota => rm_iota,
                                 :nrm_lambda => nrm_lambda,
                                 :nrm_sigmah => nrm_sigmah)
end


function get_joint_identifiers(jf)
    dd = Dict{Symbol, Float64}()
    for ft in [:st, :rt]
        # Type Prefix
        tp = ft == :st ? :s_ : :r_

        # Get Parameters
        pd = get_type_identifiers(jf, ft)

        dd[Symbol(tp, :rm_iota)] = pd[:rm_iota]
        dd[Symbol(tp, :nrm_lambda)] = pd[:nrm_lambda]
        dd[Symbol(tp, :nrm_sigmah)] = pd[:nrm_sigmah]
    end

    return dd
end


function form_cols_lists(dfnames::Array{Symbol, 1};
                         svm::Bool=true,
                         eq_type::Symbol=:misrep)

    # Common Columns ############################
    common_cols = commoncols

    # Case in which both firms are CVM
    if svm == false
        common_cols = [x for x in common_cols if x != :sigmal]
    end
    # ###########################################

    # Extra and Excluded Columns ################
    extra_cols = [:datetime, :eq_type, epmcols...]
    excluded_cols = [:type]

    rmp_cols = rmp_param_cols
    if eq_type != :full_info
        extra_cols = vcat(extra_cols, :sf_defaults_first)
        excluded_cols = vcat(excluded_cols, :vb)
    end

    if eq_type == :misrep
        extra_cols = vcat(extra_cols, :misrep_type)
        rmp_cols = jeq_rmp_param_cols
    elseif eq_type in [:pool, :pooling, :sep, :separating]
        jeq_cols = [:misrep_type, :misrep_ic_objf, :misrep_ic_objf_val, :jeq_objf]
        extra_cols = vcat(extra_cols, jeq_cols)
        rmp_cols = jeq_rmp_param_cols
    end

    # ###########################################

    # List of firm-specific variables ###########
    specific_cols = [col for col in dfnames if
                     !(col in vcat(excluded_cols,
                                   extra_cols,
                                   common_cols,
                                   rmp_cols))]

    return Dict{Symbol, Array{Symbol, 1}}(:common_cols => Symbol.(common_cols),
                                          :extra_cols => Symbol.(extra_cols),
                                          :excluded_cols => Symbol.(excluded_cols),
                                          :specific_cols => Symbol.(specific_cols),
                                          :rmp_cols => Symbol.(rmp_cols))
end


# ** DataFrame Types
function type_fun(x::Symbol)
    if any([occursin("neg", String(x)),
            occursin("defaults_first", String(x))])
        return Bool
    elseif occursin("time", String(x))
        return DateTime
    elseif any([occursin("type", String(x)),
                occursin("rmp", String(x)),
                String(x) == "jeq_objf",
                String(x) == "misrep_ic_objf"])
        return Symbol
    end

    return Float64
end


function get_df_cols_types(eqtype::Symbol; svm::Bool=true)
    all_cols = vcat(jks_cols, vb_cols, share_cols,
                    param_cols, extra_cols)

    colsl = Array{Symbol,1}()
    if eqtype in [:fi, :full_info]
        excluded_cols = [:mu_s, :fi_vb, :st_vb, :rt_vb]
        all_cols = vcat(all_cols, rmp_param_cols)
        colsl = [x for x in all_cols if !(x in excluded_cols)]
    else
        if eqtype == :misrep
            excluded_cols = [:mu_s]
            jeq_cols = []
        else
            excluded_cols = []
            jeq_cols = [:jeq_objf, :misrep_ic_objf,
                        :misrep_ic_objf_val]
        end

        tmp = Array{Symbol, 1}([x for x in all_cols if
                                !(x in excluded_cols)])
        colsd = form_cols_lists(tmp; svm=svm)
        renamer(ft, x) = (x in [:st_vb, :rt_vb]) ? x : Symbol(ft, :_, x)
        sfcols = [x for x in colsd[:specific_cols] if !(x == :rt_vb)]
        sf_cols = [renamer(:s, x) for x in sfcols]
        rfcols = [x for x in colsd[:specific_cols] if !(x == :st_vb)]
        rf_cols = [renamer(:r, x) for x in rfcols]

        colsl = vcat(colsd[:extra_cols], jeq_cols,
                     sf_cols, rf_cols, colsd[:common_cols],
                     jeq_rmp_param_cols)
    end

    coltypes = Dict([x => type_fun(x) for x in colsl])

    return coltypes
end

# ** Empty DataFrames
# function get_empty_df(; col_names::Array{Symbol,1}=fi_df_names)
#    tmp = [eval(type_fun(x)[]) for x in col_names]
#    return DataFrame(tmp, col_names)
#end
function get_empty_df(eq_type::Symbol; svm::Bool=true)
    colsd = get_df_cols_types(eq_type; svm=svm)
    tmp = [eval(type_fun(x)[]) for x in keys(colsd)]
    return DataFrame(tmp, [x for x in keys(colsd)])
end


# ** DataFrame Reshape
# Reshape Misrep, Pool and Separating Eq DFs #######################################
function reshape_sf_rf_df(df::DataFrame; svm::Bool=true) # - ADJUSTED
  colsd = form_cols_lists(Symbol.(names(df)); svm=svm,
                           eq_type=Symbol(df[1, :eq_type]))

    extradf = DataFrame(df[1, colsd[:extra_cols]])
    sfdf = DataFrame(df[1, [x for x in colsd[:specific_cols] if !(x == :rt_vb)]])
    rfdf = DataFrame(df[2, [x for x in colsd[:specific_cols] if !(x == :st_vb)]])

    # Common Parameters
    commondf =  DataFrame(df[1, colsd[:common_cols]])

    # Safe Firm Results
    renamer(ft, x) = (x in [:st_vb, :rt_vb]) ? x : Symbol(ft, :_, x)
    rename!(sfdf, [renamer(:s, x) for x in names(sfdf)])
    rename!(rfdf, [renamer(:r, x) for x in names(rfdf)])

    # RMP Cols
    rmp_params_df = DataFrame(df[1, colsd[:rmp_cols]])

    return hcat(extradf, sfdf, rfdf, commondf, rmp_params_df)
end


function reshape_joint_df(df::DataFrame; svm::Bool=true)
    rdf = DataFrame()
    for j in range(1, stop=floor(size(df, 1)/2))
        tmp = df[2*Int(j) - 1:2*Int(j), :]
        rdf = vcat([rdf, reshape_sf_rf_df(tmp; svm=svm)]...)
    end

    return rdf
end
# ###############################################################################

# ** DataFrames Loaders
function load_joint_eqdf(df_fpath_name::String; svm::Bool=true)
    eq_type = Symbol()
    if occursin(eq_type_dict["full_info"][:dfn], df_fpath_name)
        eq_type = :fi
    elseif occursin(eq_type_dict["misrep"][:dfn], df_fpath_name)
        eq_type = :misrep
    elseif occursin(eq_type_dict["pooling"][:dfn], df_fpath_name)
        eq_type = :pool
    elseif occursin(eq_type_dict["separating"][:dfn], df_fpath_name)
        eq_type = :sep
    end

    cd = get_df_cols_types(eq_type; svm=svm)
    symbol_cols = [x for x in keys(cd) if cd[x] == Symbol]

    df = DataFrame()
    if isfile(df_fpath_name)

        # Replace Symbol DataType by String
        mcd = Dict([x => (x in symbol_cols) ? String : cd[x] for x in keys(cd)])
        df = CSV.read(df_fpath_name, DataFrame; types=mcd)

        # Re-Adjust Columns
        for x in symbol_cols
            df[!, x] = Symbol.(df[:, x])
        end
    else
        println(string(eq_type), " file does not exist!")
    end

    return df
end


# ** DataFrame Parameter Matcher
function fi_param_matcher(jf, ftype::Symbol,
                       rmp::Symbol, df::DataFrame)
    fr = getfield(getfield(jf, ftype), rmp).fr
    cond = [!isnothing(fr) for i in 1:size(df, 1)]
    match = isnothing(fr) ? :skip : :matched

    if cond[1]
        # Type
        cond = .&(cond, Symbol.(df[:, :type]) .== ftype)

        # RMP
        cond = .&(cond, Symbol.(df[:, :rmp]) .== rmp)

        # Bond Contract
        for x in [:m, :c, :p]
            cond = .&(cond, abs.(df[:, x] .- getfield(jf.bc, x)) .< 1e-5)
        end

        # Firm Parameters
        for par in [x for x in param_cols if x != :delta]
            pv = getfield(fr.pm, par)
            if isnan(pv)
                cond = .&(cond, isnan.(df[:, par]))
            else
                cond = .&(cond, abs.(df[:, par] .- pv) .< 1e-5)
            end
        end

        match = sum(cond) == 0 ? :unmatched : match
    end

    return cond, match
end


function misrep_param_matcher(jf, df::DataFrame;
                              st_rmp::Symbol=:rm,
                              rt_rmp::Symbol=:rm,
                              misrep_type::Symbol=:rt)

  colsd = form_cols_lists(Symbol.(names(df));
                            svm=is_svm(jf),
                            eq_type=Symbol(df[1, :eq_type]))
    renamer(ft, x) = (x in [:st_vb, :rt_vb]) ? x : Symbol(ft, :_, x)
    rmp_renamer(ft, x) = (x == :iota) ? Symbol(ft, :_rm_, x) : Symbol(ft, :_nrm_, x)

    sf = getfield(getfield(jf, :st), st_rmp).fr
    rf = getfield(getfield(jf, :rt), rt_rmp).fr
    cond = [.&(!isnothing(sf), !isnothing(rf)) for i in 1:size(df, 1)]
    match = !cond[1] ? :skip : :matched

    if cond[1]
        # RMP
        cond = .&(cond,
                  Symbol.(df[:, renamer(:s,:rmp)]) .==  st_rmp,
                  Symbol.(df[:, renamer(:r,:rmp)]) .==  rt_rmp)

        # Equilibrium Type
        cond = .&(cond, Symbol.(df[:, :eq_type]) .== :misrep)

        # Misrep Type
        cond = .&(cond, Symbol.(df[:, :misrep_type]) .== misrep_type)

        # Bond Contract
        for x in [:m, :c, :p]
            cond = .&(cond, abs.(df[:, x] .- getfield(jf.bc, x)) .< 1e-5)
        end

        # Common Parameters
        for par in [x for x in vcat(colsd[:common_cols], :kappa)]
            pv = getfield(sf.pm, par)
            if isnan(pv)
                cond = .&(cond, isnan.(df[:, par]))
            else
                cond = .&(cond, abs.(df[:, par] .- pv) .< 1e-5)
            end
        end

        # Firm Parameters
        for fr in [sf, rf]
            ft = (fr == sf) ? :s : :r
            for par in [x for x in [:iota, :lambda, :sigmah]]
                pv = getfield(fr.pm, par)
                if isnan(pv)
                    cond = .&(cond, isnan.(df[:, renamer(ft, par)]))
                else
                    cond = .&(cond, abs.(df[:, renamer(ft, par)] .- pv) .< 1e-5)
                    cond = .&(cond, abs.(df[:, rmp_renamer(ft, par)] .- pv) .< 1e-5)
                end
            end
        end

        match = sum(cond) == 0 ? :unmatched : match
    end

    return cond, match
end


function jeq_param_matcher(jf, df::DataFrame, eq_type::Symbol;
                           st_rmp::Symbol=:rm,
                           rt_rmp::Symbol=:rm,
                           misrep_type::Symbol=:rt,
                           jeq_objf::Symbol=:exp_firm_value,
                           misrep_ic_objf::Symbol=:exp_firm_value,
                           misrep_ic_objf_val::Float64=NaN)

  colsd = form_cols_lists(Symbol.(names(df));
                            svm=true, #JointEqStructs.is_svm(jf),
                            eq_type=Symbol(df[1, :eq_type]))

    renamer(ft, x) = (x in [:st_vb, :rt_vb]) ? x : Symbol(ft, :_, x)

    sf = getfield(getfield(jf, :st), st_rmp).fr
    rf = getfield(getfield(jf, :rt), rt_rmp).fr

    cond = [.&(!isnothing(sf), !isnothing(rf)) for i in 1:size(df, 1)]
    match = !cond[1] ? :skip : :matched

    if cond[1]
        # RMP
        cond = .&(cond,
                  Symbol.(df[:, renamer(:s,:rmp)]) .==  st_rmp,
                  Symbol.(df[:, renamer(:r,:rmp)]) .==  rt_rmp)

        # Equilibrium Type
        cond = .&(cond, Symbol.(df[:, :eq_type]) .== eq_type)

        # Misrep Type
        cond = .&(cond, Symbol.(df[:, :misrep_type]) .== misrep_type)

        # Bond Contract
        for x in [:m, :c, :p]
            cond = .&(cond, abs.(df[:, x] .- getfield(jf.bc, x)) .< 1e-5)
        end

        # Common Parameters
        for par in [x for x in vcat(colsd[:common_cols], :kappa)]
            pv = getfield(sf.pm, par)
            if isnan(pv)
                cond = .&(cond, isnan.(df[:, par]))
            else
                cond = .&(cond, abs.(df[:, par] .- pv) .< 1e-5)
            end
        end

        # Firm Parameters
        for fr in [sf, rf]
            ft = (fr == sf) ? :s : :r
            for par in [x for x in [:iota, :lambda, :sigmah]]
                pv = getfield(fr.pm, par)
                if isnan(pv)
                    cond = .&(cond, isnan.(df[:, renamer(ft, par)]))
                else
                    cond = .&(cond, abs.(df[:, renamer(ft, par)] .- pv) .< 1e-5)
                end
            end
        end

        # Joint Eq Specific Columns
        cond = .&(cond, df[:, :jeq_objf] .== jeq_objf)
        cond = .&(cond, df[:, :misrep_ic_objf] .== misrep_ic_objf)
        cond = .&(cond, abs.(df[:, :misrep_ic_objf_val] .- misrep_ic_objf_val) .< 1e-5)

        match = sum(cond) == 0 ? :unmatched : match
    end

    return cond, match
end


# ** DataFrame Results Extractor
function fi_results_extractor(jf, df::DataFrame)
    cond = [false for i in 1:size(df, 1)]

#    eq_type = Symbol(df[1, :eq_type])

    # Store matching results
    matchd = Dict{Symbol, Symbol}()
    for ft in [:st, :rt]
        for rmp in [:rm, :nrm]
            tmp_cond, match = fi_param_matcher(jf, ft, rmp, df)

            cond = cond .| tmp_cond
            matchd[Symbol(ft, :_, rmp)] = match
        end
    end

    return df[cond, :], matchd
end


function misrep_results_extractor(jf, df::DataFrame,
                                  st_fi_rmp::Symbol,
                                  rt_fi_rmp::Symbol)
    cond = [false for i in 1:size(df, 1)]

    if Symbol(df[1, :eq_type]) != :misrep
        println("Wrong dataframe input. Exiting...")
        return
    end

    # Store matching results
    matchd = Dict{Symbol, Dict}(:st_misrep => Dict{Symbol, Any}(),
                                :rt_misrep => Dict{Symbol, Any}())

    # Risky Type Copies Safe Type's Capital Structure
    for rt_rmp in [:rm, :nrm]
        tmp_cond, match = misrep_param_matcher(jf, df;
                                               st_rmp=st_fi_rmp,
                                               rt_rmp=rt_rmp,
                                               misrep_type=:rt)
        cond = cond .| tmp_cond
        matchd[:rt_misrep][Symbol(:st_, st_fi_rmp, :_rt_, rt_rmp)] = match
    end

    # Safe Type Copies Risky Type's Capital Structure
    for st_rmp in [:rm, :nrm]
        tmp_cond, match = misrep_param_matcher(jf, df;
                                               st_rmp=st_rmp,
                                               rt_rmp=rt_fi_rmp,
                                               misrep_type=:st)
        cond = cond .| tmp_cond
        matchd[:st_misrep][Symbol(:st_, st_rmp, :_rt_, rt_fi_rmp)] = match
    end

    return df[cond, :], matchd
end


function jeq_results_extractor(jf, df::DataFrame,
                               jeqid::Dict, eq_type::Symbol;
                               jeq_objf::Symbol=:exp_firm_value)

    cond = [false for i in 1:size(df, 1)]

    if Symbol(df[1, :eq_type]) != eq_type
        println("Wrong dataframe input. Exiting...")
        return
    end

    # Store matching results
    matchd = Dict{Int64, Symbol}()

    # Risky Type Copies Safe Type's Capital Structure
    for row in 1:size(jeqid[:adf], 1)
        tmp_cond, match = jeq_param_matcher(jf, df, eq_type;
                                st_rmp=jeqid[:adf][row, :st_rmp],
                                rt_rmp=jeqid[:adf][row, :rt_rmp],
                                misrep_type=jeqid[:adf][row, :misrep_type],
                                jeq_objf=jeq_objf,
                                misrep_ic_objf=jeqid[:adf][row, :misrep_ic_objf],
                                misrep_ic_objf_val=jeqid[:adf][row, :misrep_ic_objf_val])

        cond = cond .| tmp_cond
        matchd[row] = match
    end

    return df[cond, :], matchd
end


# * Update Results
function type_comparison(df::DataFrame, x, y)
    if typeof(y) == Bool
        return df[:, x] .== y
    elseif typeof(y) == Float64
        cond  = isnan(y) ? isnan.(df[:, x]) : abs.(df[:, x] .- y) .< 1e-5
        return cond
    elseif typeof(y) in [Symbol, String]
        return String.(df[:, x]) .== String(y)
    end
end


function get_jks_id_cols(df::DataFrame)
    eq_type = Symbol(df[1, :eq_type])
    colsd = form_cols_lists(Symbol.(names(df)); eq_type=eq_type)

    id_cols = []
    type_specific_params = [:delta, :lambda, :iota, :sigmah, :rmp]
    if eq_type == :full_info
        id_cols = [x for x in vcat([colsd[:common_cols], [:type],
                                    colsd[:extra_cols], type_specific_params]...)
                   if !(x in [:datetime, :mu_s])]
    elseif eq_type == :misrep
        specific_params = [x for x in colsd[:specific_cols] if
                           any([occursin(string(y), string(x))
                                for y in type_specific_params])]
        fi_cols = [:s_fi_vb, :r_fi_vb]

        #if eq_type == :misrep
        id_cols = [x for x in vcat([colsd[:common_cols], fi_cols,
                                    colsd[:extra_cols],
                                    colsd[:rmp_cols],
                                    specific_params]...)
                   if !(x in [:datetime])]

    elseif eq_type in [:pool, :pooling, :sep, :separating]
        dfct = get_df_cols_types(eq_type; svm=true)
        exc = [:_eq, :leverage, :debt, 
               :yield, :yield_spd, 
               :_MBR, :st_vb, :rt_vb,
               :s_mu_b, :r_mu_b, :date]
        id_cols = [x for x in keys(dfct) if
                   all([!(occursin(String(y), String(x))) for y in exc])]
    else
        println("Error! Equilibrium type not recognized. Exiting...")
        return
    end

    return unique(id_cols)
end


function get_unique_df(jks_fpath::String, eq_type::Symbol;
                       save_filtered_df::Bool=true,
                       del_old_files::Bool=false)
    if eq_type in [:full_info, :fi]
        eqtype = "full_info"
    elseif eq_type in [:misrep, :misrepresentation]
        eqtype = "misrep"
    elseif eq_type in [:pool, :pooling]
        eqtype = "pooling"
    elseif eq_type in [:sep, :separating]
        eqtype = "separating"
    else
        println("Error! Wrong equilibrium type. Exiting...")
        return
    end

    # Get FileName Prefix
    fpath_name_prefix = string(eq_type_dict[eqtype][:dfn])

    # Get List of Files
    files = [x for x in readdir(jks_fpath) if occursin(fpath_name_prefix, x)]
  
    # Main Dataframe
    mainf = string(fpath_name_prefix, ".csv")

    # Load all results files
    df = DataFrame()
    for x in [x for x in files if (x != mainf)]
        # tmp = CSV.read(string(jks_fpath, "/", x))
        tmp = load_joint_eqdf(string(jks_fpath, "/", x); svm=true)
        df = vcat([df, tmp]...)
    end

    # Adjust column names
    for x in [:s, :r]
      col = Symbol(x, :_, x, :t_vb)  
      if col in Symbol.(names(df))
        rename!(df, Dict(col => Symbol(x, :t_vb)))
      end
    end

    # Load Main DataFrame
    if mainf in files
        println("Loading Main Dataframe")
        main = load_joint_eqdf(string(jks_fpath, "/", mainf); svm=true)
        
        # Identify missing columns
        mcols = Symbol.(names(main))
        dfcols = Symbol.(names(df))
        missing_cols = [x for x in dfcols if !(x in mcols)]
      
        # Create missing columns
        for mc in missing_cols
          # identify column type
          ctp = typeof(df[1, mc])
          if ctp  == Float64
            main[!, mc] .= NaN
          else 
            println(string("Error! Need object for column type ", ctp))
            return
          end
        end
        
        # Concatenate dataframes
        println("concatenating...")
        df = vcat([df, main]...)
    end

    # Get list of ID columns
    id_cols = get_jks_id_cols(df)

    # Form DataFrame with Unique Rows
    uqdf = unique(df, id_cols)

    # Filter Original DataFrame to extract latest unique results
    filtered_df = DataFrame()
    for row in 1:size(uqdf, 1)
        tmp_uqdf = uqdf[row, :]
        rows = .&([type_comparison(df, x, tmp_uqdf[x]) for x in id_cols]...)
        tmpdf = df[rows, :]
        filtered_df = vcat([filtered_df, DataFrame(tmpdf[argmax(tmpdf[:, :datetime]), :])]...)
    end

    if del_old_files
        res_files = [x for x in readdir(jks_fpath) if
                     .&(occursin(fpath_name_prefix, x), !occursin("fig", x),
                       # Preserve Old Files
                       [!occursin(string("_", i, ".csv"), x) for i in range(1, stop=20)]...)]

        if !isempty(res_files)
            println("Removing old files...")
            for file in res_files
                rm(string(jks_fpath, "/", file))
            end

            save_filtered_df=true
        end
    end

    if save_filtered_df
        df_fpath_name = string(jks_fpath, "/", fpath_name_prefix, ".csv")

        # Save Filtered Results
        println("Saving filtered results...")
        CSV.write(df_fpath_name, filtered_df)
    end

    return filtered_df
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

    if svm_comb_num > 0
        cvm_comb_num = 0
    end

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
        ftdf = df[Symbol.(df[:, :type]) .== ftype, :]
        cond = .&(cond, size(ftdf, 1) == 1)
        if !cond
            println(string("Error! Type", ft, "sub-dataframe has more than one row."))
            println("Exiting...")
        end

        ft_rmp = Symbol.(ftdf[1, :rmp])
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
