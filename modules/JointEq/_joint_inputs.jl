


empty_jep = JointEqParams(vcat(fill(NaN,3), 
                               FirmSpecificParams(fill(NaN,3)...),
                               FirmSpecificParams(fill(NaN,3)...),
                               FirmCommonParams(fill(NaN, 7)...))...)

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
