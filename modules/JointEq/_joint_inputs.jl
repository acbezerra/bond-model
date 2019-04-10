


empty_jep = JointEqParams(vcat(fill(NaN,3), 
                               FirmSpecificParams(fill(NaN,3)...),
                               FirmSpecificParams(fill(NaN,3)...),
                               FirmCommonParams(fill(NaN, 7)...))...)

# Directories and File Names
ep_dir = "EP"
fi_dir = "FI"
fidf_name = "fidf"
misrepdf_name = "misrepdf"
pooldf_name = "pooldf"
sepdf_name = "sepdf"


# DataFrames Column Types
fidf_col_types = vcat(fill(Float64, 9), [Bool], fill(Float64, 17))
# Misrepresentation, Pool and Separating DFs Column Types
mps_col_types = vcat(Bool, fill(Float64, 11), Bool,
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


