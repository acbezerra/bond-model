module_path = "/home/artur/BondPricing/Julia/modules/"
modlpath = joinpath(module_path, "AnalyticFunctions")
if !(modlpath in LOAD_PATH)
    push!(LOAD_PATH, modlpath)
end


# Interpolate Primitive SVM Bond Pricing Functions
module BondPrInterp

using Distributed  # @spawn macro
using Distributions
using Interpolations

# User-defined package:
using AnalyticFunctions: rdisc, rdisc_pvs, rfbond_price,
                         cvm_bond_price, rf_debt,
                         get_rgrow, get_cvm_vb, get_param,
                         psi_v_td, dv_psi_v_td,
                         cvm_F, cvm_G, zhi_vb,
                         on_default_payoff, no_vol_shock_cf_pv

include("svm_interp_funs.jl")
include("svm_interp_fin_diff_funs.jl")
include("svm_interp_mat_funs.jl")
include("bondpr_funs.jl") 
include("debtpr_funs.jl") 

end


# fun_rel_path  = "../functions/"
# modname = "BondPrInterp"
# include(joinpath(fun_rel_path, modname, "svm_interp_funs.jl")) 
# include(joinpath(fun_rel_path, modname, "svm_interp_fin_diff_funs.jl")) 
# include(joinpath(fun_rel_path, modname, "svm_interp_mat_funs.jl")) 
# include(joinpath(fun_rel_path, modname, "bondpr_funs.jl")) 
# include(joinpath(fun_rel_path, modname, "debtpr_funs.jl")) 

