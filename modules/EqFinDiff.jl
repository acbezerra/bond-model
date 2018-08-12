module EqFinDiff

# I am using functions from the ANalytic Functions module
# rfbond_price

# User-defined package:
using AnalyticFunctions: rdisc, rdisc_pvs, rfbond_price,
			 zhi_eq, get_rgrow, get_cvm_vb, get_param
using BondPrInterp: get_cvm_debt_price, get_pv_rfdebt


fun_rel_path  = "../functions/"
modname = "EqFinDiff"
include(joinpath(fun_rel_path, modname, "eq_get_methods.jl")) 
include(joinpath(fun_rel_path, modname, "eq_fin_diff_funs.jl")) 
include(joinpath(fun_rel_path, modname, "tmp_file.jl")) 

end
