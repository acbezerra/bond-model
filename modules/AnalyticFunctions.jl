# analytic functions 2
module AnalyticFunctions

using Distributions
using Interpolations


fun_rel_path  = "../functions/"
modname = "AnalyticFunctions"

include(joinpath(fun_rel_path, modname, "cvm_funs.jl")) 
include(joinpath(fun_rel_path, modname, "svm_analytic_funs.jl"))
include(joinpath(fun_rel_path, modname, "af_get_methods.jl"))

end


