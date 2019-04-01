module_path = "/home/artur/BondPricing/Julia/modules/"
modlpath = joinpath(module_path, "ModelObj")
if !(modlpath in LOAD_PATH)
    push!(LOAD_PATH, modlpath)
end

module AnalyticFunctions

using Distributions
using Interpolations

# User-defined package:
using ModelObj: extract_param

include("cvm_funs.jl")
include("svm_analytic_funs.jl")
include("af_get_methods.jl")

end
