module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["ModelObj", "AnalyticFunctions",
            "EqFinDiff", "Batch"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


module FullInfoEq 

using Distributed
using Dierckx

using Parameters
using Printf
using DataFrames
using CSV

using ModelObj: get_obj_model

using AnalyticFunctions: get_cvm_vb,
                         get_param

using EqFinDiff: eq_fd

using Batch: interp_values


include("_fi_auxiliary_functions.jl")
include("_fi_optimal_vb.jl")
include("_fi_optimal_bond_measure.jl")

end
