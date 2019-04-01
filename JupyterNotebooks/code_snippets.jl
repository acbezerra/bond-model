using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
include(string(module_path, "/", "TestFunctions.jl"))
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff"]
for modl in modls
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end

include(string(module_path, "/", "TestFunctions.jl"))

ENV["LINES"] = 750

bt = Batch.BatchObj()
bt = Batch.set_par_dict(bt, 1)

# Set Directories & Paths (Main, Batch, Maturity, Combination)
bt = Batch.set_comb_res_paths(bt)

# Construct Firm Object
svm = ModelObj.firm_constructor(bt.mi._svm_dict)


mval = unique(diagdf[:m])[1]
dfm1 = sort(diagdf[diagdf[:m] .== mval, :], [:m, :comb_num])

latest_day = maximum(Dates.Day.(DateTime.(diagdf[:, :modified])))
cond = Dates.Day.(DateTime.(dfm1[:, :modified])) .== latest_day
size(dfm1, 1) - sum(cond)



