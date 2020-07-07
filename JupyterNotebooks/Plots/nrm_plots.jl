# %% codecell
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
using PyPlot
# using PyCall
using Seaborn
using LaTeXStrings

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
# push!(LOAD_PATH, module_path)
modls = ["ModelObj", "AnalyticFunctions", "BondPrInterp",
         "EqFinDiff", "Batch", "ModelPlots",
         "FullInfoEq", "JointEqStructs", "JointEqStructs"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

ENV["LINES"] = 750
ENV["COLUMNS"] = 1000



# %% codecell
using Distributions
modls = ["Batch", "ModelObj", "AnalyticFunctions",
         "BondPrInterp", "EqFinDiff",
         "JointEqStructs", "FullInfoEq", "JointEq", "ModelPlots"]
