push!(LOAD_PATH, "/home/artur/BondPricing/Julia/modules/")

using Distributed

include(string(LOAD_PATH[end], "ModelObj.jl"))
include(string(LOAD_PATH[end], "AnalyticFunctions.jl"))
include(string(LOAD_PATH[end], "BondPrInterp.jl"))
include(string(LOAD_PATH[end], "EqFinDiff.jl"))

pardict = Dict("V0" => 100, 
               "vbl" => 60, 
               "vbh" => 62, 
                "m" => 1.,
                "c" => 6.39,
                "p" => 61.68,
                "sigmal" => .25,
                "sigmah" => .25,
                "r" => .08,
                "gross_delta" => .02,
                "iota" => .0,
                "xi" => 1.,
                "kappa" => .015,
                "alpha" => .6,
                "pi" => .27,
                "lambda" => .3)

svm = ModelObj.firm_constructor(pardict)

f0fd, f1fd, f2fd, f3fd = @time BondPrInterp.fd_interp(svm, vtN=50, vbhlN=10)
