
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
modls = ["AnalyticFunctions",  "ModelObj", "BondPrInterp",
         "EqFinDiff", "Batch", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end


#bt = Batch.BatchObj()
#pt = ModelPlots.PlotsObj(bt; firm_obj_fun=:firm_value, cvm_m =1., svm_m =1.)
#sbt = Batch.get_bt(; model="svm", m=1., m_comb_num=1)
#cbt = Batch.get_bt(; model="cvm", m=1., m_comb_num=1)


# Market Illiquidity - Transaction Costs ############################
# Parameters
cvmdict = Dict(:mu_b => [1.], 
               :m => [1.], 
               :xi => [1.], 
               :sigmal => [.15], 
               :kappa => [25., 50.] .* 1e-4, 
               :iota => [5. * 1e-4])
svmdict = deepcopy(cvmdict)
svmdict[:iota] = [0.]
svmdict[:kappa] = [10., 25.] .* 1e-4
svmdict[:lambda] = [.1]

# cvm_combs = Batch.get_batch_comb_numbers(cbt, cvmdict)[:comb_num]
# svm_combs = Batch.get_batch_comb_numbers(sbt, svmdict)[:comb_num]

cvmdf, svmdf, xvar = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                                firm_obj_fun=:firm_value)

dvar = :kappa
for yvar in [:firm_value, :MBR]
    fig = ModelPlots.cvm_vs_svm_plotfun(cvmdf, svmdf,
                                        xvar, yvar, dvar;
                                        save_fig=true)
end
# ###################################################################


# Risk-Management Costs #############################################
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                               :m  => [1.0],
                               :gross_delta => [0.02],
                               :kappa  => [0.0025],
                               :mu_b => [1.0],
                               :xi => [1.0],
                               :iota => [0.0005,  0.00075])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.1]
svmdict[:kappa] = [.001]
svmdict[:iota] = [.0]

# cvm_combs = Batch.get_batch_comb_numbers(cbt, cvmdict)[:comb_num]
# svm_combs = Batch.get_batch_comb_numbers(sbt, svmdict)[:comb_num]

cvmdf, svmdf, xvar = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                                firm_obj_fun=:firm_value)


dvar = :iota
for yvar in [:firm_value, :MBR]
    fig = ModelPlots.cvm_vs_svm_plotfun(cvmdf, svmdf,
                                        xvar, yvar, dvar;
                                        save_fig=true)
end
# ###################################################################


# Volatility Risk Intensity #########################################
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                               :m  => [1.0],
                               :gross_delta => [0.02],
                               :kappa  => [0.0025],
                               :mu_b => [1.0],
                               :xi => [1.0],
                               :iota => [0.0005])
svmdict = deepcopy(cvmdict)
svmdict[:lambda] = [.1, .3]
svmdict[:kappa] = [.001]
svmdict[:iota] = [.0]

# cvm_combs = Batch.get_batch_comb_numbers(cbt, cvmdict)[:comb_num]
# svm_combs = Batch.get_batch_comb_numbers(sbt, svmdict)[:comb_num]


cvmdf, svmdf, xvar = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                                firm_obj_fun=:firm_value)

dvar = :lambda
for yvar in [:firm_value, :MBR]
    fig = ModelPlots.cvm_vs_svm_plotfun(cvmdf, svmdf,
                                        xvar, yvar, dvar;
                                        save_fig=true)
end
# ###################################################################
