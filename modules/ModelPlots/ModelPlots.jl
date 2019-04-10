module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["Batch"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


module ModelPlots

using Distributed
using Dierckx

using Parameters
using Printf
using DataFrames
using CSV

# Plots
using PyPlot
using Seaborn
using LaTeXStrings

# User Defined Modules
using Batch: BatchStruct,
             main_dir,
             mat_dir_prefix,
    form_main_dir_path,
    load_cvm_opt_results_df,
    load_svm_opt_results_df

plots_dir="Plots"
graph_dir="HeatSurf"

xylabels = Dict{Symbol, Array{String,1}}(:m => ["m", ".0f"],
                                         :xi => ["\\xi", ".0f"],
                                         :kappa => ["\\kappa^{{EP}} \\, (b.p.)", ".0f"],
                                         :lambda => ["\\lambda", ".2f"],
                                         :sigmal => ["\\sigma_l", ".2f"],
                                         :sigmah => ["\\sigma_h", ".2f"])

zlabels = Dict{Symbol, Array{String,1}}(:c => ["Coupon", "%.2f"],
                                       :p => ["Principal", "%.2f"],
                                       :vb => ["VB", "%.1f"],
                                       :debt => ["Debt", "%.1f"],
                                       :equity => ["Equity", "%.1f"],
                                       :firm_value => ["Debt + Equity", "%1d"],
                                       :leverage => ["Leverage", "%1d"],
                                       :ROE => ["ROE", "%1d"])


@with_kw mutable struct PlotStruct
    # Batch Obj
    bt

    # Results DataFrame
   
    cvm_data::DataFrame
    svm_data::DataFrame

    # Surface Slice
    _svm_surf::DataFrame

    # Labels
    xylabels::Dict{Symbol, Array{String,1}}
    zlabels::Dict{Symbol,Array{String,1}}
end






function PlotsObj(bt;
                  load_cvm_data::Bool=true,
                  load_svm_data::Bool=true,
                  svm_m::Float64=NaN,
                  xylabels::Dict{Symbol, Array{String, 1}}=xylabels,
                  zlabels::Dict{Symbol,Array{String,1}}=zlabels)

    cvm_data = DataFrame()
    if load_cvm_data
        try
            cvm_data = load_cvm_opt_results_df()
        catch
            println("Unable to load CVM data")
        end
    end
    
    svm_data = DataFrame()
    if load_svm_data
        try
            svm_data = load_svm_opt_results_df(bt; m=svm_m)
        catch
            println("Unable to load SVM data")
        end
    end

    return PlotStruct(bt,
                      cvm_data,
                      svm_data,
                      DataFrame(),
                      xylabels, zlabels)
end


function set_cvm_data(pt::PlotStruct, df::DataFrame)
    if isempty(pt.cvm_data)
        println("Setting the Data")
    else
        println("Updating the Data")
    end

    pt.cvm_data=df
end


function plots_form_combinations(bt, fig_name_vars::Array{Symbol,1})
    value_lists = [bt.bp._param_values_dict[x] for x in fig_name_vars]

    # Remaining Params
    # rparams = [x for x in bt._params_order if !(x in fig_name_vars)]

    combs  = [[x1, x2, x3, x4, x5] for 
              x1 in bt.bp._param_values_dict[fig_name_vars[1]],
              x2 in bt.bp._param_values_dict[fig_name_vars[2]],
              x3 in bt.bp._param_values_dict[fig_name_vars[3]],
              x4 in bt.bp._param_values_dict[fig_name_vars[4]],
              x5 in bt.bp._param_values_dict[fig_name_vars[5]]]

    # combinations  = hcat(combs...)'

    return combs
end


function plot_svmdf_slicer(pt, fixed_params::Dict{Symbol, Float64})
    locs = [pt.svm_data[:, x] .== fixed_params[x] for x in keys(fixed_params)]
    srows = sum(hcat(locs...), dims=2) .== size(hcat(locs...), 2)
    return pt.svm_data[vcat(srows...), :]
end


include("_svm_plot_methods.jl")
include("_cvm_vs_svm_plot_methods.jl")

end
