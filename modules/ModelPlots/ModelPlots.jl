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
    comb_folder_dict,
    form_main_dir_path,
    load_cvm_opt_results_df,
    load_svm_opt_results_df,
    str_format_fun, get_bt,
    form_main_dir_path, main_dir,
    BatchObj, get_bt,
    get_batch_comb_numbers

include("_plots_inputs.jl")

@with_kw mutable struct PlotStruct
    # Batch Obj
    bt

    # Results DataFrame
    cvm_data::DataFrame
    svm_data::DataFrame

    # Firm Objective Function
    obj_fun::Symbol

    # Surface Slice
    _svm_surf::DataFrame

    # Labels
    xylabels::Dict{Symbol, Array{String,1}}
    zlabels::Dict{Symbol,Array{String,1}}
end


function PlotsObj(bt;
                  load_cvm_data::Bool=true,
                  load_svm_data::Bool=true,
                  firm_obj_fun::Symbol=:MBR,
                  cvm_m::Float64=NaN,
                  svm_m::Float64=NaN,
                  xylabels::Dict{Symbol, Array{String, 1}}=xylabels,
                  zlabels::Dict{Symbol,Array{String,1}}=zlabels)

    cvm_data = DataFrame()
    if load_cvm_data
        try
            cvm_data = load_cvm_opt_results_df(; m=cvm_m,
                                               firm_obj_fun=firm_obj_fun)
        catch
            println("Unable to load CVM data")
        end
    end
    
    svm_data = DataFrame()
    if load_svm_data
        try
            svm_data = load_svm_opt_results_df(bt;
                                               firm_obj_fun=firm_obj_fun,
                                               m=svm_m)
        catch
            println("Unable to load SVM data")
        end
    end

    return PlotStruct(bt,
                      cvm_data,
                      svm_data,
                      firm_obj_fun,
                      DataFrame(),
                      xylabels, zlabels)
end


include("_plots_auxiliary_methods.jl")
include("_cvm_plot_methods.jl")
include("_svm_plot_methods.jl")
include("_cvm_vs_svm_plot_methods.jl")

# Risk-Management Policy Choice
include("_RMP/_rmp_auxiliary.jl")
include("_RMP/_rmp_curves.jl")
include("_RMP/_rmp_color_regions.jl")
include("_RMP/_rmp_fi_plot_methods.jl")

# Joint Equilibrium
include("_JEQ/_jeq_auxiliary.jl")
include("_JEQ/_jeq_curves.jl")
include("_JEQ/_jeq_color_regions.jl")
include("_JEQ/_jeq_plot_methods.jl")

# Contour Plots 
include("_Contour/_contour_auxiliary.jl")
include("_Contour/_contour_plot_methods.jl")
include("_Contour/_contour_payoff_functions.jl")

end
