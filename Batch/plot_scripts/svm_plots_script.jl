

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
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end


# SVM Surface & Heatmap Inputs #######################################
z_vars = [:c, :p, :vb, :debt, :equity, :firm_value, :leverage, :ROE]
fig_name_vars = [:mu_b, :m, :xi, :lambda, :sigmal]
combinations = ModelPlots.plots_form_combinations(bt, fig_name_vars)
xy_list = [:sigmah, :kappa]
fig_name_vars = [:mu_b, :m, :xi, :lambda, :sigmal]

graph_dict1 = Dict{Symbol, Any}(:azim => 50., 
                                :heat_reverse_y => false, 
                                :heat_reverse_x => false,
                                :cbaxes => [.925, 0.15, 0.015, 0.675],
                                :axs_wspace => .2)

graph_dict2 = Dict{Symbol, Any}(:azim => -50., 
                                :heat_reverse_y => true, 
                                :heat_reverse_x => true,
                                :cbaxes => [.975, 0.15, 0.015, 0.675],
                                :axs_wspace => .1)

ax1_dist = 8.5
sup_title=.575
plt_cmap="viridis"
seaborn_style="darkgrid"
return_fig=false
# ####################################################################

bt = Batch.BatchObj()
# pt = ModelPlots.PlotsObj(bt; svm_m =1.)

for comb in combinations
    fixed_params = Dict(zip(fig_name_vars, comb))
    pt = ModelPlots.PlotsObj(bt; svm_m =1.)
    pt = ModelPlots.set_svm_surf(pt, fixed_params)
    
    for z in z_vars
        graph_dict = graph_dict1
        if z == :equity
            graph_dict = graph_dict2
        end

        
        ModelPlots.svm_plot_heatmap_surf(pt, xy_list, z, fixed_params,
                                         heat_reverse_x=graph_dict[:heat_reverse_x],
                                         heat_reverse_y=graph_dict[:heat_reverse_y],
                                         elev=25., azim=graph_dict[:azim], #zpad=10.,
                                         ax1_dist=ax1_dist, 
                                         cbaxes=graph_dict[:cbaxes], 
                                         axs_wspace=graph_dict[:axs_wspace], 
                                         sup_title_x=sup_title, 
                                         plt_cmap=plt_cmap,
                                         seaborn_style=seaborn_style,
                                         return_fig=return_fig)
        PyPlot.close()
    end
end
