module_path = "/home/artur/BondPricing/Julia/modules/"
modnames = ["ModelObj", "AnalyticFunctions",
            "BondPrInterp", "EqFinDiff",
            "Batch", "FullInfoEq"]
for modl in modnames
    if !(joinpath(module_path, modl) in LOAD_PATH)
        push!(LOAD_PATH, joinpath(module_path, modl))
    end
end


module JointEq 

using Distributed
using Dierckx

using Parameters
using Printf
using DataFrames
using CSV

using ModelObj: set_opt_k_struct,
                Firm, grid_creator,
                get_obj_model

using AnalyticFunctions: get_cvm_vb,
                         get_param,
                         rfbond_price,
                         on_default_payoff

using BondPrInterp: get_pv_rfdebt,
                    get_cvm_bond_price,
                    get_cvm_debt_price,
                    get_svm_bond_price,
                    get_svm_debt_price
    
using EqFinDiff: get_cvm_eq,
                 get_eq_Vmax,
                 eq_set_k_struct,
                 eq_fd_core,
                 eq_fd_export_results,
                 eq_fd

using Batch: BatchObj,
             get_bt,
             get_bt_svm,
             get_bt_mobj,
             get_batch_comb_num,
             load_cvm_opt_results_df,
             load_svm_opt_results_df,
             opt_k_struct_df_name,
             opt_k_struct_df_coltypes,
             BatchStruct,
             interp_values,
             svm_param_values_dict,
             common_params,
             _params_order

using FullInfoEq: set_full_information_vb!, find_optimal_bond_measure

# Structs, Objects and Constructor Methods #################
include("_joint_objects/_joint_structs.jl")

# Inputs -> in this order (need to define structs first)
include("_joint_inputs.jl")


include("_joint_objects/_joint_ep_constructor.jl")
include("_joint_objects/_joint_constructors.jl")
include("_joint_objects/_joint_k_struct_funs.jl")
# ##########################################################

# File and DataFrame Methods ###############################
include("_joint_auxiliary/_joint_functions.jl")
include("_joint_auxiliary/_joint_file_methods.jl")
include("_joint_auxiliary/_joint_dataframe_methods.jl")
# ##########################################################

include("_joint_pricing.jl")

# Equilibrium Methods -> Equity, Vb, mu_b Functions ########
include("_joint_equilibrium/_joint_eq_fin_diff.jl")
include("_joint_equilibrium/_joint_optimal_vb.jl")
include("_joint_equilibrium/_joint_optimal_bond_measure.jl")
# ##########################################################

end
