
# Batch Inputs
# Directories:
main_dir = "Julia"
res_dir = "Results"
mat_dir_prefix = "m_"
coupon_dir_prefix = "c_"

batch_file_name = "batch_obj"

common_params = Dict{Symbol, Float64}(:V0 => 100.,
                                      :alpha => .6,
                                      :pi => .27,
                                      :r => .08)


_params_order = [:mu_b, :m, :xi, :kappa, :gross_delta,
                 :lambda, :iota, :sigmal, :sigmah]

# #######################################################
# ################# Combination Folders #################
# #######################################################
# Notice that mu_b is set to 1.
comb_folder_dict = Dict{Symbol,Array{String,1}}(:mu_b => ["mub_", "%.2f"],
                                                :m => ["m_", "%.2f"],
                                                :c => ["c_", "%.2f"],
                                                :xi => ["__xi_", "%.2f"], 
                                                :kappa => ["__kappa_", "%.2f", "_bp"],
                                                :gross_delta => ["__gross_delta_", "%.2f", "_bp"], 
                                                :iota => ["__iota_", "%.2f", "_bp"],
                                                :lambda => ["__lambda_", "%.3f"],
                                                :sigmal => ["__sigmal_", "%.3f"],
                                                :sigmah => ["__sigmah_", "%.3f"])


# #######################################################


# Define Parameter Combinations:
# #######################################################
# ######################### CVM #########################
# #######################################################
# Filename
cvm_sol_fn_prefix = "cvm_res"

# Parameters
cvm_param_values_dict = Dict{Symbol,Array{Float64,1}}(:mu_b => [1.],
                                                      :m => [1.], #, 3., 5., 7., 10.],
                                                      :xi => [1.],
                                                      :kappa => [k * 1e-4 for k in [25, 50, 100, 150]],
                                                      :gross_delta => [.02],
                                                      :lambda => [NaN],
                                                      :iota => [i * 1e-4 for i in [0.0, .5, 1., 1.5, 2.0, 2.5,
                                                                                   5, 7.5, 10, 15, 20, 25, 50]],
                                                      :sigmal => [.15, .2, .25, .30],
                                                      :sigmah => [NaN])

# Coupon Values Grid:
cvm_c_step = .5
cvm_coupon_grid = vcat([.25, .5], range(1, stop=25 + cvm_c_step, step=cvm_c_step))

# #######################################################
# ######################### SVM #########################
# #######################################################
# Julia Paths & Files


# FileNames
# svm_partI_fn_prefix = "svm_data_c"
# svm_partII_fn_prefix = "svm_res_c"
debt_at_par_cp_fn_prefix = string("debt_at_par_", coupon_dir_prefix)
eq_fd_cp_fn_prefix = string("eq_fd_", coupon_dir_prefix)




# Parameters
svm_param_values_dict = Dict{Symbol,Array{Float64,1}}(:mu_b => [1.],
                                                      :m => [1.],     #, 3., 5., 7., 10.],
                                                      :xi => [1.],
                                                      :kappa => [k * 1e-4 for k in [10, 25, 30, 40, 50]],
                                                      :gross_delta => [.02],
                                                      :lambda => [.1, .2, .3, .5, .75],
                                                      :iota => [.0],
                                                      :sigmal => [.15],
                                                      :sigmah => [.2, .225, .25, .275, .30, .35])


# Coupon Values Grid:
svm_c_step = .5
svm_coupon_grid = Array(range(svm_c_step, stop=11., step=svm_c_step))


# ### Set Tolerance Levels for Filtering p values###
# p_interp_fun -> eq_deriv_tol/debt_tol  -> p_tol_search
# eq_deriv_tol = vcat([1, 5] .* 1e-6, 
#                     kron([1e-5, 1e-4, 1e-3, 1e-2],[1, 2.5, 5, 7.5]),
#                     [1, 2.5, 5] .* 1e-2)
# debt_tol = kron([1, 5], [10^i for i in range(0, stop=3)] .* 1e-4)
eq_deriv_tol1 = kron([1, 5], [1e-6, 1e-5, 1e-4])
debt_tol1 =  kron([1, 5], [10^i for i in [0, 1.]] .* 1e-4)
eq_deriv_tol2 = vcat(kron([1e-3, 1e-2],[1, 2.5, 5, 7.5]),
                    [1, 2.5, 5] .* 1e-2)
debt_tol2 = kron([1, 2.5, 5.], [10^i for i in range(0, stop=2)] .* 1e-4)

function toldf_fun(debt_tol::Array{Float64,1}, eq_deriv_tol::Array{Float64,1})
    # Form all possible combinations and sort them so that
    # debt is relaxed first.
    value_lists = [sort(debt_tol), sort(eq_deriv_tol)]
    cols = [:debt_diff, :eq_deriv]

    A = vcat([(i, j) for i=debt_tol, j=eq_deriv_tol]...)
    return sort(DataFrame(map(idx -> getindex.(A, idx),
                              eachindex(first(A))), cols),
                cols)
end
toldf = sort(vcat(toldf_fun(debt_tol1, eq_deriv_tol1),
                  toldf_fun(debt_tol2, eq_deriv_tol2)),
             [:eq_deriv, :debt_diff])



# Filter (P, VB) values -> discard candidates for which
# |(Debt - P)/P| > tol:
pvb_tol_vec = Array([1, 10, 10^2, 10^3]) * 10^-6



# DataFrame Columns
main_params = [:gross_delta, :delta, :iota, :kappa, :lambda, :sigmah, :sigmal]
k_struct_params = [:mu_b, :m, :c, :p, :vb]
fixed_params = [:V0, :r, :alpha , :pi, :xi]
debt_vars = [:debt_diff,
             :abs_debt_diff,
             :debt_per_diff,
             :abs_debt_per_diff]
equity_vars = [:eq_deriv,       
               :eq_deriv_min_val, 
               :eq_min_val,    
               :eq_negative,     
               :eq_vb]
share_values = [:debt, :equity, :firm_value, :leverage, :MBR]
dfcols = vcat(main_params, k_struct_params, debt_vars[1],
              equity_vars, debt_vars[2:end], share_values, fixed_params)



# #######################################################
# ####################### RESULTS #######################
# #######################################################
soldf_name = "soldf"
optdf_name = "optdf"
opt_k_struct_df_name="opt_k_structs"
#comb_opt_k_struct_df_coltypes=vcat(String, fill(Float64, 33), Bool, Float64)
opt_k_struct_df_coltypes=vcat(Int64, Float64, Int64, String, fill(Float64, 32), Bool, Float64)
cvm_opt_k_struct_df_coltypes=vcat([Int64, Float64, Int64, String],
                                  fill(Float64, 32))
