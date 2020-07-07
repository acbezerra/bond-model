
# Transaction Costs and Measure of safe firms
st_rt_comb[:k_ep] = 25 * 1e-4
st_rt_comb[:k_otc] = 50 * 1e-4


# ** Form Joint Firms Objects
# ** Extract list with all possible combinations
# Set Objective Functions
firm_obj_fun = :firm_value

# Transaction Costs and Volatility Risk Parameters
pardict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :mu_b => [1.0],
                                        :xi => [1.0])

# Risk Type Parameters
st_iota_vec = [.0]
rt_iota_vec = [.0]
iota_vec = unique(vcat(st_iota_vec, rt_iota_vec))
lambda_vec = Batch.svm_param_values_dict[:lambda]
sigmah_vec = Batch.svm_param_values_dict[:sigmah]
mu_s_vec = [.2]

