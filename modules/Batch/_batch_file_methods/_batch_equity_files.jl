
function conditionally_load_df(df, path_fname::String,
                               types::Array{DataType,1}=vcat(fill(Float64, 8), [Int64, Int64]))
    try 
        return eval(df)
    catch
        return CSV.read(path_fname; types=types)
    end
end


function save_eq_results(bt, res, df::DataFrame, eq_fd_cp_path_fname::String)
    sol = sort!(vcat(df[bt.dfc.dfcols], res[bt.dfc.dfcols]), [:p])

    # Identify solution
    sol[:opt] = .&(abs.(sol[:abs_debt_diff] .- 
                         minimum(sol[:abs_debt_diff])) .< 1e-4, 
                        abs.(sol[:eq_deriv]) .< 1e-4)
    
    CSV.write(eq_fd_cp_path_fname, sol)
end


# function load_eq_results(bt, path_file_name; 
#                          col_types=vcat(fill(Float64, 16), [Bool], 
#                                         fill(Float64, 14), [Bool]))
#     dfs_list = [x for x in readdir(path_file_name) if 
#                 occursin(bt.dfn.eq_fd_cp_fn_prefix, x)]

#     LL = fetch(@spawn [CSV.read(string(path_file_name, "/", df), 
#                    types=col_types) for df in dfs_list])

#     return sort!(vcat(LL...), [:c])
# end


function load_eq_results(bt, svm, dfn; use_all_eqdf::Bool=true)
    coltypes = vcat(fill(Float64, 16), [Bool], fill(Float64, 14))

    if .&(use_all_eqdf,  string("all_", dfn) in readdir(bt.mi.comb_res_path))
        all_eqdf = CSV.read(string(bt.mi.comb_res_path, "/all_", dfn))
        if size(all_eqdf, 2) == size(coltypes, 1)
            all_eqdf = CSV.read(string(bt.mi.comb_res_path, "/all_", dfn),
                            types=coltypes)
        end
        
       # Compute candidate vb for each p value
        eqdf_final = Batch.eq_fd_processing(bt, svm, all_eqdf)

        # Compute (p, vb) solution for c value
        res = Batch.eq_fd_sol(bt, svm, eqdf_final)

        # Save results
        # Batch.save_eq_results(bt, res, eqdf_final,
        #                      string(bt.mi.comb_res_path, "/", dfn))
        sol = sort!(vcat(eqdf_final[bt.dfc.dfcols], res[bt.dfc.dfcols]), [:p])

        # Identify solution
        sol[:opt] = .&(abs.(sol[:abs_debt_diff] .- 
                            minimum(sol[:abs_debt_diff])) .< 1e-4, 
                       abs.(sol[:eq_deriv]) .< 1e-4)

        # Drop duplicates
        unique!(sol, [:p])
    
        CSV.write(string(bt.mi.comb_res_path, "/", dfn), sol)

        return sol
    elseif !use_all_eqdf   
        sol = CSV.read(string(bt.mi.comb_res_path, "/", dfn); 
                       types=vcat(coltypes, [Bool]))

        return unique!(sol, [:p, :vb])
    else
        return DataFrame()
    end
end

