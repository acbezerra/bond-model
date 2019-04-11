
function diag_df(bt, i::Integer)
    combdf = DataFrame(bt.bp.df[i, :])
    
    bt = Batch.set_par_dict(bt; comb_num=combdf[:comb_num][1], display_msgs=false)
    bt = Batch.set_comb_res_paths(bt)
    
    combdf = hcat(combdf, DataFrame(folder_name = bt.mi.comb_res_dir,
                                    folder_created = isdir(bt.mi.comb_res_path),
                                    batch_obj=false, 
                                    modified = " ", 
                                    bond_surfs = false, 
                                    count = 0))

    # Count solutions:
    if combdf[:folder_created][1]
        dfs_list = [x for x in readdir(bt.mi.comb_res_path) if 
                    occursin(bt.dfn.eq_fd_cp_fn_prefix, x)]
        combdf[:count] = size(dfs_list, 1)
        combdf[:batch_obj] = sum([occursin(Batch.batch_file_name, x) for 
                                  x in readdir(bt.mi.comb_res_path)]) > 0

        if combdf[:batch_obj][1]
            batch_obj_file = string(bt.mi.comb_res_path, "/",
                                    Batch.batch_file_name, ".jld")
            combdf[:modified] = string(Dates.unix2datetime(stat(batch_obj_file).mtime))
        end
        combdf[:bond_surfs] = sum([occursin("bond_surfaces", x) for 
                                   x in readdir(bt.mi.comb_res_path)]) > 0
    end
    
    return combdf
end


function diagnosis(bt)
    # combDFs = fetch(@spawn [diag_df(bt, comb) for comb in 
    #                         range(1, stop=size(bt._params_combs[:], 1))])

    combDFs = fetch(@spawn [diag_df(bt, i) for i in 1:nrow(bt.bp.df)])
    cols1 = [:comb_num, :m, :m_comb_num, :folder_created, :batch_obj,
             :modified, :bond_surfs, :count, :folder_name]
    cols2 = [x for x in names(combDFs[1]) if !(x in cols1)]
    return sort!(vcat(combDFs...)[vcat(cols1, cols2)], [:m, :m_comb_num])
end