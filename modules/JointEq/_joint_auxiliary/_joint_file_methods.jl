
function get_jeq_mus_fname(jks_fpath::String; eq_type::String="pooling", mu_s::Float64=NaN)
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end

    musn="_tmp"
    if !isnan(mu_s)
        musn = string(jeq_comb_folder_dict[:mu_s][1], 
                      str_format_fun(jeq_comb_folder_dict[:mu_s][2], mu_s))
    end
    return string(jks_fpath, "/", eq_type_dict[eq_type][:fn_prefix], musn, ".csv")
end


function get_jeq_contour_fname(jks_fpath::String, comb_num::Int64; eq_type::String="full_info")
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end

    return string(jks_fpath, "/", eq_type_dict[eq_type][:fn_prefix], "_", comb_num, ".csv")
end


function collect_joint_eq_files(jks_fpath::String; eq_type::String="full_info")
    if !(eq_type in keys(eq_type_dict))
        println("Unrecognized equilibrium type.")
        println("Please enter 'full_info', 'misrep', 'pooling' or 'separating'.")
        println("Exiting...")
        return
    end
    
    return [CSV.read(string(jks_fpath, "/", x)) for x in readdir(jks_fpath) 
                if occursin(string(eq_type_dict[eq_type][:fn_prefix], "_"), x)]
end


function process_results(jks_fpath::String, eq_type::String)
    if !(eq_type in keys(eq_type_dict))
        println("Please set equilibrium_type to 'full_info', 'misrep', 'pooling' or 'separating'. Exiting...")
        return
    end
    
    df_l = [CSV.read(string(jks_fpath, "/", x)) for x in readdir(jks_fpath) 
            if occursin(string(eq_type_dict[eq_type][:fn_prefix], "_"), x)]
    df = vcat(df_l...)
    if eq_type == "full_info"
        df = df[nonunique(df[[:iota, :sigmah]]) .== false, :]
        df = vcat(sort(df[isnan.(df[:sigmah]), :], [:iota]),
                  sort(df[df[:iota] .== .0, :], [:sigmah]))   
    else
        df = vcat(sort(df[isnan.(df[:r_sigmah]), :], [:r_iota]),
                  sort(df[df[:r_iota] .== .0, :], [:r_sigmah]))
    end            
end



# ###########################################################################
# TRASH #####################################################################
# ###########################################################################
# function create_ep_full_info_path(jf; ep_dir::String=ep_dir,
#                                   fi_dir::String=fi_dir)
#     ep_ks_path = get_ep_results_path(jf; ep_dir=ep_dir)
#     ep_ks_fi_path = string(ep_ks_path, "/", fi_dir)
#     if !isdir(ep_ks_fi_path)
#         mkdir(ep_ks_fi_path)
#     end
# end


# function get_ep_full_info_path(jf; ep_dir::String=ep_dir,
#                                fi_dir::String=fi_dir)
#     create_ep_full_info_path(jf; ep_dir=ep_dir, fi_dir=fi_dir)
#     return ep_ks_fi_path = string(get_ep_results_path(jf), "/", fi_dir)
# end


# function exists_ep_df(jf, dfname; ep_dir::String=ep_dir)
#     #    ep_ks_fi_path = get_ep_full_info_path(jf; ep_dir=ep_dir, fi_dir=fi_dir)
#     ep_ks_path = get_ep_results_path(jf; ep_dir=ep_dir)
#     return string(dfname, ".csv") in readdir(ep_ks_path)
# end


# function get_ep_ks_dir(jf)
#     return string(jeq_comb_folder_dict[:m][1], str_format_fun(jeq_comb_folder_dict[:m][2], jf.jks.m), 
#                   jeq_comb_folder_dict[:c][1], str_format_fun(jeq_comb_folder_dict[:c][2], jf.jks.c), 
#                   jeq_comb_folder_dict[:p][1], str_format_fun(jeq_comb_folder_dict[:p][2], jf.jks.p))
# end


# function get_ep_res_path(jf; ep_dir::String=ep_dir)
#     main_dir_path = jf.svm_bt.mi.main_dir_path
#     res_path = string(main_dir_path, "/", jf.svm_bt.dfn.res_dir)
#     return string(res_path, "/", ep_dir)
# end


# function get_ep_ks_path(jf; ep_dir::String=ep_dir)
#     ep_path = get_ep_res_path(jf; ep_dir=ep_dir)
#     ep_ks_dir = get_ep_ks_dir(jf)
#     return string(ep_path, "/", ep_ks_dir) 
# end


# function create_ep_results_dir(jf; ep_dir::String=ep_dir)
#     ep_path = get_ep_res_path(jf; ep_dir=ep_dir)
#     ep_ks_path = get_ep_ks_path(jf) 

#     # Check if folders exist. If not, create them.
#     for xdir in [ep_path, ep_ks_path]
#         if !isdir(xdir)
#             mkdir(xdir)
#         end
#     end
# end


# function get_ep_results_path(jf; ep_dir::String=ep_dir)
#     create_ep_results_dir(jf; ep_dir=ep_dir)
#     return get_ep_ks_path(jf)
# end


# function exists_ep_df(fpath, dfname)
#     return string(dfname, ".csv") in readdir(fpath)
# end

