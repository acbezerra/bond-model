

function get_ep_ks_dir(jf)
    return string("m_", @sprintf("%.1f", jf.jks.m), 
                  "__c_", @sprintf("%.2f", jf.jks.c), 
                  "__p_", @sprintf("%.2f", jf.jks.p))
end


function get_ep_res_path(jf; ep_dir::String=ep_dir)
    main_dir_path = jf.svm_bt.mi.main_dir_path
    res_path = string(main_dir_path, "/", jf.svm_bt.dfn.res_dir)
    return string(res_path, "/", ep_dir)
end


function get_ep_ks_path(jf; ep_dir::String=ep_dir)
    ep_path = get_ep_res_path(jf; ep_dir=ep_dir)
    ep_ks_dir = get_ep_ks_dir(jf)
    return string(ep_path, "/", ep_ks_dir) 
end


function create_ep_results_dir(jf; ep_dir::String=ep_dir)
    ep_path = get_ep_res_path(jf; ep_dir=ep_dir)
    ep_ks_path = get_ep_ks_path(jf) 

    # Check if folders exist. If not, create them.
    for xdir in [ep_path, ep_ks_path]
        if !isdir(xdir)
            mkdir(xdir)
        end
    end
end


function get_ep_results_path(jf; ep_dir::String=ep_dir)
    create_ep_results_dir(jf; ep_dir=ep_dir)
    return get_ep_ks_path(jf)
end


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


function exists_ep_df(fpath, dfname)
    return string(dfname, ".csv") in readdir(fpath)
end
# function exists_ep_df(jf, dfname; ep_dir::String=ep_dir)
#     #    ep_ks_fi_path = get_ep_full_info_path(jf; ep_dir=ep_dir, fi_dir=fi_dir)
#     ep_ks_path = get_ep_results_path(jf; ep_dir=ep_dir)
#     return string(dfname, ".csv") in readdir(ep_ks_path)
# end