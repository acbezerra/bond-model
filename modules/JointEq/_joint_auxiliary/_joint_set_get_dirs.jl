
# Get Directories ##############################################
function get_jeq_common_params_dir(jf)
    common_param_dir = ""
    for par in common_dir_par_list
        if par in [:kappa, :gross_delta]
            par_string =  string(str_format_fun(jeq_comb_folder_dict[par][2],
                                                1e4 * get_param(jf.sf, par)),
                                 jeq_comb_folder_dict[par][3])
        else
            par_string =  string(str_format_fun(jeq_comb_folder_dict[par][2],
                                                get_param(jf.sf, par)))
        end
        common_param_dir = string(common_param_dir, jeq_comb_folder_dict[par][1], par_string)
    end

    return common_param_dir
end


function get_jks_dir(jf)
    pcr = jf.jks.p/jf.jks.c
    
    return string(jeq_comb_folder_dict[:m][1], str_format_fun(jeq_comb_folder_dict[:m][2], jf.jks.m), 
                  jeq_comb_folder_dict[:pcr][1], str_format_fun(jeq_comb_folder_dict[:pcr][2], pcr))
end


function get_jeq_jks_fpath(jf)
    jeq_fpath = make_jeq_res_fpath()
    common_params_dir = get_jeq_common_params_dir(jf)
    
    common_params_fpath = string(jeq_fpath, "/", common_params_dir)
    jks_dir = get_jks_dir(jf)

    return string(common_params_fpath, "/", jks_dir)
end
# ##############################################################


# Get Paths ####################################################
function get_res_fpath(jf)
    return string(jf.svm_bt.mi.main_dir_path, "/", res_dir)
end


function get_jeq_res_fpath(jf)
    res_dir_path = make_res_fpath(jf)

    return string(res_dir_path, "/", jeq_dir)
end


function get_jeq_common_params_fpath(jf)
    jeq_fpath = make_jeq_res_fpath(jf)
    common_params_dir = get_jeq_common_params_dir(jf)
    
    return string(jeq_fpath, "/", common_params_dir)
end


function get_jeq_jks_fpath(jf)
    common_params_fpath =  get_jeq_common_params_fpath(jf)
    jks_dir = get_jks_dir(jf)

    return string(common_params_fpath, "/", jks_dir)
end
# ##############################################################


# Make Directories #############################################
function make_res_fpath(jf)
    res_dir_path = get_res_fpath(jf)
    if !isdir(res_dir_path)
        mkdir(res_dir_path)
    end

    return res_dir_path
end


function make_jeq_res_fpath(jf)
    jeq_fpath = get_jeq_res_fpath(jf)
    if !isdir(jeq_fpath)
        mkdir(jeq_fpath)
    end

    return jeq_fpath
end


function make_jeq_common_params_fpath(jf)
    common_params_fpath = get_jeq_common_params_fpath(jf)
    if !isdir(common_params_fpath)
        mkdir(common_params_fpath)
    end

    return common_params_fpath
    
end


function make_jeq_jks_fpath(jf)
    common_params_fpath = make_jeq_common_params_fpath(jf)
    jks_fpath = get_jeq_jks_fpath(jf)
    if !isdir(jks_fpath)
        mkdir(jks_fpath)
    end

    return jks_fpath 
end
# ##############################################################
