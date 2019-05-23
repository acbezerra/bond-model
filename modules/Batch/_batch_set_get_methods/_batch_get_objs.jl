
# ############## Get Batch & SVM Objects ##############
function get_bt(; model::String="svm", comb_num::Integer=0,
                m::Float64=0., m_comb_num::Integer=0,
                bt=nothing, display_msgs=true)
    # Create Model Object'
    if bt==nothing
        bt = BatchObj(;model=model)
    end
    
    # Set Combination
    bt = set_par_dict(bt; comb_num=comb_num,
                      m=m, m_comb_num=m_comb_num,
                      display_msgs=display_msgs)
    
    # Set Directories & Paths (Main, Batch, Maturity, Combination)
    bt = set_comb_res_paths(bt)

    # Set Directories & Paths (Main, Batch, Maturity, Combination)
    bt = set_comb_res_paths(bt)

    return bt
end


function get_bt_cvm(;comb_num::Integer=0,
                    m::Float64=0., m_comb_num::Integer=0,
                    bt=nothing, display_msgs::Bool=true)

    bt = get_bt(;model="cvm",
                comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                bt=bt, display_msgs=display_msgs)
    cvm = firm_constructor(bt.mi._svm_dict; model="cvm")
    
    return bt, cvm
end    


function get_bt_svm(;comb_num::Integer=0,
                    m::Float64=0., m_comb_num::Integer=0,
                    bt=nothing, display_msgs::Bool=true,
                    compute_surfs::Bool=false)
    bt = get_bt(;model="svm",
                comb_num=comb_num,
                m=m, m_comb_num=m_comb_num,
                bt=bt, display_msgs=display_msgs)

    # Create Directories
    mk_comb_res_dirs(bt)
    
    local batch_obj_exists = false
    try
        # batch_obj_exists = check_batch_results(bt)
        # batch_obj_exists = check_batch_file(bt)
        println("Loading SVM object...")
        svm = load_bpr_surfs(bt)
        
        # Check if Bond Pricing Surfaces exist
        surfs_exist = check_bpr_surfaces(svm)
        batch_obj_exists = true
    catch
        # println("Unable to load Batch object file. Recomputing... ")
        println("Unable to locate Batch object file. Recomputing...")
        batch_obj_exists = false
    end

    
    println(string("Batch object exists: ", batch_obj_exists))
    # Check if results exist
    if !batch_obj_exists
    #     # Load SVM object
    #     println("Loading SVM object...")
    #     # svm = load_batch_results(bt)["svm"]
    #     svm = load_bpr_surfs(bt)
        
    #     # Check if Bond Pricing Surfaces exist
    #     surfs_exist = check_bpr_surfaces(svm)
    # else
        # Create SVM object
        svm = firm_constructor(bt.mi._svm_dict; model="svm")
        surfs_exist = false
    end

    if !surfs_exist | compute_surfs
        # Compute Bond Price Inputs
        svm = @time bpr_surfs(svm)

        # Save results
        save_svm_surfaces(svm, bt.mi.comb_res_path)
    end

    # Interpolate Bond Pricing Surfaces
    println("Interpolating bond pricing surfaces...")
    svm = @time bpr_interp(svm)
    
    return bt, svm
end


function get_bt_mobj(;model::String="svm",
                     comb_num::Integer=0,
                     m::Float64=0., m_comb_num::Integer=0,
                     bt=nothing,
                     display_msgs::Bool=true,
                     compute_surfs::Bool=false)

    if model == "svm"
        return get_bt_svm(;comb_num=comb_num,
                          m=m, m_comb_num=m_comb_num,
                          bt=nothing, display_msgs=display_msgs,
                          compute_surfs=compute_surfs)
    elseif model == "cvm"
        return get_bt_cvm(;comb_num=comb_num,
                          m=m, m_comb_num=m_comb_num,
                          bt=nothing, display_msgs=display_msgs)
    else
        println("Please choose model cvm or svm. Exiting...")
        return
    end
end
