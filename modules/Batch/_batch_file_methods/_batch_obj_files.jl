
function save_svm_surfaces(svm, file_path)
    bond_path_fname = string(file_path, "/bond_surfaces", ".jld")
    bs = Dict("vtgrid" => svm.bs.vtgrid,
              "ttmgrid" => svm.bs.ttmgrid,
              "vbhlgrid" => svm.bs.vbhlgrid,
              "f11_surf" => svm.bs.f11_surf, 
              "f12_surf" => svm.bs.f12_surf, 
              "f13_surf" => svm.bs.f13_surf, 
              "f21_surf" => svm.bs.f21_surf,
              "f22_surf" => svm.bs.f22_surf)
    # save(bond_path_fname, "bs", bs)
    jldopen(bond_path_fname, "w") do file
        write(file, "bs", bs)
    end
end


function save_batch_results(bt; svm=nothing, fname::String=batch_file_name)
    bt = mk_comb_res_dirs(bt)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")
    

    if !(svm==nothing)
#        save(path_fname, "bt", bt, "svm", svm)
        jldopen(path_fname, "w") do file
            write(file, "bt", bt)
            write(file, "svm", svm)
        end

        save_svm_surfaces(svm, bt.mi.comb_res_path)
    else
        # save(path_fname, "bt", bt)
        jldopen(path_fname, "w") do file
            write(file, "bt", bt) 
        end
    end
end


function check_batch_file(bt; fname::String=batch_file_name)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")

    return .&(isfile(path_fname), stat(path_fname).size > 0) 
end


function load_batch_results(bt; fname::String=batch_file_name)
    path_fname = string(bt.mi.comb_res_path, "/", fname, ".jld")

    exist_cond = false
    if check_batch_file(bt; fname=batch_file_name)
        try
            res = load(path_fname)
            # Check if both bt and svm objects are present:
            exist_cond = .&([x in keys(res) for x in ["bt", "svm"]]...)

            if !exist_cond
                println("Batch or SVM object missing. ")
            end
        catch
            println("Unable to load batch object file.")
        end
    else
        println("Batch object file not found.")
    end

    if exist_cond
        return load(path_fname)
    end
end


function check_batch_results(bt; fname::String=batch_file_name)
    # if batch object exists, check that surfaces have
    # been computed.
    svm_exists = false

    if check_batch_file(bt; fname=fname)
        # Load Batch Object:
        res = load_batch_results(bt)

        # Check that SVM object is in results object:
        return ("svm" in keys(res))
    end

    return svm_exists
end


function check_bpr_surfaces(svm)
    # Check that SVM object contains bond pricing surfaces:
    return !(any(isnan, svm.bs.f11_surf) || 
             any(isnan, svm.bs.f12_surf) ||
             any(isnan, svm.bs.f13_surf) ||
             any(isnan, svm.bs.f21_surf) || 
             any(isnan, svm.bs.f22_surf))
end


function load_bpr_surfs(bt; batch_file_name::String=batch_file_name, interp_surfs::Bool=true)
    try check_batch_results(bt; fname=batch_file_name)
        svm = load_batch_results(bt, fname=batch_file_name)["svm"]
    catch
        println("Batch Object is incompatible! Loading surfaces only instead.")
        svm = firm_constructor(bt.mi._svm_dict)
        svmbs = load(string(bt.mi.comb_res_path, "/bond_surfaces.jld"))["bs"]
        svm.bs.vtgrid = svmbs["vtgrid"]
        svm.bs.ttmgrid = svmbs["ttmgrid"]
        svm.bs.vbhlgrid = svmbs["vbhlgrid"]
        svm.bs.f11_surf = svmbs["f11_surf"]
        svm.bs.f12_surf = svmbs["f12_surf"]
        svm.bs.f13_surf = svmbs["f13_surf"]
        svm.bs.f21_surf = svmbs["f21_surf"]
        svm.bs.f22_surf = svmbs["f22_surf"]
    end
    
    # Interpolate Bond Pricing Surfaces
    if interp_surfs
        println("Interpolating Bond Pricing Surfaces...")
        svm = @time bpr_interp(svm)
    end
    
    return svm
end
