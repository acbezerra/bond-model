function grid_creator(z0, z1, n)
    dz = (z1 - z0) / float(n)
    if z0^2 < 1e-6
        return dz, range(z0+dz, stop=z1, length=n)
    else
        return dz, range(z0, stop=z1, length=n)
    end
    # return dz, range(z0, stop=z1, length=n) 
end


function set_bpr_inputs_struct(bprf)
    return BPrInputs(bprf[:vtmax], 
                     bprf[:vtN],
                     bprf[:ttm_max],
                     bprf[:ttmN],
                     bprf[:vbhlmin],
                     bprf[:vbhlmax],
                     bprf[:vbhlN],
                     bprf[:vmax],
                     bprf[:vN], 
                     bprf[:uN], 
                     bprf[:vtN_ref],
                     bprf[:ttmN_ref],
                     bprf[:vbhlN_ref])
end


function set_bpr_surfs_struct(bprf)
    return BPrSurfs(bprf[:vtgrid],
                    bprf[:ttmgrid],
                    bprf[:vbhlgrid],
                    bprf[:f11_surf],
                    bprf[:f12_surf], 
                    bprf[:f13_surf], 
                    bprf[:f21_surf], 
                    bprf[:f22_surf])
end


function set_bpr_grids(svm)
    _, svm.bs.vtgrid = grid_creator(0.0, svm.bi.vtmax, svm.bi.vtN)
    _, svm.bs.ttmgrid = grid_creator(0.0, svm.bi.ttm_max, svm.bi.ttmN)
    _, svm.bs.vbhlgrid = grid_creator(svm.bi.vbhlmin, svm.bi.vbhlmax, svm.bi.vbhlN)

    return svm
end


function set_bpr_surfs(svm, bprf)
    svm.bs.f11_surf = bprf[:f11_surf]
    svm.bs.f12_surf = bprf[:f12_surf]
    svm.bs.f13_surf = bprf[:f13_surf]
    svm.bs.f21_surf = bprf[:f21_surf]
    svm.bs.f22_surf = bprf[:f22_surf]

    return svm
end


function set_bpr_funs_struct(svm, bprf)
    svm.bf.f11 = bprf[:f11]
    svm.bf.f12 = bprf[:f12]
    svm.bf.f13 = bprf[:f13]
    svm.bf.f21 = bprf[:f21]
    svm.bf.f22 = bprf[:f22]

    return svm
end


function set_opt_k_struct(mobj, df::DataFrame; m::Float64=NaN, tol::Float64=1e-6)

    if mobj.model == "svm"
        vars = fieldnames(ModelObj.FirmParams)
    elseif mobj.model == "cvm"
        vars = [x for x in fieldnames(ModelObj.FirmParams) if
                !(x in [:lambda, :sigmah])]
    else
        println("Please enter either a SVM or a CVM object. Exiting...")
        return
    end

    # Find Index
    LL = []
    for var in vars 
        append!(LL, [abs.(df[var] .- getfield(mobj.pm, var)) .< tol])
    end

    if !isnan(m)
        append!(LL, [abs.(df[:m] .- m) .< tol])
    end
    index = .&(LL...)


    # Set Capital Structure
    for var in [:mu_b, :m, :c, :p]
        setfield!(mobj.optKS, var, df[index, var][1])
    end
    mobj.optKS.vbl = df[index, :vb][1]
    
    return mobj
end

