
function firm_constructor(dict::Dict{Symbol,Float64};
                          bpr_inputs::Dict{Symbol, Real}=bpr_inputs_dict,
                          model::String="svm")
    for par in [:vbl, :vbh, :c, :p]
        if !(par in keys(dict))
            println(string("Setting initial ", par, " value to ", NaN))
            dict[par] = NaN
        end
    end

    if !(:mu_b in keys(dict)) | isnan(dict[:mu_b])
        mu_b = 1.0
        println(string("Setting measure of bonds to ", mu_b))
        dict[:mu_b] = mu_b
    end

    if model == "cvm"
        println("Constant Volatility Model: setting vbh to vbl, lambda to NaN")
        dict[:lambda] = NaN
        dict[:vbh] = dict[:vbl]
    end

    # Set Firm Parameters (dynamics, default, liquidity)
    firm_params =  FirmParams(dict[:V0], dict[:alpha], dict[:pi],
	                      dict[:r], dict[:gross_delta], dict[:iota],
	                      dict[:xi], dict[:kappa],
	                      dict[:lambda], dict[:sigmal], dict[:sigmah])

    # Bond Pricing Inputs
    ttm_max = dict[:m]
    bpr_surfs_inputs = BPrInputs(bpr_inputs[:vtmax], 
                                 bpr_inputs[:vtN],
                                 ttm_max,
                                 bpr_inputs[:ttmN],
                                 bpr_inputs[:vbhlmin],
                                 bpr_inputs[:vbhlmax],
                                 bpr_inputs[:vbhlN],
                                 bpr_inputs[:vmax],
                                 bpr_inputs[:vN], 
                                 bpr_inputs[:uN], 
                                 bpr_inputs[:vtN_ref],
                                 bpr_inputs[:ttmN_ref],
                                 bpr_inputs[:vbhlN_ref])


    # Bond Pricing Grids & Surfaces
    _, vtgrid = grid_creator(0.0, bpr_surfs_inputs.vtmax, bpr_surfs_inputs.vtN)
    _, ttmgrid = grid_creator(0.0, bpr_surfs_inputs.ttm_max, bpr_surfs_inputs.ttmN)
    _, vbhlgrid = grid_creator(bpr_surfs_inputs.vbhlmin, bpr_surfs_inputs.vbhlmax, bpr_surfs_inputs.vbhlN)
        
    mat1 = fill(NaN, 1, 1)
    mat2 = fill(NaN, 1, 1, 1)
        
    bpr_surfs = BPrSurfs(vtgrid, ttmgrid, vbhlgrid,
                     mat1, mat2, mat2, mat1, mat1)

    # Bond Pricing Interpolated Functions 
    bpr_interp_funs = BPrInterpFuns(nothing, nothing, nothing,
                                    nothing, nothing)

    bpr_fixed_ttm_inputs = BPrFixedTTMInputs(NaN, 0, NaN, NaN, NaN, 0)

    bpr_fixed_ttm_interp_funs = BPrInterpFuns(nothing, nothing, nothing,
                                nothing, nothing)

    # Optimal Capital Structs
    optKS = KStruct(NaN, NaN, NaN, NaN, NaN) 
  

    # Construct Firm Object
    return Firm(dict[:mu_b], dict[:m],
                dict[:c], dict[:p],
	    	dict[:vbl], dict[:vbh],
	        firm_params, model,
                bpr_surfs_inputs,
                bpr_surfs, 
                bpr_interp_funs,
                bpr_fixed_ttm_inputs,
                bpr_fixed_ttm_interp_funs,
                optKS)

end


function create_misrep_type(svmi, svmj)
    # Extract Type-i parameters
    ijdict = Dict{String, Float64}()
    for par in vcat(obj_params, firm_params)
        ijdict[string(par)] = extract_param(svmi, string(par))
    end

    # Set iota to zero
    ijdict[:iota] = .0

    # Copy capital structure from Type-j
    for par in [:mu_b, :c, :p]
        ijdict[par] = extract_param(svmj, par)
    end

    return firm_constructor(ijdict)
end
