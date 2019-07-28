module_path = "/home/artur/BondPricing/Julia/modules/"
push!(LOAD_PATH, module_path)
# modnames = ["ModelObj", "AnalyticFunctions",
#             "EqFinDiff", "Batch"]
# for modl in modnames
#     if !(joinpath(module_path, modl) in LOAD_PATH)
#         push!(LOAD_PATH, joinpath(module_path, modl))
#     end
# end

module FullInfoEq 

using Distributed
using Dierckx

using Parameters
using Printf
using DataFrames
using CSV

using ModelObj: get_obj_model

using AnalyticFunctions: get_cvm_vb,
                         get_param

using EqFinDiff: eq_fd

using Batch: interp_values


# * Auxiliary Functions ##############################
# include("_fi_auxiliary_functions.jl")
function full_info_eq_deriv_root_search(svm, df; N::Int64=10^5,
                                        k::Int64=3, bc::String="extrapolate")
    
    # Interpolate Values
    eq_deriv_fun = Dierckx.Spline1D(df[:, :vb], df[:, :eq_deriv], k=3, bc="extrapolate")
    
    # Compute optimal VB:
    res = DataFrame(vb = [NaN])
    vbroots = roots(eq_deriv_fun; maxn=8)

    if !isempty(vbroots)
        eq_min_val_interp = Dierckx.Spline1D(df[:, :vb], df[:, :eq_min_val],
                                             k=3, bc="extrapolate")
        abs_eq_min_val = abs.(eq_min_val_interp(vbroots))
        res[!, :vb] .= vbroots[argmin(abs_eq_min_val)]
    else
        ref_vbgrid = range(minimum(df[:, :vb]), stop=maximum(df[:, :vb]), length=N)
        res[!, :vb] .= ref_vbgrid[argmin(abs.(eq_deriv_fun(ref_vbgrid)))]
    end
    
    # Equity Values
    res[!, :eq_deriv] .= eq_deriv_fun(res[:, :vb])

    # Interpolate Functions
    interp_cols = vcat([:debt, :equity],
                       [:eq_deriv_min_val, 
                        :eq_min_val, :eq_vb])
    res = interp_values(res, df, :vb, interp_cols; k=k, bc=bc)

    res[!, :eq_negative] .= (res[:, :eq_min_val] .< -.005)

    # Fixed Values
    return  res #non_interp_values(svm, res)
end


# * Optimal VB Functions #############################
# include("_fi_optimal_vb.jl")
function find_full_info_vb(svm, jks;
                           mu_b::Float64=NaN,
                           lb::Float64=.75, ub::Float64=1.25,
                           vbN::Int64=15, N::Int64=10^5)

    if isnan(mu_b)
        mu_b = jks.mu_b
    end
    
    if get_obj_model(svm) == "cvm"
        return get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    else
        vbl = get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=jks.m,
                          c=jks.c, p=jks.p)
    
        vbh = get_cvm_vb(svm, svm.pm.sigmal;
                         mu_b=mu_b, m=jks.m,
                         c=jks.c, p=jks.p)

        vbgrid = range(.75 * minimum([vbl, vbh]),
                       stop=1.25 * maximum([vbl, vbh]), length=vbN)

        res = @time fetch(@spawn [eq_fd(svm; vbl=vbl, mu_b=mu_b,
                                        m=jks.m, c=jks.c, p=jks.p)
                                  for vbl in vbgrid])

        eqdf = full_info_eq_deriv_root_search(svm, vcat(res...); N=N)

        return eqdf[1, :vb]
    end
end


function set_full_information_vb!(jf, jks;
                                  rerun_fi_vb::Bool=false,
                                  fi_sf_vb::Float64=NaN,
                                  fi_rf_vb::Float64=NaN,
                                  lb::Float64=.75,
                                  ub::Float64=1.25,
                                  vbN::Int64=20)

    if any([.&(isnan(jks.fi_sf_vb), isnan(fi_sf_vb)), rerun_fi_vb])
        fi_sf_vb = find_full_info_vb(jf.sf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_sf_vb, fi_sf_vb)
    elseif !isnan(fi_sf_vb)
         setfield!(jks, :fi_sf_vb, fi_sf_vb)       
    end

    if any([.&(isnan(jks.fi_rf_vb), isnan(fi_rf_vb)), rerun_fi_vb])
        fi_rf_vb = find_full_info_vb(jf.rf, jks; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_rf_vb, fi_rf_vb)
    elseif !isnan(fi_rf_vb)
         setfield!(jks, :fi_rf_vb, fi_rf_vb)       
    end

    return jks
end


# * Optimal Bond Measure #############################
# include("_fi_optimal_bond_measure.jl")
function bond_measure_full_info_vb(svm, jks, mu_b::Float64;
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5)

    jks2 = deepcopy(jks)
    setfield!(jks2, :mu_b, mu_b)

    return find_full_info_vb(svm, jks2;
                             lb=lb, ub=ub,
                             vbN=15, N=vbN2)
end


function form_mu_b_vb_pairs(svm, jks;
                            mu_b_min::Float64=.5,
                            mu_b_max::Float64=5.,
                            mu_bN::Int64=20,
                            mu_bN2::Int64=10^5,
                            lb::Float64=.75,
                            ub::Float64=1.25,
                            vbN::Int64=15,
                            vbN2::Int64=10^5,
                            spline_k::Int64=3,
                            spline_bc::String="extrapolate")
    tmp=DataFrame()
    exit = false
    count = 0
    while !exit
        println(string("COUNT: ", count))
        # Form mu_b grid
        mu_b_grid = range(mu_b_min, stop=mu_b_max; length=mu_bN)

        # Get Optimal VB for each mu_b value
        if get_obj_model(svm) == "cvm"
            vbl_list = [get_cvm_vb(svm, svm.pm.sigmal;
                               mu_b=mu_b, m=jks.m, c=jks.c, p=jks.p)
                        for mu_b in mu_b_grid]
        
        else
            vbl_list = fetch(@spawn [bond_measure_full_info_vb(svm, jks, mu_b;
                                                               lb=lb, ub=ub,
                                                               vbN=vbN, vbN2=vbN2)
                                     for mu_b in mu_b_grid])
        end

        # Run Equity Finite Differences Method ###############
        tmp = DataFrame(mu_b = mu_b_grid, vbl = vbl_list)

        # Evaluate exit condition
        loc1 = tmp[:vbl] .<= 1.05 * get_param(svm, :V0)
        cond1 = sum(loc1) <= .9 * mu_bN

        loc2 = tmp[:vbl] .> .9 * get_param(svm, :V0)
        cond2 = sum(loc2) < maximum([2, .1 * mu_bN])

        if cond1
            ref_mu_b_grid = range(mu_b_min, stop=mu_b_max; length=10^5)
            vblf = Dierckx.Spline1D(tmp[:mu_b], tmp[:vbl]; k=spline_k, bc=spline_bc)

            mu_b_max = ref_mu_b_grid[argmin(abs.(vblf(ref_mu_b_grid) .- 1.05 * get_param(svm, :V0)))]
        elseif cond2
            mu_b_max = 1.1 * mu_b_max
        else
            exit = true
        end
        count += 1
    end

    return tmp
end


function find_optimal_bond_measure(svm;
                                   mu_b_min::Float64=.5,
                                   mu_b_max::Float64=5.,
                                   mu_bN::Int64=20,
                                   mu_bN2::Int64=10^5,
                                   jks=JointKStruct(fill(NaN, 10)...),
                                   m::Float64=NaN,
                                   c::Float64=NaN,
                                   p::Float64=NaN,
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5,
                                   spline_k::Int64=3,
                                   spline_bc::String="extrapolate",
                                   tol::Float64=2.5 * 1e-3)

    # Set Capital Structure #########################
    if !isnan(m)
        setfield!(jks, :m, m)
    elseif isnan(jks.m)
        setfield!(jks, :m, jf.jks.m)
    end
    
    if !isnan(c)
        setfield!(jks, :c, c)
    elseif isnan(jks.c)
        setfield!(jks, :c, jf.jks.c)
    end

    if !isnan(p)
        setfield!(jks, :p, p)
    elseif isnan(jks.p)
        setfield!(jks, :p, jf.jks.p)
    end
    # ###############################################

    
    # Get Optimal VB for each mu_b value ############
    tmp = form_mu_b_vb_pairs(svm, jks;
                             mu_b_min=mu_b_min,
                             mu_b_max=mu_b_max,
                             mu_bN=mu_bN, mu_bN2=mu_bN2,
                             lb=lb, ub=ub,
                             vbN=vbN, vbN2=vbN2)
    # ###############################################

    
    # Run Equity Finite Differences Method ##########
    res = fetch(@spawn [eq_fd(svm; vbl=tmp[r, :vbl],
                              mu_b=tmp[r, :mu_b],
                              m=jks.m, c=jks.c, p=jks.p)
                        for r in 1:size(tmp, 1)])
    df = vcat(res...)
    # ###############################################

    
    # Interpolate VB and Firm Value in mu_b #########
    funs = Dict()
    for var in [:vb, :firm_value]
        funs[var] = Dierckx.Spline1D(df[:mu_b], df[var];
                                     k=spline_k, bc=spline_bc)
    end
    # ###############################################

    # Form mu_b refined grid
    mu_b_grid_ref = range(minimum(df[:mu_b]),
                          stop=maximum(df[:mu_b]), length=mu_bN2)

    # Get optimal mu_b
    opt_mu_b = mu_b_grid_ref[argmax(funs[:firm_value](mu_b_grid_ref))]

    # Set optimal mu_b
    setfield!(jks, :mu_b, opt_mu_b)

    # Get Optimal VB
    opt_vbl = bond_measure_full_info_vb(svm, jks, jks.mu_b)
    setfield!(jks, :vbl, opt_vbl)
    
    eqdf = eq_fd(svm; vbl=jks.vbl,
                 mu_b=jks.mu_b,
                 m=jks.m, c=jks.c, p=jks.p)

    
    # Filter to eliminate mu_b candidates for which there is not optimal vbl
    count = 1
    while .&(abs.(eqdf[1, :eq_deriv_min_val]) > tol, count<=10)
        funs[:eq_deriv_min_val] = Dierckx.Spline1D(df[:mu_b], abs.(df[:eq_deriv_min_val]);
                                                   k=spline_k, bc=spline_bc)
        # Filter mu_b grid
        filtered_mu_b_grid = mu_b_grid_ref[funs[:eq_deriv_min_val](mu_b_grid_ref) .< tol]
        
        # Get optimal mu_b
        opt_mu_b = filtered_mu_b_grid[argmax(funs[:firm_value](filtered_mu_b_grid))]
        
        # Set optimal mu_b
        setfield!(jks, :mu_b, opt_mu_b)
        
        # Get Optimal VB
        opt_vbl = bond_measure_full_info_vb(svm, jks, jks.mu_b)
        setfield!(jks, :vbl, opt_vbl)
        
        eqdf = eq_fd(svm; vbl=jks.vbl,
                     mu_b=jks.mu_b,
                     m=jks.m, c=jks.c, p=jks.p)

        df = sort!(vcat(df, eqdf), :mu_b)
        count += 1
    end

    return eqdf
end

# * END MODULE
end
