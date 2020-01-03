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

using Dates

using ModelObj: Firm, get_obj_model

using AnalyticFunctions: get_cvm_vb,
                         get_param

using EqFinDiff: eq_fd

using Batch: interp_values, dfcols

# using JointEq: JointKStruct
using JointEqStructs: BondContract, JointFirms,
                      JointKStruct, type_fun

# * Containers
types_order = vcat([Symbol(x, :_, y) for x in [:st, :rt], 
                    y in [:iota, :lambda, :sigmah]]...)

types_names = [:iota, :lambda, :sigmah]
cp_names = [:V0, :alpha, :gross_delta, :pi, :r, :xi, :sigmal]
fi_df_names =  vcat([x for x in dfcols if !occursin("diff", String(x))],
                    [:eq_type, :datetime, :type, :rmp])

typeval = Dict{DataType, Any}(Float64 => NaN,
                              Bool => false,
                              String => "",
                              Symbol => :full_info,
                              DateTime => DateTime(0))


# * Auxiliary Functions
function get_types_comb_df(types_dict::Dict{Symbol, Array{Float64, N} where N})
    types_combs =  [Array([x1, x2, x3, x4, x5, x6]) for x1 in types_dict[types_order[1]], 
                                                        x2 in types_dict[types_order[2]], 
                                                        x3 in types_dict[types_order[3]],
                                                        x4 in types_dict[types_order[4]], 
                                                        x5 in types_dict[types_order[5]], 
                                                        x6 in types_dict[types_order[6]]]

    df = DataFrame(hcat(types_combs...)')
    rename!(df, types_order)
    df[!, :pair_num] .= 1:size(df, 1)
    
    return df
end


function get_empty_fidf(; col_names::Array{Symbol,1}=fi_df_names)
    tmp = [eval(type_fun(x)[]) for x in fi_df_names]
    return DataFrame(tmp, fi_df_names)
end


function check_common_params(jf, df::DataFrame)
    for x in [:V0, :alpha, :gross_delta, :r, :xi]
        z = unique(df[:,x])
        if !.&(length(z) == 1, abs.(z[1] .- getfield(jf.cp, x)) .< 1e-5)
            println("Error! Multiple values for the same param in common params list.")
            println("Exiting...")
            return false
        end
    end
    println("All common parameters match!")
    return true
end


function set_missing_values(jf, ft::Symbol)
    reminder = [x for x in fi_df_names if 
                !(x in vcat(cp_names, types_names))]
    tmp = Dict{Symbol, Any}([x => getfield(jf.cp, x) for x in cp_names])
    for x in types_names
        tmp[x] = getfield(jf.td, Symbol(ft,:_, x))
    end
    for x in reminder
        tmp[x] = typeval[type_fun(x)] 
    end

    return DataFrame(tmp)
end


# Check that type can choose not to manage the volatility risk   
function check_nrm(jf, ft::Symbol)
    return any([[!isnan(getfield(jf.td, Symbol(ft, :_, x))) for x in [:lambda, :sigmah]]...,
                abs(getfield(jf.td, Symbol(ft, :_lambda))) .< 1e-5,
                abs(getfield(jf.td, Symbol(ft, :_sigmah))- jf.cp.sigmal) .< 1e-5])
end


function extract_fi_results(jf, df::DataFrame)
    if !check_common_params(jf, df)
        return
    else
        dflist = []
        for ft in [:st, :rt]
            # Risk-Management
            sigmal_cond = abs.(df[:,:sigmal] .- jf.cp.sigmal) .< 1e-5
            iota = getfield(jf.td, Symbol(ft,:_iota))
            iota_cond = isnan(iota) ? isnan.(df[:, :iota]) : abs.(df[:, :iota] .- iota) .< 1e-5
            cond = .&(sigmal_cond, iota_cond)
            if !isempty(df[cond, :])
                rmdf = df[cond, :]
            else
                rmdf = set_missing_values(jf, ft)
            end
            rmdf[!, :type] .= ft
            rmdf[!, :rmp] .= :rm
            push!(dflist, rmdf)
            

            # No-Risk-Management
            if check_nrm(jf, ft)
                cond =  abs.(df[:,:sigmal] .- jf.cp.sigmal) .< 1e-5
                for x in [:lambda, :sigmah]
                    xval = getfield(jf.td, Symbol(ft,:_, x))
                    tmp_cond = isnan(xval) ? isnan.(df[:,x]) : abs.(df[:,x] .- xval) .< 1e-5
                    cond = .&(cond, tmp_cond)
                end
                if !isempty(df[cond, :])
                    nrmdf = df[cond, :]
                else
                    nrmdf = set_missing_values(jf, ft)
                end
                nrmdf[!, :type] .= ft
                nrmdf[!, :rmp] .= :nrm
                push!(dflist, nrmdf)
            end
        end
        
        return vcat(dflist...)
    end
end


# * Root Search Functions
# include("_fi_auxiliary_functions.jl")
function full_info_eq_deriv_root_search(svm::Firm, df::DataFrame; N::Int64=10^5,
                                        k::Int64=3, bc::String="extrapolate")
    
    # Interpolate Values
    eq_deriv_fun = Dierckx.Spline1D(df[:, :vb], df[:, :eq_deriv], k=3, bc="extrapolate")
    
    # Compute optimal VB:
    res = DataFrame(:vb => NaN)
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


# * Optimal VB Functions
# include("_fi_optimal_vb.jl")
function find_full_info_vb(svm::Firm, bc, mu_b::Float64;
                           lb::Float64=.75, ub::Float64=1.25,
                           vbN::Int64=15, N::Int64=10^5)

    if get_obj_model(svm) == "cvm"
        return get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=bc.m,
                          c=bc.c, p=bc.p)
    else
        vbl = get_cvm_vb(svm, svm.pm.sigmal;
                          mu_b=mu_b, m=bc.m,
                          c=bc.c, p=bc.p)
    
        vbh = get_cvm_vb(svm, svm.pm.sigmal;
                         mu_b=mu_b, m=bc.m,
                         c=bc.c, p=bc.p)

        vbgrid = range(.75 * minimum([vbl, vbh]),
                       stop=1.25 * maximum([vbl, vbh]), length=vbN)

        res = @time fetch(@spawn [eq_fd(svm; vbl=vbl, mu_b=mu_b,
                                        m=bc.m, c=bc.c, p=bc.p)
                                  for vbl in vbgrid])

        eqdf = full_info_eq_deriv_root_search(svm, vcat(res...); N=N)

        return eqdf[1, :vb]
    end
end


# * Optimal Bond Measure
# include("_fi_optimal_bond_measure.jl")
function bond_measure_full_info_vb(svm::Firm, bc,
                                   mu_b::Float64;
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5)

    return find_full_info_vb(svm, bc, mu_b;
                             lb=lb, ub=ub,
                             vbN=15, N=vbN2)
end


function form_mu_b_vb_pairs(svm::Firm, bc;
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

    tmp = DataFrame()
    exit = false
    count = 0
    while !exit
        println(string("COUNT: ", count))

        # Form mu_b grid
        mu_b_grid = range(mu_b_min, stop=mu_b_max; length=mu_bN)

        # Get Optimal VB for each mu_b value
        if get_obj_model(svm) == "cvm"
            vbl_list = [get_cvm_vb(svm, svm.pm.sigmal;
                                   mu_b=mu_b, m=bc.m, c=bc.c, p=bc.p)
                        for mu_b in mu_b_grid]
        
        else
            vbl_list = fetch(@spawn [bond_measure_full_info_vb(svm, bc, mu_b;
                                                               lb=lb, ub=ub,
                                                               vbN=vbN, vbN2=vbN2)
                                     for mu_b in mu_b_grid])
        end

        # Run Equity Finite Differences Method ###############
        tmp = DataFrame(:mu_b => mu_b_grid, :vbl => vbl_list)

        # Evaluate exit condition
        loc1 = tmp[:, :vbl] .<= 1.05 * get_param(svm, :V0)
        cond1 = sum(loc1) <= .9 * mu_bN

        loc2 = tmp[:, :vbl] .> .9 * get_param(svm, :V0)
        cond2 = sum(loc2) < maximum([2, .1 * mu_bN])

        if cond1
            ref_mu_b_grid = range(mu_b_min, stop=mu_b_max; length=10^5)
            vblf = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, :vbl]; k=spline_k, bc=spline_bc)

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


function find_optimal_bond_measure(svm::Firm, bc;
                                   mu_b_min::Float64=.5,
                                   mu_b_max::Float64=5.,
                                   mu_bN::Int64=20,
                                   mu_bN2::Int64=10^5,
                                   lb::Float64=.75,
                                   ub::Float64=1.25,
                                   vbN::Int64=15,
                                   vbN2::Int64=10^5,
                                   spline_k::Int64=3,
                                   spline_bc::String="extrapolate",
                                   tol::Float64=2.5 * 1e-3)

    # Get Optimal VB for each mu_b value ############
    tmp = form_mu_b_vb_pairs(svm, bc;
                             mu_b_min=mu_b_min,
                             mu_b_max=mu_b_max,
                             mu_bN=mu_bN, mu_bN2=mu_bN2,
                             lb=lb, ub=ub,
                             vbN=vbN, vbN2=vbN2)
    # ###############################################

    
    # Run Equity Finite Differences Method ##########
    res = fetch(@spawn [eq_fd(svm; vbl=tmp[r, :vbl],
                              mu_b=tmp[r, :mu_b],
                              m=bc.m, c=bc.c, p=bc.p)
                        for r in 1:size(tmp, 1)])
    df = vcat(res...)
    # ###############################################

    
    # Interpolate VB and Firm Value in mu_b #########
    funs = Dict()
    for var in [:vb, :firm_value]
        funs[var] = Dierckx.Spline1D(df[:, :mu_b], df[:, var];
                                     k=spline_k, bc=spline_bc)
    end
    # ###############################################

    # Form mu_b refined grid
    mu_b_grid_ref = range(minimum(df[:, :mu_b]),
                          stop=maximum(df[:, :mu_b]), length=mu_bN2)

    # Get optimal mu_b
    opt_mu_b = mu_b_grid_ref[argmax(funs[:firm_value](mu_b_grid_ref))]

    # Get Optimal VB
    opt_vbl = bond_measure_full_info_vb(svm, bc, opt_mu_b)

    
    eqdf = eq_fd(svm; vbl=opt_vbl,
                 mu_b=opt_mu_b,
                 m=bc.m, c=bc.c, p=bc.p)

    
    # Filter to eliminate mu_b candidates for which there is not optimal vbl
    count = 1
    while .&(abs.(eqdf[1, :eq_deriv_min_val]) > tol, count<=10)
        funs[:eq_deriv_min_val] = Dierckx.Spline1D(df[:, :mu_b], abs.(df[:, :eq_deriv_min_val]);
                                                   k=spline_k, bc=spline_bc)
        # Filter mu_b grid
        filtered_mu_b_grid = mu_b_grid_ref[funs[:eq_deriv_min_val](mu_b_grid_ref) .< tol]
        
        # Get optimal mu_b
        opt_mu_b = filtered_mu_b_grid[argmax(funs[:firm_value](filtered_mu_b_grid))]
        
        # Get Optimal VB
        opt_vbl = bond_measure_full_info_vb(svm, bc, opt_mu_b)
        
        eqdf = eq_fd(svm; vbl=opt_vbl,
                     mu_b=opt_mu_b,
                     m=bc.m, c=bc.c, p=bc.p)

        df = sort!(vcat(df, eqdf), :mu_b)
        count += 1
    end

    return eqdf
end


# * Full Information Equilibrium Functions
function compute_fi_eq(jf, df::DataFrame, ft::Symbol, rmp::Symbol)
    cond = .&(df[:, :type] .== ft, df[:, :rmp] .== rmp)
    fv(df) = .&(!isempty(df), :firm_value in names(df)) ? df[1, :firm_value] : NaN
    tmp = df[cond, :]
    
    fr = getfield(getfield(jf, ft), rmp).fr
    if .&(isnan(fv(tmp)), !isnothing(fr))
        tmp = find_optimal_bond_measure(fr, jf.bc)
        tmp[!, :eq_type] .= :full_info
        tmp[!, :datetime] .= Dates.now()
        tmp[!, :type] .= ft
        tmp[!, :rmp] .= rmp
        df = vcat([df[cond .== false, :], tmp]...)
    end
        
    return df
end


function get_fi_results(jf, fi_fpath_name::String;
                        rerun_full_info::Bool=false,
                        save_results::Bool=true)
    
    if .&(isfile(fi_fpath_name), !rerun_full_info)
        fidf = extract_fi_results(jf, CSV.read(fi_fpath_name))
    else
        fidf=get_empty_fidf()
        rerun_full_info = true
    end
    
    if rerun_full_info
        for ft in [:st, :rt], rmp in [:rm, :nrm]
           fidf = compute_fi_eq(jf, fidf, ft, rmp)
        end
    end

    if save_results
        CSV.write(fi_fpath_name, fidf)
    end
    
    return fidf
end


function get_fi_eq(jf; fidf::DataFrame=DataFrame())
    fidf = isempty(fidf) ? get_empty_fidf() : fidf
    fieqdf = get_empty_fidf()
    for ft in [:st, :rt]
        tmp = fidf[fidf[:, :type] .== ft, :]
        if size(tmp, 1) < 1
            for rmp in [:rm, :nrm]
                tmp = compute_fi_eq(jf, tmp, ft, rmp)
            end
            fidf = vcat([fidf, DataFrame(tmp)]...)
        end
        fieqdf = vcat([fieqdf, DataFrame(tmp[argmax(tmp[:, :firm_value]), :])]...)
    end
    
    return fidf, fieqdf
end


# * Joint Capital Structure
function set_full_information_vb!(sf, rf, jks;
                                  rerun_fi_vb::Bool=false,
                                  fi_st_vb::Float64=NaN,
                                  fi_rt_vb::Float64=NaN,
                                  lb::Float64=.75,
                                  ub::Float64=1.25,
                                  vbN::Int64=20)

    # Form Bond Contract
    bc = BondContract(jks.m, jks.c, jks.p)

    if any([.&(isnan(jks.fi_st_vb), isnan(fi_st_vb)), rerun_fi_vb])
        fi_st_vb = find_full_info_vb(sf, bc, jks.mu_b; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_st_vb, fi_st_vb)
    elseif !isnan(fi_st_vb)
         setfield!(jks, :fi_st_vb, fi_st_vb)       
    end

    if any([.&(isnan(jks.fi_rt_vb), isnan(fi_rt_vb)), rerun_fi_vb])
        fi_rt_vb = find_full_info_vb(rf, bc, jks.mu_b; lb=lb, ub=ub, vbN=vbN)
        setfield!(jks, :fi_rt_vb, fi_rt_vb)
    elseif !isnan(fi_rt_vb)
         setfield!(jks, :fi_rt_vb, fi_rt_vb)       
    end

    return jks
end


# * END MODULE
end
