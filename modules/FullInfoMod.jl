module_path = "/home/artur/BondPricing/Julia/modules/"
push!(LOAD_PATH, module_path)

module FullInfoMod

using Printf
using DataFrames
using CSV
using Dates

using Batch: dfcols

using FullInfoEq: find_optimal_bond_measure,
                  JointKStruct

using JointEq: JointFirms, store_ep_params




# * Optimal Bond Measure



# * Full Information Equilibrium Functions
function recompute_fi(jf, df::DataFrame, ft::Symbol, rmp::Symbol, ep_jks)
    cond = .&(df[:, :type] .== ft, df[:, :rmp] .== rmp)
    fv(df) = .&(!isempty(df), :firm_value in names(df)) ? df[1, :firm_value] : NaN
    tmp = df[cond, :]
    
    fr = getfield(getfield(jf, ft), rmp).fr
    if .&(isnan(fv(tmp)), !isnothing(fr))
        tmp = find_optimal_bond_measure(fr, jf.bc; jks=ep_jks)
        tmp[!, :eq_type] .= :full_info
        tmp[!, :datetime] .= Dates.now()
        tmp[!, :type] .= ft
        tmp[!, :rmp] .= rmp
        df = vcat([df[cond .== false, :], tmp]...)
    end
        
    return df
end


# Electronic Platform Misrepresentation ###########################
    # Do risky firms have an incentive to copy the capital structure
    # of the safe firms?
    if run_misrep
        misrep_jks = deepcopy(ep_jks)
        setfield!(misrep_jks, :mu_s, 1.)
        if !isnan(fi_sf_mu_b)
            setfield!(misrep_jks, :mu_b, fi_sf_mu_b)
        end
        
        ep_misrep_eqdf = find_joint_optimal_vb(ep_jf, misrep_jks;
                                               mu_b=fi_sf_mu_b,
                                               rerun_fi_vb=true)

        # Add Objective Function columns
        ep_misrep_eqdf[!, :obj_fun] .= String(sf_obj_fun)
        ep_misrep_eqdf[2, :obj_fun] = "misrep"
        ep_misrep_eqdf[!, :eq_type] .= "misrep"
        
        # Reshape
        ep_misrep_eqdf = reshape_sf_rf_df(ep_misrep_eqdf)
    else
        ep_misrep_eqdf = DataFrame()
    end
    # ##################################################################


    

end
