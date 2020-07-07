using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
using PyPlot
# using PyCall
using Seaborn
using LaTeXStrings

module_path = string(main_path, "/", "Julia/modules/")
modls = ["Batch", "ModelObj",
         "AnalyticFunctions", "BondPrInterp",
         "EqFinDiff", "JointEqStructs", "FullInfoEq",
         "ModelPlots", "JointEq"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# * Parameters
firm_obj_fun = :firm_value
s_iota = 0.0
k_ep = 25 * 1e-4
k_otc = 32.5

min_lambda = .0
max_lambda = .5

# ** Joint Equilibrium
pool_jeq_objf = :st_firm_value
sep_jeq_objf = :st_firm_value

# ** DataFrame files Inputs
del_old_files = true
save_filtered_df = true

# * Unconstrained (OTC) Results
# ###############################################################################
# Get OTC Results ###############################################################
# ###############################################################################
# Transaction Costs and Volatility Risk Parameters
cvmdict = Dict{Symbol,Array{Float64,1}}(:sigmal => [0.15],
                                        :m  => [1.],
                                        :gross_delta => [0.02],
                                        :mu_b => [1.0],
                                        :xi => [1.0],
                                        :iota => [x for x in Batch.cvm_param_values_dict[:iota]
                                                  if .&(x >= s_iota, x <= 20. * 1e-4)])
svmdict = deepcopy(cvmdict)
svmdict[:kappa] = [k_ep]
svmdict[:iota] = [.0]
svmdict[:sigmah] = Batch.svm_param_values_dict[:sigmah]
svmdict[:lambda] = [Batch.svm_param_values_dict[:lambda][4]]

# Get Safe and Risky Firms' Full Info Optimal Results ##########################
cvmdf, svmdf, _ = ModelPlots.get_cvm_svm_dfs(cvmdict, svmdict;
                                             firm_obj_fun=firm_obj_fun)


# Safe Type's Firm Value as a function of Transaction Costs Kappa
cond = abs.(cvmdf[:, :iota] .- s_iota) .< 1e-5
fi_fv_fun = Dierckx.Spline1D(cvmdf[cond, :kappa] .* 1e4, cvmdf[cond, :firm_value];
                             k=3, bc="extrapolate")

# * Form Joint Firm Object and Set Path
k_otc = k_ep
mu_s = .2
_, jfep = JointEq.otc_ep_jfs(k_otc, k_ep;
                             mu_s=mu_s,
                             st_iota=s_iota,
                             st_lambda=NaN,
                             st_sigmah=NaN,
                             rt_iota=s_iota,
                             rt_lambda=svmdict[:lambda][1],
                             rt_sigmah=svmdict[:sigmah][1])

# jks_fpath = string(main_path, "/Julia/Results/JEQ/kappa_25.00_bp__sigmal_0.150/m_1.00_pcr_12.00")
jks_fpath = JointEq.make_jeq_jks_fpath(jfep)

# * Load Results
fidf = JointEqStructs.get_unique_df(jks_fpath, :fi;
                                    save_filtered_df=save_filtered_df,
                                    del_old_files=del_old_files)
fidf = ModelPlots.nrmp_remove_extra_vals(fidf)

misrepdf = JointEqStructs.get_unique_df(jks_fpath, :misrep;
                                        save_filtered_df=save_filtered_df,
                                        del_old_files=del_old_files)
misrepdf = ModelPlots.nrmp_remove_extra_vals(misrepdf)

pooldf = JointEqStructs.get_unique_df(jks_fpath, :pool;
                                      save_filtered_df=save_filtered_df,
                                      del_old_files=del_old_files)
pooldf = ModelPlots.nrmp_remove_extra_vals(pooldf)

sepdf = JointEqStructs.get_unique_df(jks_fpath, :sep;
                                     save_filtered_df=save_filtered_df,
                                     del_old_files=del_old_files)
sepdf = ModelPlots.nrmp_remove_extra_vals(sepdf)


# ** Convert Risk-Management and Transaction Costs to Basis Points
# %% Adjust kappa and iota to raw values
for col in [:kappa, :iota]
    cond = fidf[:, col] .> 1
    if any(cond)
        fidf[cond, col] .= fidf[cond, col] .* 1e-4
    end
end
for df in [misrepdf, pooldf, sepdf]
    for col in [:kappa, :r_iota, :s_iota]
        cond = df[:, col] .> 1
        if any(cond)
            df[cond, col] .= df[cond, col] .* 1e-4
        end
    end
end

# Iota and kappa in basis points
fidf[!, :iota] .= fidf[:, :iota] .* 1e4
fidf[!, :kappa] .= fidf[:, :kappa] .* 1e4
for df in [misrepdf, sepdf, pooldf]
    for ft in [:s_, :r_]
        df[!, Symbol(ft, :iota)] .= df[:, Symbol(ft, :iota)] .* 1e4
    end

    df[!, :kappa] .= df[:, :kappa] .* 1e4
end

# * Filter FI & Misrep Results
fidf[isnan.(fidf[:, :nrm_lambda]), :sigmah] .= unique(fidf[:, :sigmal])[1]
fidf[isnan.(fidf[:, :nrm_lambda]), :nrm_sigmah] .= unique(fidf[:, :sigmal])[1]
fidf[isnan.(fidf[:, :nrm_lambda]), :lambda] .= svmdict[:lambda][1]
fidf[isnan.(fidf[:, :nrm_lambda]), :nrm_lambda] .= svmdict[:lambda][1]
fidf = fidf[fidf[:, :nrm_lambda] .== svmdict[:lambda][1], :]

misrepdf = misrepdf[misrepdf[:, :r_nrm_lambda] .== svmdict[:lambda][1], :]

# * Filter JEQ Results
# ** Pooling Payoffs
pool_eq_cond = .&(pooldf[:, :r_nrm_lambda] .== svmdict[:lambda][1],
                  pooldf[:, :jeq_objf] .== pool_jeq_objf)
pooldf = pooldf[pool_eq_cond, :]

# ** Misrep and Separating Payoffs
# if misrep_type == :st, then mu_s = 0
# if misrep_type == :rt, then mu_s = 1
sep_eq_cond = .&(sepdf[:, :r_nrm_lambda] .== svmdict[:lambda][1],
                  sepdf[:, :jeq_objf] .== sep_jeq_objf)
sepdf = sepdf[sep_eq_cond, :]
for df in [misrepdf, sepdf]
    rt_misrep_cond = .&(abs.(df[:, :mu_s] .- 1.) .< 1e-5,
                        df[:, :misrep_type] .== :rt)
    st_misrep_cond = .&(abs.(df[:, :mu_s]) .< 1e-5,
                        df[:, :misrep_type] .== :st)
    df = df[(st_misrep_cond .| rt_misrep_cond), :]
end
sepdf = sepdf[sepdf[:, :jeq_objf] .== sep_jeq_objf, :]


# * Equilibrium Type Indicators and Equilibrium-Specific Payoffs
# ** EP v.s. OTC Firm Value
fi_fv=fidf[1, :firm_value]
println(string("s_fi_fv:  ", fi_fv))
k_ep = 25.
k_otc=32.5


# Safe Type's Firm Value in the EP Full Information Equilibrium:
s_fi_fv = fi_fv_fun(k_ep) # fidf[abs.(fidf[:iota] .- s_iota * 1e4) .< 1e-5, :firm_value][1]
println(string("s_fi_fv 2: ", s_fi_fv))

# Safe Type's Firm Value in the EP Full Information Equilibrium:
s_otc_fv = fi_fv_fun(k_otc)
println(string("s_otc_fv: ", s_otc_fv))

max_mu_s = .75
pooldf = pooldf[pooldf[:, :mu_s] .<= max_mu_s, :]
mu_s_vec = unique(pooldf[:, :mu_s])

# ** Full Information Plots
min_sigmah = minimum(fidf[:, :sigmal])
fid = ModelPlots.nrmp_contour_plot_prepare_df(fidf;
                                              min_lambda=min_lambda,
                                              max_lambda=max_lambda,
                                              min_sigmah=min_sigmah,
                                              xvar=plot_xvar,
                                              yvar=plot_yvar,
                                              yvals=vcat(.0, mu_s_vec))

fid = ModelPlots.nrmp_surface_interpolator(fid)


# ** Misrepresentation
mpd = ModelPlots.nrmp_contour_plot_prepare_df(misrepdf;
                                              misrep_type=:rt,
                                              min_lambda=min_lambda,
                                              max_lambda=max_lambda,
                                              min_sigmah=min_sigmah,
                                              xvar=plot_xvar,
                                              yvar=plot_yvar,
                                              yvals=mu_s_vec)

mpd = ModelPlots.nrmp_contour_plot_additional_zvars(fid, mpd, :misrep)


# ** Pooling Payoffs
poold = ModelPlots.nrmp_contour_plot_prepare_df(pooldf;
                                                misrep_type=:rt,
                                                min_lambda=min_lambda,
                                                max_lambda=max_lambda,
                                                min_sigmah=min_sigmah,
                                                xvar=plot_xvar,
                                                yvar=plot_yvar)

poold = ModelPlots.nrmp_contour_plot_additional_zvars(fid, poold, :pool)


# ** Separating Payoffs
sepd = ModelPlots.nrmp_contour_plot_prepare_df(sepdf;
                                               misrep_type=:rt,
                                               min_lambda=min_lambda,
                                               max_lambda=max_lambda,
                                               min_sigmah=min_sigmah,
                                               xvar=plot_xvar,
                                               yvar=plot_yvar,
                                               yvals=mu_s_vec)

sepd = ModelPlots.nrmp_contour_plot_additional_zvars(fid, sepd, :sep)

# ** JEQ Payoffs
fun_dict = ModelPlots.get_nrmp_contour_equilibria_funs(fid[:interp], mpd[:interp],
                                                       poold[:interp], sepd[:interp],
                                                       k_ep, k_otc, fi_fv_fun)

pltd = ModelPlots.get_eq_contour_mesh_grid(poold[:xvals], poold[:yvals], fun_dict, N=10^3)


# * Old
# pooldf = JointEq.get_jeq_results(jfep, jks,
#                                 :pool, :exp_firm_value,
#                                 jeqid; load_df=true,
#                                 recompute_df=false,
#                                 save_df=true)
# sepdf = JointEq.get_jeq_results(jfep, jks,
#                                 :sep, :st_firm_value,
#                                 jeqid; load_df=true,
#                                 recompute_df=false,
#                                 save_df=true)

# Old
# fidf = JointEq.process_results(jks_fpath, "full_info")
# misrepdf = JointEq.process_results(jks_fpath, "misrep")
# pooldf = JointEq.process_results(jks_fpath, "pooling")
# sepdf = JointEq.process_results(jks_fpath, "separating")
