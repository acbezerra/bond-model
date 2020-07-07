# %% codecell
using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
# using PyPlot
using Seaborn

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
# include(string(module_path, "/", "TestFunctions.jl"))
modls = ["Batch", "ModelObj", "AnalyticFunctions",
         "BondPrInterp", "EqFinDiff",
         "JointEqStructs", "FullInfoEq", "JointEq", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

ENV["LINES"] = 750
ENV["COLUMNS"] = 10000

# %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# %% codecell
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
plot_script_name = "nrmp_fixed_mu_s_contour_plots.jl"
include(string(plot_script_path, "/", plot_script_name))




# %% codecell
xvar = :sigmah
yvar = :mu_s
plot_xvar = :nrm_sigmah
plot_yvar = :mu_s
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
plot_script_name = "nrmp_fixed_lambda_contour_plots_init.jl"
include(string(plot_script_path, "/", plot_script_name))


# %%
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
plot_script_name = "nrmp_fixed_lambda_contour_plots.jl"
include(string(plot_script_path, "/", plot_script_name))



# %%
ModelPlots.contour_tlabels




# %%
# %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

zvar = :mbr
    ft = :rt
    ft_zvar = ft == :st ? Symbol(:s_, zvar) : Symbol(:r_, zvar)
    eq_type = :ep_market
    zvar_diff = false
    diff_percent=false

    # s_fi_fv and s_otc_fv generated in the init file.
    # plot_title = ModelPlots.get_final_contour_plot_title(poold[:df], zvar, ft)
    plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                                        xvar=iso_xlabel, yvar=iso_ylabel,
                                                        diff=zvar_diff, diff_percent=diff_percent,
                                                        k_otc=k_otc,
                                                        s_fi_fv=s_fi_fv,
                                                        s_otc_fv=s_otc_fv)

    # function will back out the inverse coupon ratio and the misrepresenting type from the dataframe.
    file_path_name = ModelPlots.get_nrmp_contour_plot_path_name(poold[:df], ft_zvar;
                                                                xvar=iso_xlabel, yvar=iso_ylabel,
                                                                eq_type=eq_type)

    ep_r_mbr_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[ft_zvar],
                                                                 pltd[:bool_Z];
                                                                 sup_title_sep=true,
                                                                 iso_xlabel=iso_xlabel,
                                                                 iso_ylabel=iso_ylabel,
                                                                 fig_title=plot_title,
                                                                 file_path_name=file_path_name,
                                                                 iso_cmap=iso_cmap,
                                                                 iso_levels=iso_levels,
                                                                 cat_cmap=cat_cmap,
                                                                 fig_aspect=.4)


™# %%
PyPlot.figure(constrained_layout=false)

# %%
ModelPlots.iso_plt_inputs[:fig_aspect]

# %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                                    xvar=iso_xlabel, yvar=iso_ylabel,
                                                    diff=zvar_diff, diff_percent=diff_percent,
                                                    k_otc=k_otc, s_otc_fv=s_otc_fv)

# %%
s_fi_fv

# %%
eq_type

# %%
plot_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, ft_zvar; ft=ft,
                                                    xvar=iso_xlabel, yvar=iso_ylabel,
                                                    diff=zvar_diff, diff_percent=diff_percent,
                                                    s_otc_fv=s_otc_fv,
                                                    fig_aspect=.5)

# %%
fig_title = ModelPlots.get_nrmp_contour_plot_title(poold[:df], eq_type, zvar; ft=ft,
                                        xvar=iso_xlabel, yvar=iso_ylabel,
                                        diff=zvar_diff, diff_percent=diff_percent)

# %%
    plot_title = ModelPlots.get_final_contour_plot_title(poold[:df], zvar, ft)

# %%
unique(pltd[:bool_OTC_EP])

# %%
misrep_incentive = (mpd[:r_mbr][:zvals] .> fid[:mbr][:zvals][2:end, 2:end])
sep_eq_cond = sepd[:s_fv][:zvals].> poold[:s_fv][:zvals]
(sepd[:s_fv][:zvals]).*(misrep_incentive).*(sep_eq_cond)
# %%
(sepd[:s_fv][:zvals]).*(misrep_incentive).*(sep_eq_cond .==false)
# %%
(poold[:s_fv][:zvals]).*(misrep_incentive).*(sep_eq_cond .==false)

# %%
(sepd[:s_fv][:zvals]).*((sepd[:s_fv][:zvals] .- poold[:s_fv][:zvals]) .> 0)
max.(sepd[:s_fv][:zvals], poold[:s_fv][:zvals])
# %%
s_otc_fv

# %%
(max.(sepd[:s_fv][:zvals], poold[:s_fv][:zvals]) .- s_otc_fv ).*(misrep_incentive).*(sep_eq_cond .==false)

# %%
fid[:df][fid[:df][:, :sigmah].== .15, :]



# %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end


# %% codecell
plot_xvar = :nrm_sigmah
plot_yvar = :nrm_lambda
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
plot_script_name = "nrmp_fixed_mu_s_contour_plots_init.jl"
include(string(plot_script_path, "/", plot_script_name))

# %%
ModelPlots.contour_plots_title_params_order


# %% codecell
plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
plot_script_name = "nrmp_fixed_mu_s_contour_plots.jl"
include(string(plot_script_path, "/", plot_script_name))

# %%


# %%
# # ###################################################################################
# * TRASH #############################################################################
# # ###################################################################################
# %%
xvar = :sigmah
yvar = :mu_s
fixed_vars_vec = [:mu_s, :lambda]
fixed_var = [x for x in fixed_vars if !(x in [xvar, yvar])][1]


# %%
fid[:df][3, :nrm_lambda]

# %%
fid

# %%


# %%
function get_nrmp_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs)
                                         # fi_fv::Float64, fi_fv_fun, k_otc::Float64)
    jeq_ind = (x, y) -> mp_funs[:r_mbr](x, y) >= fi_funs[:mbr](x, y)
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:s_fv](x, y) >= sep_funs[:s_fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false

    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :r_mbr => Dict{Symbol, Any}())

    fun_dict[:r_mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[:mbr](x, y) : .0
    fun_dict[:r_mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[:r_mbr](x, y) : .0
    fun_dict[:r_mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[:r_mbr](x, y) : .0

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) +
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) +
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:r_mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:r_mbr](x, y))

    fi_fv = 100.
    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv + #fi_funs[:fv](x, y) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:s_fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:s_fv](x, y))

    # fun_dict[:bool_otc] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : 0
    # fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : fun_dict[:eq_bool](x, y)
    #
    # catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)],
    #                 [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
    # fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]

    return fun_dict
end


# %%
function get_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs,
                                     fi_fv::Float64, fi_fv_fun, k_otc::Float64)
    jeq_ind = (x, y) -> mp_funs[:mbr](x, y) >= fi_funs[:mbr](x, y)
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:safe][:fv](x, y) >= sep_funs[:safe][:fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false

    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :mbr => Dict{Symbol, Any}())
    for zvar in [:fv, :mbr, :lev]
        fun_dict[:mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[zvar](x, y) : .0
        fun_dict[:mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[:risky][zvar](x, y) : .0
        fun_dict[:mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[:risky][zvar](x, y) : .0
    end

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) +
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) +
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:risky][:mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:risky][:mbr](x, y))
    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv + #fi_funs[:fv](x, y) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:safe][:fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:safe][:fv](x, y))
    fun_dict[:bool_otc] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : 0
    fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(k_otc)) ? 1 : fun_dict[:eq_bool](x, y)

    catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)],
                    [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
    fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]

    return fun_dict
end


# %%
function jeq_payoff_functions(fi_funs::Dict{Symbol, Any}, jfd::Dict;
                              eq_type::String="pooling",
                              xvar::Symbol=:iota, yvar::Symbol=:sigmah)


    if any([fi_funs[:xvar] != xvar, fi_funs[:yvar] != yvar])
        println("Full Info and Misrep x and y variables do not coincide. Exiting...")
        return
    end

    jeq_funs = Dict{Symbol, Any}(:xvar => xvar, :yvar => yvar,
                                 :safe => Dict{Symbol, Any}(),
                                 :risky => Dict{Symbol, Any}())

    r_obj_fun = :MBR
    if eq_type == "separating"
        r_obj_fun = :firm_value
    end


    # Set objective function symbol
    zsym = contour_zvars_sym[r_obj_fun]

    # Risky Firm's Objective Function Payoff in case of Risk-Management v.s. No Risk-Management
    jeq_funs[:risky][Symbol(:rm_, zsym)], jeq_funs[:risky][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar, Symbol(:r_, r_obj_fun))

    # Choose Risk-Management if it maximizes Payoff
    jeq_funs[:risky][:rm_cond] = (x, y) -> jeq_funs[:risky][Symbol(:rm_, zsym)](x, y) >= jeq_funs[:risky][Symbol(:nrm_, zsym)](x, y)

    # Does this RMP differ from the optimal RMP in Full Info Eq?
    jeq_funs[:risky][:jeq_fi_rm_diff] = (x, y) -> (jeq_funs[:risky][:rm_cond](x, y) != fi_funs[:rm_cond](x, y)) ? 1. : 0.

    # Risky Firm's Objective Function Payoff
    jeq_funs[:risky][zsym] = (x, y) -> maximum([jfd[Symbol(:r_, r_obj_fun)][xvar](x),
                                                jfd[Symbol(:r_, r_obj_fun)][yvar](y)])

    # Difference in RMP between Joint and Full Information Equilibrium
    jeq_funs[:risky][Symbol(:jeq_fi_, zsym, :_diff)] = (x, y) -> jeq_funs[:risky][zsym](x, y) - fi_funs[zsym](x, y)

    # Safe Firm's Payoff Depends on what Risky Firm chooses
    jeq_funs[:safe][Symbol(:rm_, zsym)], jeq_funs[:safe][Symbol(:nrm_, zsym)] = get_rm_payoff_funs(jfd, xvar, yvar,
                                                                                                   Symbol(:s_, r_obj_fun))
    jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[:safe][Symbol(:rm_, zsym)](x, y) : jeq_funs[:safe][Symbol(:nrm_, zsym)](x, y)

    # jeq_funs[:safe][zsym] = (x, y) -> jeq_funs[:safe][zsym](x, y) - Need FI payoff #fi_funs[zsym](x, y)

    for zvar in [z for z in contour_zvars if z != r_obj_fun]
        zsym2 = contour_zvars_sym[zvar]

        for ft in keys(contour_firm_types)
            ft_z = contour_firm_types[ft]
            jeq_funs[ft][Symbol(:rm_, zsym2)], jeq_funs[ft][Symbol(:nrm_, zsym2)] = get_rm_payoff_funs(jfd, jfd[:xvar],
                                                                                                     jfd[:yvar],
                                                                                                     Symbol(ft_z, zvar))
            jeq_funs[ft][zsym2] = (x, y) -> jeq_funs[:risky][:rm_cond](x, y) ? jeq_funs[ft][Symbol(:rm_, zsym2)](x, y) : jeq_funs[ft][Symbol(:nrm_, zsym2)](x, y)
        end

        jeq_funs[:risky][Symbol(:jeq_fi_, zsym2, :_diff)] = (x, y) -> jeq_funs[:risky][zsym2](x, y) - fi_funs[zsym2](x, y)
    end


    return deepcopy(jeq_funs)
end
