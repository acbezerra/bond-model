# %%
using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using Dierckx
using Dates
using PyPlot
using Seaborn
using LaTeXStrings
using Revise

dir_vec = rsplit(pwd(), "/")
pos = findall(x -> x .== "bond-model", dir_vec)[1]
main_path = join(dir_vec[1:pos], "/")
module_path=string(main_path, "/modules")
modls = ["2period", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

ENV["LINES"] = 750
ENV["COLUMNS"] = 1000

# %%
dir_vec = rsplit(pwd(), "/")
pos = findall(x -> x .== "bond-model", dir_vec)[1]
main_path = join(dir_vec[1:pos], "/")

# %%
# include(string(joinpath(module_path, "2period"), ".jl"))
# p2m_init_path = string(main_path, "/Julia/Batch/p2m/")
# include(string(joinpath(p2m_init_path, "p2m_initializer"), ".jl"))

# %%
p2m_init_path = string(main_path, "/Batch/p2m/")
res = include(string(joinpath(p2m_init_path, "p2m_initializer"), ".jl"))


# %%
df[:, :s_mbr] = df[!, :s_eq]./(sf.V0 .- df[!, :s_debt])
df[:, :r_mbr] = df[!, :r_eq]./(rf.V0 .- df[!, :r_debt])


# %%
?ax1.plot()

# %%
fig, ax1 = PyPlot.subplots()


# ax1.set_xlabel('time (s)')
# ax1.set_ylabel('exp', color=color)
#
Seaborn.set_style("darkgrid")
s_color = "blue"
r_color = "purple"
s_cond = df[:, :rq] .== .0
r_cond = df[:, :rq] .== .4
PyPlot.plot(df[s_cond, :mu_b], df[s_cond, :s_fv], color=s_color)
PyPlot.plot(df[r_cond, :mu_b], df[r_cond, :r_fv], color=r_color)
# # ax1.tick_params(axis='y', labelcolor=color)
#
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
color = "red"
# ax2.set_ylabel('sin', color=color)  # we already handled the x-label with ax1
PyPlot.plot(df[s_cond, :mu_b], df[s_cond, :s_mbr], linestyle="-.", color=s_color)
PyPlot.plot(df[r_cond, :mu_b], df[r_cond, :r_mbr], linestyle="-.", color=r_color)
# ax2.tick_params(axis='y', labelcolor=color)
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
PyPlot.show()

# %%
mu_s = .7
pgrid = range(1e-3, 1.; length=10)
smuhat = (mu_s, p_s, p_r) -> (p_s * mu_s)/(p_s * mu_s + p_r * (1 - mu_s))
poolbpr(mu_s, p_s, p_r) = p2m.get_pool_bpr(sf, rf, smuhat(mu_s, p_s, p_r); rq=.5)

# %%
maximum([smuhat(mu_s, p_s, p_r) for p_s in pgrid, p_r in pgrid])
minimum([smuhat(mu_s, p_s, p_r) for p_s in pgrid, p_r in pgrid])

# %%
[smuhat(mu_s, p_s, p_r) for p_s in pgrid, p_r in pgrid]

# %%
p2m.setq(rf, q=.3)




# %% POOLING -> adjust get_pool_res to maximize joint firm value
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end
qgrid = [.5]
rf.q = qgrid[1]
mu_s_grid = .015:.015:.975
fidf = p2m.get_opt_lev(sf, rf; )
pdff = vcat([p2m.get_pool_res(sf, rf, qgrid, mu_s; fidf=fidf) for mu_s in mu_s_grid]...)

# %%
# mu_s = vcat(.0, mu_s_grid, 1.)
mu_b = vcat(fidf[:, :r_fi_mu_b], pdff[:, :pool_mu_b], fidf[:, :s_fi_mu_b])
s_pool_mbr = vcat(fidf[:, :r_fi_mu_b], pdff[:, :pool_mu_b], fidf[:, :s_fi_mu_b])
opt_mub = Dierckx.Spline1D(mu_s, mu_b; k=3, bc="extrapolate")

# %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end
fg = p2m.plot_pool_mub_mbr(pdff)


# %%
findfirst("bond-model", fpath)[end]



# %%




# %%
[split(pwd(), "/") .== "bond-model"]

#  %%
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end
plot_follder = "/home/artur/BondPricing/bond-model/Plots/M2P/"


file_path_name = p2m.get_2pm_contour_plot_path_name(:mbr)

# %%
# plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z;
#                                             ft=Symbol(ft, :t))

                                            # %%
fig.suptitle("Optimal Pooling Capital Structure and Shareholders' Payoff \n ",
             L"$q_s=0.0$, $q_r=", pdff[1, :q])

# %%
string(L"$q_s=0.0$, $q_r=", pdff[1, :q], L"$")

# %%
println(L"$q_s=0.0$, $q_r= df[1, :q]")

# %%
println(L"$\pi$", " \n", "dgfkhjd")

# %%
sf.sig

# %%
sf.alpha

# %%

title_eq_type = [eq_type_2pm_title[k][2] for k in keys(eq_type_2pm_title)
                 if (eq_type_2pm_title[k][1] == eq_type)][1]

# %%
sf.D


# %%
p2m.eq_type_2pm_title

# %%
p2m.q_pool_res(sf, rf, .2; fidf=fidf)

# %%
p2m.get_bpr(sf)

# %%
p2m.get_pool_bpr(sf, rf, .3)

# %%
p2m.get_pool_bpr(sf, rf, smuhat(mu_s, .3, .4); rq=.5)

# %%
bpr1 = [poolbpr(mu_s, p_s, p_r) for p_s in pgrid, p_r in pgrid]

# %%
sfv = Dierckx.Spline1D(df[s_cond, :mu_b], df[s_cond, :s_fv], k=3, bc="extrapolate")
mu_b_grid = .0:.005:1.1
s_mu_b_opt = mu_b_grid[argmax(sfv(mu_b_grid))]
sfv(s_mu_b_opt)

# %%
s_mu_b_opt

# %%
rfv = Dierckx.Spline1D(df[r_cond, :mu_b], df[r_cond, :r_fv], k=3, bc="extrapolate")
mu_b_grid = .0:.005:1.1
r_mu_b_opt = mu_b_grid[argmax(rfv(mu_b_grid))]
rfv(r_mu_b_opt)

# %%
r_mu_b_opt

# %%
mu_b_grid[argmax(sfv(mu_b_grid))]

# %%
fidf[1, :]

# %%
s_cond = df[:, :rq] .== .0
PyPlot.plot(df[s_cond, :mu_b], df[s_cond, :s_fv])

# %%
PyPlot.plot(df[s_cond, :mu_b], df[s_cond, :s_mbr])


# %%
names(df[:, :])

# %%
tmp2 = df[df[:, :rq] .== 0.3, [:mu_b, :r_debt, :r_eq, :r_fv]]
tmp2[!, :r_mbr] = tmp2[:, :r_eq]./(sf.V0 .- tmp2[:, :r_debt])

# %%
PyPlot.plot(tmp2[!, :mu_b], tmp2[!, :r_fv])

# %%
PyPlot.plot(tmp2[!, :mu_b], tmp2[!, :r_mbr])

# %%
fidf_list = fetch(@spawn [p2m.get_opt_lev(sf, rf; df=df, rq=rq) for rq in rqgrid])
fidf = vcat(fidf_list...)
