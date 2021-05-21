##
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

# dir_vec = rsplit(pwd(), "/")
# pos = findall(x -> x .== "bond-model", dir_vec)[1]
# main_path = join(dir_vec[1:pos], "/")
main_path = "/home/artur/BondPricing/bond-model"
module_path=string(main_path, "/modules")
modls = ["2period", "ModelPlots"]
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

ENV["LINES"] = 750
ENV["COLUMNS"] = 1000
##

# %%
# include(string(joinpath(module_path, "2period"), ".jl"))
# p2m_init_path = string(main_path, "/Julia/Batch/p2m/")
# include(string(joinpath(p2m_init_path, "p2m_initializer"), ".jl"))

##
p2m_init_path = string(main_path, "/Batch/p2m/")
res = include(string(joinpath(p2m_init_path, "p2m_initializer"), ".jl"))
##

##
keys(res)
##


# %%
include(string(joinpath(p2m_init_path, "p2m_plots"), ".jl"))

##
sf = res[:sf]
rf = res[:rf]
otcdf = res[:otcdf]
fidf = res[:ext_fid][:df]
ilc_otc = res[:ilc_otc]
theta = 1.

tmp = p2m.otc_misrep_mbr(sf, rf, fidf, otcdf, ilc_otc; theta = 1.)
p2m.plot_r_otc_misrep_mbr(rf, tmp)
##


##
row = 100
mus = pooldf[row, :mu_s]
rq = pooldf[row, :q]
pool_mu_b = pooldf[row, :pool_mu_b]
pool_debt = pooldf[row, :s_pool_debt] + (pooldf[row, :s_pool_debt] - pooldf[row, :r_pool_debt])^10
s_pool_eq = p2m.get_epr(sf; q = 0.0, mu_b = pool_mu_b)
r_pool_eq = p2m.get_epr(rf; q=rq, mu_b = pool_mu_b)
s_fv = s_pool_eq + pool_debt
r_fv = r_pool_eq + pool_debt
pooldf[row, :exp_pool_fv] - (mus * s_fv + (1 - mus) * r_fv)
# pooldf[row, :s_pool_fv] - p2m.get_fv(sf; q = 0.0, mu_b = pool_mu_b)
##

# %%
pdff[1, [:mu_s, :q, :s_pool_mbr]]


# %%
qgrid = [.5]
rf = p2m.setq(rf, q=qgrid[1])
mu_s_grid = .015:.015:.975
fidf = p2m.get_opt_lev(sf, rf; )
pdff = vcat([p2m.get_pool_res(sf, rf, qgrid, mu_s; fidf=fidf) for mu_s in mu_s_grid]...)

qr = qgrid[1]
s_pool_fv = res[:poold][:df][abs.(res[:poold][:df][:, :q] .- qr) .< 1e-5, :s_pool_fv][1]
s_pool_mbr = res[:poold][:df][abs.(res[:poold][:df][:, :q] .- qr) .< 1e-5, :s_pool_mbr][1]
s_sep_fv = res[:sepd][:df][abs.(res[:sepd][:df][:, :q] .- qr) .< 1e-5, :s_sep_fv][1]
s_sep_mbr = res[:sepd][:df][abs.(res[:sepd][:df][:, :q] .- qr) .< 1e-5, :s_sep_mbr][1]
println("s_pool_fv: ", s_pool_fv)
println("s_pool_mbr: ", s_pool_mbr)
println("s_sep_fv: ", s_sep_fv)
println("s_sep_mbr: ", s_sep_mbr)

fg = p2m.plot_pool_mub_mbr(pdff)

# %%
res[:sepd][:df][abs.(res[:sepd][:df][:, :q] .- qr) .< 1e-5, :r_sep_mbr][1]

# %%
mu_s = .3
pool_bpr = [mu_s * p2m.get_bpr(sf; mu_b=mu_b) + (1 - mu_s) * p2m.get_bpr(rf; mu_b=mu_b) for mu_b in mu_b_grid]
pool_debt = mu_b_grid .* pool_bpr
# pool_debt_fun = [p2m.get_pool_debt(sf, rf, mu_s; mu_b=mu_b) for mu_b in mu_b_grid]
# maximum(abs.(pool_debt_fun .- pool_debt))
s_pool_fv = [p2m.get_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid]  .+ pool_debt
s_pool_mbr = [p2m.get_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid] ./ (sf.V0 .- pool_debt)
r_pool_fv = [p2m.get_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid]  .+ pool_debt
r_pool_mbr = [p2m.get_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid] ./ (rf.V0 .- pool_debt)


# %%
p2m.get_pool_res(sf, rf, [.5], .3)

# %%
p2m.get_bpr(sf; mu_b=mu_b)


# %%
# Compute FI, MISREP and POOL Results for Shortened List of mub values
rf = p2m.setq(rf; q=.5)
mu_s = .3
mu_b_grid = range(.05, stop=1.6, length=100)
s_fv = fetch(@spawn [p2m.get_fv(sf; mu_b=mu_b) for mu_b in mu_b_grid])
s_mbr = fetch(@spawn [p2m.get_mbr(sf; mu_b=mu_b) for mu_b in mu_b_grid])
r_fv = fetch(@spawn [p2m.get_fv(rf; mu_b=mu_b) for mu_b in mu_b_grid])
r_mbr = fetch(@spawn [p2m.get_mbr(rf; mu_b=mu_b) for mu_b in mu_b_grid])
r_misrep_mbr = fetch(@spawn [p2m.get_misrep_mbr(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid])
# pdf = p2m.get_pool_res(sf, rf, mu_b_grid, mu_s = .3)

pool_bpr = [mu_s * p2m.get_bpr(sf; mu_b=mu_b) + (1 - mu_s) * p2m.get_bpr(rf; mu_b=mu_b) for mu_b in mu_b_grid]
pool_debt = mu_b_grid .* pool_bpr
# pool_debt_fun = [p2m.get_pool_debt(sf, rf, mu_s; mu_b=mu_b) for mu_b in mu_b_grid]
# maximum(abs.(pool_debt_fun .- pool_debt))
s_pool_fv = [p2m.get_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid]  .+ pool_debt
s_pool_mbr = [p2m.get_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid] ./ (sf.V0 .- pool_debt)
r_pool_fv = [p2m.get_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid]  .+ pool_debt
r_pool_mbr = [p2m.get_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid] ./ (rf.V0 .- pool_debt)

# %%
# Interpolate Functions
s_fv_interp = Dierckx.Spline1D(mu_b_grid, s_fv; k=3, bc="extrapolate")
s_mbr_interp = Dierckx.Spline1D(mu_b_grid, s_mbr; k=3, bc="extrapolate")
r_fv_interp = Dierckx.Spline1D(mu_b_grid, r_fv; k=3, bc="extrapolate")
r_mbr_interp = Dierckx.Spline1D(mu_b_grid, r_mbr; k=3, bc="extrapolate")
r_misrep_mbr_interp = Dierckx.Spline1D(mu_b_grid, r_misrep_mbr; k=3, bc="extrapolate")
s_pool_fv_interp = Dierckx.Spline1D(mu_b_grid, s_pool_fv; k=3, bc="extrapolate")
s_pool_mbr_interp = Dierckx.Spline1D(mu_b_grid, s_pool_mbr; k=3, bc="extrapolate")
r_pool_fv_interp = Dierckx.Spline1D(mu_b_grid, r_pool_fv; k=3, bc="extrapolate")
r_pool_mbr_interp = Dierckx.Spline1D(mu_b_grid, r_pool_mbr; k=3, bc="extrapolate")

exp_pool_fv_interp(x) = mu_s * s_pool_fv_interp(x) + (1-mu_s) * r_pool_fv_interp(x)

# Refine Grid and Back-out Results
mu_b_grid_ref = range(minimum(mu_b_grid), stop=maximum(mu_b_grid), length=10^5)
s_fi_fv = maximum([s_fv_interp(mu_b) for mu_b in mu_b_grid_ref])
s_fi_mub = mu_b_grid_ref[argmax([s_fv_interp(mu_b) for mu_b in mu_b_grid_ref])]
s_fi_mbr = s_mbr_interp(s_fi_mub)
r_fi_fv = maximum([r_fv_interp(mu_b) for mu_b in mu_b_grid_ref])
r_fi_mub = mu_b_grid_ref[argmax([r_fv_interp(mu_b) for mu_b in mu_b_grid_ref])]
r_fi_mbr = r_mbr_interp(r_fi_mub)

# %%
pool_mub = mu_b_grid_ref[argmax([exp_pool_fv_interp(mu_b) for mu_b in mu_b_grid_ref])]
s_pool_fv = s_pool_fv_interp(pool_mub)
s_pool_mbr = s_pool_mbr_interp(pool_mub)
r_pool_fv = r_pool_fv_interp(pool_mub)
r_pool_mbr = r_pool_mbr_interp(pool_mub)

# %% Create 2-way plot
# INTERESTING: to discourage pooling, safe type may have to drastically reduce its leverage.
# PyPlot.plot(mu_b_grid_ref, [s_fv_interp(mu_b) for mu_b in mu_b_grid_ref])
PyPlot.plot(mu_b_grid_ref, [r_misrep_mbr_interp(mu_b) for mu_b in mu_b_grid_ref])

# %%


# %%
r_fi_mub

# %%
mu_b_sep = [mu_b for mu_b in mu_b_grid_ref if r_misrep_mbr_interp(mu_b) <= r_fi_mbr]
s_sep_fv = maximum([s_fv_interp(mu_b) for mu_b in mu_b_sep])
s_sep_mub = mu_b_sep[argmax([s_fv_interp(mu_b) for mu_b in mu_b_sep])]


# %%
resd[:fid][:df][abs.(resd[:fid][:df][:, :rq] .- .5) .< 1e-5, :]


# %%
s_sep_fv

# %%
resd[:sepd][:df][abs.(resd[:sepd][:df][:, :q] .- .5) .< 1e-5, :]

# %%
?p2m.get_misrep_mbr(rf; mu_b=mu_b_grid[1])


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
# mu_s = vcat(.0, mu_s_grid, 1.)
mu_b = vcat(fidf[:, :r_fi_mu_b], pdff[:, :pool_mu_b], fidf[:, :s_fi_mu_b])
s_pool_mbr = vcat(fidf[:, :r_fi_mu_b], pdff[:, :pool_mu_b], fidf[:, :s_fi_mu_b])
opt_mub = Dierckx.Spline1D(mu_s, mu_b; k=3, bc="extrapolate")


# %% CONTOUR PLOTS
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
