# vim: set fdm=marker : 

using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using JLD
using CSV
using Dierckx
using Dates
using PyPlot
using Seaborn
using LaTeXStrings

using Revise

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
# include(string(joinpath(module_path, "2period"), ".jl"))
modls = ["2period"] # "ModelPlots",
for modl in modls
    include(string(joinpath(module_path, modl), ".jl"))
end

# Create Firms  {{{1
sf = p2m.firm_initializer(:st, p2m.pardict)
rf = p2m.firm_initializer(:rt, p2m.pardict)

# Inputs {{{1
rqgrid = .0:0.05:1.
mu_s_grid = .1:.1:.8

mu_b_min=1e-3
mubN=30
mu_b_max=NaN
mu_b_max = isnan(mu_b_max) ? sf.V0/sf.D : mu_b_max
mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)

# Illiquidity Cost
ilc_grid = 2*1e-2:5*10^-4:6*1e-2

# Compute OTC Results {{{1
ilcdf_list = fetch(@spawn [p2m.get_otc_prices(sf, ilc; mu_b_grid=Array(mu_b_grid))
                            for ilc in ilc_grid])
ilcdf = vcat(ilcdf_list...)

otcdf_list = fetch(@spawn [p2m.get_otc_opt_lev(sf; ilc=ilc, df=ilcdf) 
                            for ilc in unique(ilcdf[:, :ilc])])

otcdf = vcat(otcdf_list...)

# Reset Firms
sf = p2m.firm_initializer(:st, p2m.pardict)
rf = p2m.firm_initializer(:rt, p2m.pardict)

# Full Information and Misrepresentation  {{{1
df_list = fetch(@spawn [p2m.get_prices(sf, rf; rq=rq) for rq in rqgrid])
df = vcat(df_list...)

fidf_list = fetch(@spawn [p2m.get_opt_lev(sf, rf; df=df, rq=rq) for rq in rqgrid])
fidf = vcat(fidf_list...)

# Pooling Market Outcomes {{{1
# pooldf = p2m.get_pool_res(sf, rf, fidf[:, :rq], mu_s_grid[1]; fidf=fidf)  
pooldf_list = fetch(@spawn [p2m.get_pool_res(sf, rf, fidf[:, :rq], 
                                             mu_s; fidf=fidf)  
                            for mu_s in mu_s_grid])
pooldf = vcat(pooldf_list...)

# Separating Market Outcomes  {{{1
sepdf = p2m.get_sep_res(sf, rf, fidf)

# Form Interpolated Surfaces {{{1

## Restricted q grid {{{2
xvals = unique(sepdf[:, :q])
yvals = Array(mu_s_grid)
fid = p2m.form_eq_dict(fidf; xvals=xvals, yvals=yvals)
poold = p2m.form_eq_dict(pooldf; xvals=xvals, yvals=yvals)
sepd = p2m.form_eq_dict(sepdf; xvals=xvals, yvals=yvals)

mpf = [zf for zf in keys(fid[:interp]) if occursin("mp", String(zf))]
fi_funs = Dict([k => fid[:interp][k] for k in keys(fid[:interp]) if !(k in mpf)])
mp_funs = Dict([k => fid[:interp][k] for k in mpf])
pool_funs = poold[:interp]
sep_funs = sepd[:interp]
fi_fv = fidf[1, :s_fi_fv]

# OTC Function
ilc_otc = otcdf[:, :ilc][39]
s_otc_fv = otcdf[abs.(otcdf[:, :ilc] .- ilc_otc) .< 1e-5, :s_fi_fv][1]
s_otc_mbr = otcdf[abs.(otcdf[:, :ilc] .- ilc_otc) .< 1e-5, :s_fi_mbr][1]


fi_fv_fun = Dierckx.Spline1D(otcdf[:, :ilc], otcdf[:, :s_fi_fv]; 
                             k=3, bc="extrapolate")
fun_dict = p2m.get_contour_equilibria_funs(fi_funs, mp_funs, 
                                       pool_funs, sep_funs,
                                       fi_fv, fi_fv_fun, ilc_otc)
pltd = p2m.get_eq_contour_mesh_grid(xvals, yvals, fun_dict)

# Reset Firms
sf = p2m.firm_initializer(:st, p2m.pardict)
rf = p2m.firm_initializer(:rt, p2m.pardict)

# Full q grid {{{2 
xvals = unique(fidf[:, :rq])
yvals = Array(mu_s_grid)
ext_fid = p2m.form_eq_dict(fidf; xvals=xvals, yvals=yvals)
ext_poold = p2m.form_eq_dict(pooldf; xvals=xvals, yvals=yvals)
sepd = p2m.form_eq_dict(sepdf; xvals=unique(sepdf[:, :q]), yvals=yvals)

ext_mpf = [zf for zf in keys(ext_fid[:interp]) if occursin("mp", String(zf))]
ext_fi_funs = Dict([k => ext_fid[:interp][k] for k in keys(ext_fid[:interp]) if !(k in mpf)])
ext_mp_funs = Dict([k => ext_fid[:interp][k] for k in mpf])
ext_pool_funs = ext_poold[:interp]
ext_sep_funs = Dict{Symbol, Any}([k => nothing for k in keys(sepd[:interp])]) #sepd[:interp]
xmax = maximum(unique(sepdf[:, :q]))
for zfun in keys(ext_sep_funs)
    sep_tmp(x, y) = (x > xmax) ? -999. : sepd[:interp][zfun](x,y)
    ext_sep_funs[zfun] = deepcopy(sep_tmp)
end
fi_fv = fidf[1, :s_fi_fv]

# OTC Function
fi_fv_fun = Dierckx.Spline1D(otcdf[:, :ilc], otcdf[:, :s_fi_fv]; 
                             k=3, bc="extrapolate")

ext_fun_dict = p2m.get_contour_equilibria_funs(ext_fi_funs, ext_mp_funs, 
                                               ext_pool_funs, ext_sep_funs,
                                               fi_fv, fi_fv_fun, ilc_otc)
ext_pltd = p2m.get_eq_contour_mesh_grid(xvals, yvals, ext_fun_dict)

# Reset Firms
sf = p2m.firm_initializer(:st, p2m.pardict)
rf = p2m.firm_initializer(:rt, p2m.pardict)


# Store (q, mu_s)-point results
ptr = Dict{Symbol, Dict{Symbol, Any}}(:p1 => Dict{Symbol, Any}(:q => .3, 
                                                               :mu_s => .7),
                                      :p2 => Dict{Symbol, Any}(:q => .3, 
                                                               :mu_s => .3),
                                      :p3 => Dict{Symbol, Any}(:q => .75, 
                                                               :mu_s => .7),
                                      :p4 => Dict{Symbol, Any}(:q => .75, 
                                                               :mu_s => .3))
for pt in keys(ptr)
    # ptr[pt][:s_fi_fv] = ext_fi_funs[:s_fi_fv](ptr[pt][:q], ptr[pt][:mu_s])
    # ptr[pt][:r_fi_fv] = ext_fi_funs[:r_fi_fv](ptr[pt][:q], ptr[pt][:mu_s])

    ptr[pt][:r_mp_fv] = ext_mp_funs[:r_mp_fv](ptr[pt][:q], ptr[pt][:mu_s])
    ptr[pt][:r_mp_mbr] = ext_mp_funs[:r_mp_mbr](ptr[pt][:q], ptr[pt][:mu_s])

    for ft in [:s, :r]
        for zfun in [:fv, :mbr]
            zvar = Symbol(ft, :_fi_, zfun)
            ptr[pt][zvar] = ext_fi_funs[zvar](ptr[pt][:q], ptr[pt][:mu_s])

            zvar = Symbol(ft, :_pool_, zfun)
            ptr[pt][zvar] = ext_pool_funs[zvar](ptr[pt][:q], ptr[pt][:mu_s])

            zvar = Symbol(ft, :_sep_, zfun)
            ptr[pt][zvar] = ext_sep_funs[zvar](ptr[pt][:q], ptr[pt][:mu_s])
        end
    end
end

ptdf = vcat([DataFrame(ptr[pt]) for pt in keys(ptr)]...)
rptdf = p2m.df_reshaper(ptdf)

resd = Dict(:sf => sf, :rf => rf, :ptr => ptr, 
            :s_otc_fv => s_otc_fv, :s_otc_mbr => s_otc_mbr,
            :ilc_otc => ilc_otc,
            :fid => fid, :ext_fid => ext_fid,
            :poold => poold, :ext_poold => ext_poold,
            :sepd => sepd,
            :pltd => pltd, :ext_pltd => ext_pltd)
