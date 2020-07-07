# vim: set fdm=marker : 

module p2m

using Distributions
using Distributed
using DataFrames
using Dierckx
using Seaborn
using Printf
using LaTeXStrings

# Inputs {{{1
pardict=Dict{Symbol, Float64}(:mu_s => .2,
                              :V0 => 100.,
                              :Delta => 80.,
                              :D => 90.,
                              :pi => .27,
                              :alpha => .4,
                              :mu_factor => 1.,
                              :rf=>0.08,
                              :rbdisc=>0.1,
                              :qs => 0.,
                              :qr => .3,
                              :sigma=>.3,
                              :shock_factor=>.5)

zfuns = [:bpr, :yield, :debt, :eq, :fv, :mbr]

# Model Objects {{{1

  # nd = Distributions.Normal(0,1)
mutable struct Rates
    rf::Float64
    rbdisc::Float64
end

mutable struct Firm
   ft::Symbol

   mu_b::Float64
   V0::Float64
   D::Float64
   alpha::Float64
   pi::Float64

   q::Float64
   muxh::Float64
   muxl::Float64
   sig::Float64
   nd::Normal{Float64}
   rt::Rates
end

function firm_initializer(ft::Symbol, pardict::Dict{Symbol, Float64}; mu_b::Float64=0.)
        
    nd = Distributions.Normal(0., 1)
    rt = Rates(pardict[:rf], pardict[:rbdisc])
   
    muxh = pardict[:mu_factor] * rt.rf - .5 * pardict[:sigma]^2
    smuxl = muxh
    rmuxl = muxh - pardict[:shock_factor] * pardict[:sigma]


    q = ft == :st ? pardict[:qs] : pardict[:qr]
    muxl = ft == :st ? smuxl : rmuxl
    
    return Firm(ft, mu_b, pardict[:V0], 
                pardict[:D], pardict[:alpha],
                pardict[:pi], q, muxh, muxl, 
                pardict[:sigma], nd, rt)
end


# Auxiliary Functions {{{1
function aux_f(fr, mu_b::Float64, shock::Bool)
  mux = shock == false ? fr.muxh : fr.muxl
  return (1/fr.sig) * log(((1 - fr.pi) * mu_b * fr.D)/(fr.V0 * exp(mux)))
end


function df_add_pars(fr, df)
    for par in [:alpha, :pi, :D, :V0, :sig, :muxl, :muxh,]
        df[!, par] .= getfield(fr, par)
    end
    
    df[!, :rf] .= fr.rt.rf
    df[!, :rbdisc] .= fr.rt.rbdisc
    df[!, :ilc] = df[:, :rbdisc] - df[:, :rf]
    
    return df
end


function df_reshaper(df::DataFrame; 
                     cols::Array{Symbol,1}=[:mu_s, :q, :ft, :zf, :fi, :mp, :pool, :sep])

    rdf = DataFrame()
    for row in 1:size(df, 1)
        tmp = DataFrame(df[row, :])

        for ft in [:s, :r]
            for zf in [:fv, :mbr]
                tmpd = Dict(:mu_s => tmp[1, :mu_s], :ft => ft, :q => tmp[1, :q], :zf => zf)
                for eq in [:fi, :mp, :pool, :sep]
                    zvar = Symbol(ft, :_, eq, :_, zf)

                    tmpd[eq] = (zvar in Symbol.(names(tmp))) ? tmp[1, zvar] : NaN
                end

                rdf = vcat(rdf, DataFrame(tmpd))
            end
        end
    end

    for row in 1:size(rdf,1) 
        for col in [:fi, :mp, :pool, :sep]
            rdf[row, col] = parse(Float64, @sprintf("%.2f", rdf[row, col]))
        end
    end
        
    return sort(rdf[:, cols], [:mu_s, :q], rev=true)
end

# Analytic Functions {{{1
# Post-Shock Arrival Time Functions {{{2
function ps_bpr(fr, mu_b::Float64; shock::Bool=false)
  mux = shock == false ? fr.muxh : fr.muxl
  if mu_b < 1e-5
    return 0.
  end

  return (exp(-fr.rt.rbdisc) * (fr.D * cdf(fr.nd, - aux_f(fr, mu_b, shock)) +
           (fr.alpha * fr.V0/mu_b) * exp(mux + .5 * fr.sig^2) * cdf(fr.nd, aux_f(fr, mu_b, shock) - fr.sig)))
end

function ps_debtpr(fr, mu_b::Float64; shock::Bool=false)
  return mu_b * ps_bpr(fr, mu_b; shock=shock)
end

function ps_eqpr(fr, mu_b::Float64; shock::Bool=false) 
  mux = shock == false ? fr.muxh : fr.muxl

  return exp(-fr.rt.rf) * ((exp(mux + .5 * fr.sig^2) * fr.V0 * cdf(fr.nd, - aux_f(fr, mu_b, shock) + fr.sig))
                           - (1 - fr.pi) * mu_b * fr.D * cdf(fr.nd, - aux_f(fr, mu_b, shock)))
end

function ps_fv(fr, mu_b::Float64; shock::Bool=false)
    return ps_debt(fr, mu_b; shock=shock) + ps_eqpr(fr, mu_b; shock=shock)
end

# Pre-Shock Arrival Time Functions {{{2
function bpr(fr, mu_b::Float64)
  return (1 - fr.q) * ps_bpr(fr, mu_b; shock=false) + fr.q * ps_bpr(fr, mu_b; shock=true)
end

function debtpr(fr, mu_b::Float64)
  return (1 - fr.q) * ps_debtpr(fr, mu_b; shock=false) + fr.q * ps_debtpr(fr, mu_b; shock=true)
end

function eqpr(fr, mu_b::Float64)
  return (1 - fr.q) * ps_eqpr(fr, mu_b; shock=false) + fr.q * ps_eqpr(fr, mu_b; shock=true)
end

function fv(fr, mu_b::Float64)
  return (1 - fr.q) * ps_fv(fr, mu_b; shock=false) + fr.q * ps_fv(fr, mu_b; shock=true)
end

function mbr(fr, mu_b::Float64)
    return eqpr(fr, mu_b) / (fr.V0 - debtpr(fr, mu_b))
end


# Get Methods {{{2
function get_bpr(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return bpr(fr, fr.mu_b)
end

function get_yield(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    yd = any([isnan(fr.mu_b), abs.(fr.mu_b) <1e-5]) ? NaN : - log(bpr(fr, fr.mu_b)/fr.D)
    return yd
end

function get_dpr(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return debtpr(fr, fr.mu_b)
end

function get_epr(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return eqpr(fr, fr.mu_b) 
end

function get_fv(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return (debtpr(fr, fr.mu_b) + eqpr(fr, fr.mu_b))
end

function get_mbr(fr; q::Float64=NaN, mu_b::Float64=NaN) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return mbr(fr, fr.mu_b)
end

function get_misrep_fv(sf, rf; 
                       sq::Float64=NaN, rq::Float64=NaN, 
                       mu_b::Float64=NaN)
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)
    
    
    rE1 = get_epr(rf; mu_b=mu_b)
    sD1 = get_dpr(sf; mu_b=mu_b)
    
    return  rE1 + sD1
end

function get_misrep_mbr(sf, rf; 
                        sq::Float64=NaN, rq::Float64=NaN, 
                        mu_b::Float64=NaN)
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)
    
    
    rE1 = get_epr(rf; mu_b=mu_b)
    sD1 = get_dpr(sf; mu_b=mu_b)
    
    return  rE1 / (rf.V0 - sD1)
end


function get_pool_bpr(sf, rf, mu_s::Float64; 
                      sq::Float64=NaN, rq::Float64=NaN,
                      mu_b::Float64=NaN)
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)

    return mu_s * get_bpr(sf) + (1 - mu_s) * get_bpr(rf)
end


function get_pool_debt(sf, rf, mu_s::Float64; 
                       sq::Float64=NaN, rq::Float64=NaN,
                       mu_b::Float64=NaN)
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)

    return mu_s * get_dpr(sf) + (1 - mu_s) * get_dpr(rf)
end


function get_pool_yield(sf, rf, mu_s::Float64; 
                        sq::Float64=NaN, rq::Float64=NaN,
                        mu_b::Float64=NaN)
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)

    if any([isnan(sf.mu_b), abs.(sf.mu_b) < 1e-5])
      return NaN
    end
        
    return - log(get_pool_bpr(sf, rf, mu_s)/sf.D)
end


function get_fi_prices(fr; ft::Symbol=Symbol(""),
                    q::Float64=NaN,
                    mu_b_grid::Array{Float64,1}=Array{Float64,1}(),
                    mu_b_min::Float64=1e-3, 
                    mu_b_max::Float64=NaN, mubN::Int64=30)
 
  fr = setq(fr; q=q)   
  ft = isempty(String(ft)) ? (abs(fr.q) > 1e-3 ? :r : :s) : ft

  if isempty(mu_b_grid)  
    mu_b_max = isnan(mu_b_max) ? fr.V0/fr.D : mu_b_max
    mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)
  end

  bpr_vec = [get_bpr(fr; mu_b=mu_b) for mu_b in mu_b_grid]
  yd_vec = [get_yield(fr; mu_b=mu_b) for mu_b in mu_b_grid]
  dpr_vec = [get_dpr(fr; mu_b=mu_b) for mu_b in mu_b_grid]
  epr_vec = [get_epr(fr; mu_b=mu_b) for mu_b in mu_b_grid]
  fv_vec = [get_fv(fr; mu_b=mu_b) for mu_b in mu_b_grid]
  mbr_vec = [get_mbr(fr; mu_b=mu_b) for mu_b in mu_b_grid]


  return DataFrame(Symbol(ft, :q) => fr.q,
                   :mu_b => mu_b_grid, 
                   Symbol(ft, :_bpr) => bpr_vec,
                   Symbol(ft, :_yield) => yd_vec,  
                   Symbol(ft, :_debt) => dpr_vec, 
                   Symbol(ft, :_eq) => epr_vec, 
                   Symbol(ft, :_fv) => fv_vec, 
                   Symbol(ft, :_mbr) => mbr_vec)
end


function get_otc_prices(fr, ilc;
                        mu_b_grid::Array{Float64,1}=Array{Float64,1}(),
                        mu_b_min::Float64=1e-3, 
                        mu_b_max::Float64=NaN, mubN::Int64=30)
    
    if isempty(mu_b_grid)  
        mu_b_max = isnan(mu_b_max) ? fr.V0/fr.D : mu_b_max
        mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)
    end

    # Set Bond Investor's Discount Rate
    fr = set_ilc(fr, ilc)
    
    ft = abs(fr.q) > 1e-3 ? :r : :s 
    
    otcdf = get_fi_prices(fr; ft=ft, mu_b_grid=Array(mu_b_grid))
    
    return df_add_pars(fr, otcdf)
end


function get_prices(sf, rf; 
                    sq::Float64=NaN, rq::Float64=NaN, 
                    mu_b_min::Float64=1e-3, 
                    mu_b_max::Float64=NaN, mubN::Int64=30)
    sf = setq(sf; q=sq)
    rf = setq(rf; q=rq)

    mu_b_max = isnan(mu_b_max) ? max(sf.V0/sf.D, rf.V0/rf.D) : mu_b_max
    mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)

  
    sfid = get_fi_prices(sf; ft=:s, mu_b_grid=Array(mu_b_grid))
    rfid = get_fi_prices(rf; ft=:r, mu_b_grid=Array(mu_b_grid))

    df = innerjoin(sfid, rfid, on=:mu_b)

    df[!, :r_mp_fv] = [get_misrep_fv(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid]
    df[!, :r_mp_mbr] = [get_misrep_mbr(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid]
        
    return df
end


#= function get_prices(sf, rf; =# 
#=                     sq::Float64=NaN, rq::Float64=NaN, =# 
#=                     mu_b_min::Float64=1e-3, =# 
#=                     mu_b_max::Float64=NaN, mubN::Int64=30) =#
#=     sf = setq(sf; q=sq) =#
#=     rf = setq(rf; q=rq) =#

#=     mu_b_max = isnan(mu_b_max) ? max(sf.V0/sf.D, rf.V0/rf.D) : mu_b_max =#
#=     mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN) =#

#=     sbpr_vec = [get_bpr(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     syd_vec = [get_yield(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     sdpr_vec = [get_dpr(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     sepr_vec = [get_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     sfv_vec = [get_fv(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     smbr_vec = [get_mbr(sf; mu_b=mu_b) for mu_b in mu_b_grid] =#

#=     sfid = get_fi_prices(sf; ft=:s, mu_b_grid=Array(mu_b_grid)) =#
#=     rfid = get_fi_prices(rf; ft=:r, mu_b_grid=Array(mu_b_grid)) =#

#=     rbpr_vec = [get_bpr(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#   
#=     ryd_vec = [get_yield(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     rdpr_vec = [get_dpr(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     repr_vec = [get_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     rfv_vec = [get_fv(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     rmbr_vec = [get_mbr(rf; mu_b=mu_b) for mu_b in mu_b_grid] =#

#=     r_misrep_fv_vec = [get_misrep_fv(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid] =#
#=     r_misrep_mbr_vec = [get_misrep_mbr(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid] =#

#=     df = DataFrame(:sq => sf.q, :rq => rf.q, =#
#=                    :mu_b => mu_b_grid, =# 
#=                    :s_bpr => sbpr_vec, =#
#=                    :s_yield => syd_vec, =#  
#=                    :s_debt => sdpr_vec, =# 
#=                    :s_eq => sepr_vec, =# 
#=                    :s_fv => sfv_vec, :s_mbr => smbr_vec, =#
#=                    :r_bpr => rbpr_vec, =#
#=                    :r_yield => ryd_vec, =#
#=                    :r_debt => rdpr_vec, =# 
#=                    :r_eq => repr_vec, =# 
#=                    :r_fv => rfv_vec, :r_mbr => rmbr_vec, =# 
#=                    :r_mp_fv => r_misrep_fv_vec, =#
#=                    :r_mp_mbr => r_misrep_mbr_vec) =#
    
#=     return df =#
#= end =#

function get_opt_lev(sf, rf; df::DataFrame=DataFrame(), 
                     sq::Float64=NaN, rq::Float64=NaN,
                     zfuns::Array{Symbol, 1}=zfuns,
                     mu_b_min::Float64=1e-3, 
                     mu_b_max::Float64=NaN, mubN::Int64=30,
                     obj_fun::Symbol=:fv, mubN_ref::Int64=10^5)

    sq = isnan(sq) ? sf.q : sq
    rq = isnan(rq) ? rf.q : rq

    if isempty(df)
      df = get_prices(sf, rf; sq=sq, rq=rq,
                      mu_b_min=mu_b_min, mu_b_max=mu_b_max,
                      mubN=mubN)
    end
    
    # Filter DataFrame
    cond = .&(abs.(df[:, :sq] .- sq) .< 1e-5, 
              abs.(df[:, :rq] .- rq) .< 1e-5)
    tmp = df[cond, :]
    
    # Refine mu_b Grid
    mu_b_grid_ref = range(minimum(tmp[:, :mu_b]), stop=maximum(tmp[:, :mu_b]), length=mubN_ref)
   
    fidf = DataFrame(:sq => sq, :rq => rq)
    #= mpdf = DataFrame(:sq => sq, :rq => rq) =#
    s_opt_mu_b = NaN
    for ft in [:s_, :r_]
        # Interpolate Objective Function
        objf = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, Symbol(ft, obj_fun)]; k=3)
        opt_mu_b = mu_b_grid_ref[argmax([objf(mub) for mub in mu_b_grid_ref])]
        col = Symbol(ft, :fi_mu_b)
        fidf[!, col] .= opt_mu_b 
        s_opt_mu_b = ft == :s_ ? opt_mu_b : s_opt_mu_b

        # Interpolate Prices 
        prvec = zfuns

        # Full Information Payoffs
        for fun in prvec
          col = Symbol(ft, :fi_, fun)
          sitp = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, Symbol(ft, fun)]; k=3)
          fidf[!, col] .= sitp(opt_mu_b)
        end

        # Misrepresentation
        if ft == :r_
          for fun in [:mp_fv, :mp_mbr]
            col = Symbol(ft, fun)
            sitp = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, Symbol(ft, fun)]; k=3)
            fidf[!, col] .= sitp(s_opt_mu_b) #opt_mu_b)
          end
        end
    end

    for z in [:fv, :mbr]
      zdf = Symbol(:r_mp_fi_, z, :_diff)
      fidf[!, zdf] = fidf[:, Symbol(:r_mp_, z)] .- fidf[:, Symbol(:r_fi_, z)]

      zpdf = Symbol(:r_mp_fi_, z, :_perc_diff)
      fidf[!, zpdf] = (fidf[!, zdf] ./ fidf[:, Symbol(:r_fi_, z)]) .* 100.
    end
   
    fidf[!, :eq_type] .= :fi

    return fidf 
end


function get_otc_opt_lev(fr; ilc::Float64=NaN,
                         df::DataFrame=DataFrame(), q::Float64=NaN,
                         zfuns::Array{Symbol, 1}=p2m.zfuns,
                         mu_b_grid::Array{Float64,1}=Array{Float64,1}(),
                         mu_b_min::Float64=1e-3, 
                         mu_b_max::Float64=NaN, mubN::Int64=30,
                         obj_fun::Symbol=:fv, mubN_ref::Int64=10^5)
    
    # Set Bond Investor's Discount Rate
    if !isnan(ilc)
        fr = set_ilc(fr, ilc)
    else
        ilc = fr.rt.rbdisc - fr.rt.rf
    end
    
    if isempty(mu_b_grid)
        mu_b_max = isnan(mu_b_max) ? fr.V0/fr.D : mu_b_max
        mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)
    end
        
    q = isnan(q) ? fr.q : q
    
    if isempty(df)
        df = get_otc_prices(fr, ilc; mu_b_grid=Array(mu_b_grid))
    end
    
    # Filter DataFrame
    cond = abs.(df[:, :ilc] .- ilc) .< 1e-5
    tmp = df[cond, :]
        
    # Refine mu_b Grid
    mu_b_grid_ref = range(minimum(tmp[:, :mu_b]), 
                          stop=maximum(tmp[:, :mu_b]), length=mubN_ref)

    otcdf = DataFrame(:eq_type => :otc, :ilc => ilc)

    ft = abs(fr.q) > 1e-3 ? :r_ : :s_
    # Interpolate Objective Function
    objf = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, Symbol(ft, obj_fun)]; k=3)
    opt_mu_b = mu_b_grid_ref[argmax([objf(mub) for mub in mu_b_grid_ref])]
    col = Symbol(ft, :fi_mu_b)
    otcdf[!, col] .= opt_mu_b 

    # Interpolate Prices 
    prvec = zfuns

    # Full Information Payoffs
    for fun in prvec
      col = Symbol(ft, :fi_, fun)
      sitp = Dierckx.Spline1D(tmp[:, :mu_b], tmp[:, Symbol(ft, fun)]; k=3)
      otcdf[!, col] .= sitp(opt_mu_b)
    end

    return df_add_pars(fr, otcdf)
end


# Set Methods  {{{1
# Set Shock Probability
function setq(fr; q::Float64=NaN)
    if !isnan(q)
        setfield!(fr, :q, q)
    end
    
    return fr
end

# Set Measure of Outstanding Bonds
function set_mu_b(fr; mu_b::Float64=NaN)
    if !isnan(mu_b)
        setfield!(fr, :mu_b, mu_b)
    end
    
    return fr
end

# Set Secondary Market Illiquidity 
function set_ilc(fr, ilc::Float64)
  if !isnan(ilc)
    setfield!(fr.rt, :rbdisc, fr.rt.rf + ilc)
  end
    
  return fr
end


# Full Information & Misrepresentation {{{1
function get_fi_misrep_res(df::DataFrame, sf, rf)
    rq_vals = unique(df[:, :rq])

    misrep_cols = [:mp_fv, :mp_mbr]
    cols = vcat([:mu_b, :debt, :equity, :fv, :mbr], misrep_cols...)
    rcols(x) = occursin("mp", string(x)) ? Symbol(:r_, x) : Symbol(:r_fi_, x)
    coln(x) = x == :mu_b ? x : rcols(x)
    
    df2 = DataFrame(:eq_type => :fi)
    for rq in rq_vals
        cond = abs.(df[:, :rq] .- rq) .< 1e-5
        
        tmpd = Dict([coln(col) => df[cond, Symbol(:r_opt_, col)][1] for col in cols])
        tmpd[:q] = rq
        df2 = vcat(df2, DataFrame(tmpd))
    end

    for col in [col for col in cols if !(col in misrep_cols)]
        df2[!, Symbol(:s_fi_, col)] .= df[1, Symbol(:s_opt_, col)]
    end
    
    common_cols = [:q, :mu_b]
    s_cols = [Symbol(:s_fi_, col) for col in cols if !(col in vcat(common_cols, misrep_cols))]
    r_cols = [rcols(col) for col in cols if !(col in common_cols)]
    cols2 = vcat(:eq_type, common_cols, s_cols, r_cols)

    return df2[:, cols2]
end
  
# JEQ Functions {{{1

# Pooling Functions {{{2
function mu_b_interp_vals(df, objf::Symbol, mu_b_grid)    
    # Interpolate Functions
    funs = [y for y in Symbol.(names(df)) if !(y in [:mu_b, :q])]
    fd = Dict()
    for fun in funs
        fd[fun] = Dierckx.Spline1D(df[:, :mu_b], df[:, fun]; k=3)
    end
    
    # Optimal Leverage
    opt_mu_b = mu_b_grid[argmax(fd[objf](mu_b_grid))]
    dt = Dict([fun => fd[fun](opt_mu_b) for fun in funs])
    dt[:pool_mu_b] = opt_mu_b
    dt[:q] = unique(df[:, :q])[1]
    
    return DataFrame(dt)
end


function q_pool_res(sf, rf, mu_s::Float64; 
                    fidf::DataFrame=DataFrame(),
                    mu_b_min::Float64=1e-3, mu_b_max::Float64=NaN, 
                    mubN::Int64=15, mubN_ref::Int64=10^4)
    mu_b_max = isnan(mu_b_max) ? maximum([sf.V0/sf.D, rf.V0/rf.D]) : mu_b_max
    mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)

    fd = Dict([Symbol(ft, :_fi_, zfun) =>  NaN 
                for ft in [:s, :r], zfun in p2m.zfuns])    
    if !isempty(fidf)
       cond = abs.(fidf[:, :rq] .- rf.q) .< 1e-5
       fd = Dict([Symbol(ft, :_fi_, zfun) =>  fidf[cond, Symbol(ft, :_fi_, zfun)][1] 
                 for ft in [:s, :r], zfun in zfuns])
    end
    
    mubdf = DataFrame()
    for mu_b in mu_b_grid
        sf = set_mu_b(sf, mu_b=mu_b)
        rf = set_mu_b(rf, mu_b=mu_b)

        #= pool_debt = mu_s * get_dpr(sf) + (1 - mu_s) * get_dpr(rf) =#
        pool_bpr = get_pool_bpr(sf, rf, mu_s)
        pool_yield = get_pool_yield(sf, rf, mu_s)
        pool_debt = get_pool_debt(sf, rf, mu_s)

        s_eq = get_epr(sf)
        s_fv = s_eq + pool_debt 
        s_mbr = s_eq / (sf.V0 - pool_debt) 

        r_eq = get_epr(rf)
        r_fv = r_eq + pool_debt
        r_mbr = r_eq / (rf.V0 - pool_debt)


        tmpd = Dict{Symbol, Float64}(:mu_s => mu_s,
                                     :q => rf.q,
                                     :mu_b => mu_b,
                                     :s_pool_bpr => pool_bpr,
                                     :s_pool_yield => pool_yield,
                                     :s_pool_debt => pool_debt,
                                     :s_pool_eq => s_eq,
                                     :s_pool_fv => s_fv,
                                     :s_pool_mbr => s_mbr,
                                     :r_pool_bpr => pool_bpr,
                                     :r_pool_yield => pool_yield,
                                     :r_pool_debt => pool_debt,
                                     :r_pool_eq => r_eq,
                                     :r_pool_fv => r_fv,
                                     :r_pool_mbr => r_mbr)

        mubdf = vcat(mubdf, DataFrame(tmpd))
    end
    
    mu_b_grid_ref = range(mu_b_min, stop=mu_b_max, length=mubN_ref)
    pooldf = mu_b_interp_vals(mubdf, :s_pool_fv, mu_b_grid_ref)

    # Compute POOL-FI Differences
    for ft in [:s, :r]
      for zfun in zfuns
          zdiff = NaN
          perc_zdiff = NaN
          if !isnan(fd[Symbol(ft, :_fi_, zfun)])
             pool_col = Symbol(ft, :_pool_, zfun)
             zdiff = pooldf[:, pool_col] .- fd[Symbol(ft, :_fi_, zfun)]
             perc_zdiff = (zdiff./fd[Symbol(ft, :_fi_, zfun)]) .* 100.
          end

          col_prefix = Symbol(ft, :_pool_fi_, zfun)
          pooldf[!, Symbol(col_prefix, :_diff)] .= zdiff
          pooldf[!, Symbol(col_prefix, :_perc_diff)] .= perc_zdiff
      end
    end
    pooldf[!, :eq_type] .= :pool

    return pooldf
end

function get_pool_res(sf, rf, qgrid, mu_s::Float64;
                      fidf::DataFrame=DataFrame(),
                      sq::Float64=.0)
    sf = setq(sf, q=sq)
    pooldf = DataFrame()
    for rq in qgrid
        rf = setq(rf, q=rq)
        dt = q_pool_res(sf, rf, mu_s; fidf=fidf)
        pooldf = vcat(pooldf, dt)
    end
    
    return pooldf
end

# Separating Functions {{{2
function q_sep_res(fidf::DataFrame, sf, rf, rq::Float64; 
                   mu_b_min::Float64=1e-3, mu_b_max::Float64=NaN, 
                   mubN::Int64=15, mubN_ref::Int64=10^4)

    sf = setq(sf, q=.0)
    rf = setq(rf, q=rq)
    
    mu_b_max = isnan(mu_b_max) ? maximum([sf.V0/sf.D, rf.V0/rf.D]) : mu_b_max
    mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)
    mu_b_grid_ref = range(mu_b_min, stop=mu_b_max, length=mubN_ref)
     
    cond = abs.(fidf[:, :rq] .- rf.q) .< 1e-5
    r_fi_mbr = fidf[cond, :r_fi_mbr][1]
    
    mubdf = DataFrame()
    for mu_b in mu_b_grid
        sf = set_mu_b(sf, mu_b=mu_b)
        rf = set_mu_b(rf, mu_b=mu_b)

        dt = Dict{Symbol, Float64}(:mu_b => mu_b,
                                   :s_bpr => get_bpr(sf),
                                   :s_yield => get_yield(sf),
                                   :s_debt => get_dpr(sf),
                                   :s_eq => get_epr(sf), 
                                   :s_fv => get_fv(sf),
                                   :s_mbr => get_mbr(sf),
                                   :r_mp_mbr => get_misrep_mbr(sf, rf))
        
        mubdf = vcat(mubdf, DataFrame(dt))
    end
    
    fd = Dict()
    for fun in [key for key in Symbol.(names(mubdf)) if key != :mu_b]
        fd[fun] = Dierckx.Spline1D(mubdf[:, :mu_b], mubdf[:, fun]; k=3)
    end
    
    # Find separating mu_b value
    r_misrep_mbr_vec = fd[:r_mp_mbr](mu_b_grid_ref)
    mu_b_vec = mu_b_grid_ref[(r_misrep_mbr_vec .<= r_fi_mbr)] 
    s_sep_mu_b = mu_b_vec[argmax(fd[:s_fv](mu_b_vec))]


    sepd = Dict{Symbol, Any}(:q => rf.q,
                                 :s_sep_mu_b => s_sep_mu_b,
                                 :r_sep_mu_b => fidf[cond, :r_fi_mu_b][1])

    s_diff_fun(zfun, mu_b) = fd[Symbol(:s_, zfun)](mu_b) - fidf[cond, Symbol(:s_fi_, zfun)][1]
    s_perc_diff_fun(zfun, mu_b) = (s_diff_fun(zfun, mu_b)/fidf[cond, Symbol(:s_fi_, zfun)][1]) .* 100.
    for zvar in zfuns
      # Safe Firm
      sepd[Symbol(:s_sep_, zvar)] = fd[Symbol(:s_, zvar)](s_sep_mu_b)

      # Differences and Percentage Differences
      col_prefix = Symbol(:s_sep_fi_, zvar)
      sepd[Symbol(col_prefix, :_diff)] = s_diff_fun(zvar, s_sep_mu_b)
      sepd[Symbol(col_prefix, :_perc_diff)] = s_perc_diff_fun(zvar, s_sep_mu_b)

      # Risky Firm
      sepd[Symbol(:r_sep_, zvar)] = fidf[cond, Symbol(:r_fi_, zvar)][1]
    end
    sepd[:eq_type] = :sep
    
    return DataFrame(sepd)
end

function get_sep_res(sf, rf, fidf::DataFrame)
    sf = setq(sf, q=.0)
    sepdf = DataFrame()
    for rq in fidf[:, :rq]
        try 
            rf = setq(rf, q=rq)
            dt = q_sep_res(fidf, sf, rf, rf.q)
            sepdf = vcat(sepdf, dt)
        catch 
            return sepdf
        end
    end
    
    return sepdf
end
 

# Plot Functions {{{1

# Inputs {{{2
xylabels = Dict{Symbol, Array{String,1}}(:mu_b => ["\\mu_b", "%.2f"],
                                         :m => ["m", "%.2f"],
                                         :q => ["q", "%.2f"],
                                         :mu_s => ["\\mu_s", "%.2f"])


zlabels = Dict{Symbol, Array{String,1}}(:c => ["Coupon", "%.2f"],
                                        :p => ["Principal", "%.2f"],
                                        :vb => ["VB", "%.1f"],
                                        :debt => ["Debt", "%.1f"],
                                        :eq => ["Equity", "%.1f"],
                                        :firm_value => ["Debt + Equity", "%1d"],
                                        :fv => ["Debt + Equity", "%1d"],
                                        :leverage => ["Leverage", "%1d"],
                                        :mbr => ["Market-to-Book Ratio", "%1d"], 
                                        :yield_spd => ["Bond Spreads (b.p.)", "%1d"])

svm_plots_title_params_order = [:mu_b, :m]

iso_cmaps = Dict{Symbol, Any}(:fi => Seaborn.get_cmap("YlGnBu_r"),
                              :mp => Seaborn.palplot("Reds"),
                              :pool => "BuPu",
                              :sep => "RdPu")

rmp_cmaps = Dict{String, Any}("misrep" => "GnBu",
                              "pooling" => "BuGn")

iso_plt_inputs = Dict{Symbol,Any}(:seaborn_style => "white",
                                  :iso_levels => 20,
                                  :heat_levels => 25,
                                  :fig_aspect => .4,
                                  :iso_fontsize => 9.,
                                  :use_subgrid => true,
                                  :subgrid_rows => 3,
                                  :iso_cols => 6,
                                  :heat_cols => 4,
                                  :title_font_size => 14.5,
                                  :fig_dpi => 300,
                                  :tight_pad => 3.,
                                  :h_pad => .75,
                                  :w_pad => .75)

eq_cat_dict = Dict{Symbol, Array{Any,1}}(:fi => [4, "FI"],
                                         :sep => [3, "SEP"],
                                         :pool => [2, "POOL"],
                                         :otc => [1, "OTC"])

contour_2pm_tlabels = Dict{Symbol, Array{String,1}}(:mu_s => ["\\mu_s", "%.2f"],
                                                :q => ["q", "%.2f"],
                                                :m => ["m", "%.2f"],
                                                :rbdisc => ["r^{b, EP}_{disc}", "%.2f"],
                                                :rbdisc_otc => ["r^{b, OTC}_{disc}", "%.2f"],
                                                :sig => ["\\sigma", "%.3f"],
                                                :D => ["D", "%.2f"],
                                                :V0 => ["V_0", "%.2f"],
                                                :debt => ["Debt", "%.1f"],
                                                :eq => ["Equity", "%.1f"],
                                                :fv => ["Firm Value", "%1d"],
                                                :leverage => ["Leverage", "%1d"],
                                                :mbr => ["Market-to-Book Ratio", "%1d"],
                                                :yield => ["Bond Yield", "%.2f"],
                                                :yield_spd => ["Bond Spread (b.p.)", "%.f"])

contour_2pm_plots_title_params_order = [:m, :D, :sigmal]

eq_type_2pm_title = Dict{Symbol, Array{Any,1}}(:fi => [:fi, "Full Information"],
                                               :mp => [:mp, "Misrepresentation"],
                                               :pool => [:pool, "Pooling"],
                                               :sep => [:sep, "Separating"],
                                               :ep => [:ep, "Prevailing EP Market Equilibria"],
                                               :dual => [:dual, "Prevailing Dual Market Equilibria"])

# Auxiliary Functions {{{2
function str_format_fun(a::String,b::Float64)
    return @eval @sprintf($a, $b)
end

function extract_val(df::DataFrame, xvar::Symbol, x::Float64, 
                     yvar::Symbol, y::Float64, zvar::Symbol; tol=1e-5)
    cond = true
    if xvar in Symbol.(names(df))
        cond = .&(cond, abs.(df[:, xvar] .- x) .< tol)
    end
    if yvar in Symbol.(names(df))
       cond = .&(cond, abs.(df[:, yvar] .- y) .< tol)
    end
    
    return df[cond, zvar][1]
end


function form_2pm_mesh_grid1(xvals::Array{Float64,1},
                             yvals::Array{Float64,1},
                             zfun; N::Int64=200)
    xgrid = range(minimum(xvals), stop=maximum(xvals), length=N)
    ygrid = range(minimum(yvals), stop=maximum(yvals), length=N)

    X = Array(repeat(xgrid, 1, N)')
    Y = repeat(ygrid, 1, N)
    Z = Array([zfun(x,y) for x in xgrid, y in ygrid]')

    return X, Y, Z
end


function form_2pm_mesh_grid(df::DataFrame, xvar::Symbol, 
                            yvar::Symbol, zvar::Symbol;
                            xvals::Array{Float64,1}=Array{Float64,1}(),
                            yvals::Array{Float64,1}=Array{Float64,1}())

    cols = [zvar]
    if xvar in Symbol.(names(df))
        cols = vcat(cols, xvar)
    end
    if yvar in Symbol.(names(df))
        cols = vcat(cols, yvar)
    end
    
    tmp = df[:, cols]
    
    if isempty(xvals)
        xvals = unique(tmp[:, xvar])
    end
    if isempty(yvals)
        yvals = unique(tmp[:, yvar])
    end

    X = Array(repeat(xvals, 1, size(yvals, 1))')
    Y = Array(repeat(yvals, 1, size(xvals, 1)))
    Z = Array([extract_val(tmp, xvar, x, yvar, y, zvar) 
                for x in xvals, y in yvals]')
    
    return X, Y, Z
end


function get_eq_contour_mesh_grid(xvals::Array{Float64,1}, yvals::Array{Float64,1},
                                  fun_dict; N::Int64=10^3)

    # eq_bool = (x, y) -> fun_dict[:fi_ind](x, y) + 2 * fun_dict[:sep_ind](x, y) + 3 * fun_dict[:pool_ind](x, y)
    # eq_vals = (x, y) -> (fun_dict[:fi_ind](x,y) * fun_dict[:mbr][:fi](x, y) +
    #                      fun_dict[:sep_ind](x,y) * fun_dict[:mbr][:sep](x, y) +
    #                      fun_dict[:pool_ind](x,y) * fun_dict[:mbr][:pool](x, y))

    X, Y, bool_Z = form_2pm_mesh_grid1(xvals, yvals, fun_dict[:eq_bool], N=N)
    _, _, bool_OTC_EP = form_2pm_mesh_grid1(xvals, yvals, fun_dict[:bool_otc_ep], N=N)
    _, _, r_mbr = fetch(@spawn form_2pm_mesh_grid1(xvals, yvals, fun_dict[:r_mbr], N=N))
    _, _, s_fv = fetch(@spawn form_2pm_mesh_grid1(xvals, yvals, fun_dict[:s_fv], N=N))

    return Dict{Symbol, Any}(:X => X, :Y => Y,
                             :bool_Z => bool_Z,
                             :bool_OTC_EP => bool_OTC_EP,
                             :r_mbr => r_mbr,
                             :s_fv => s_fv)
end


function form_eq_dict(df::DataFrame;
                      xvar::Symbol=:q, 
                      xvals::Array{Float64, 1}=Array{Float64, 1}(),
                      yvar::Symbol=:mu_s,
                      yvals::Array{Float64, 1}=Array{Float64, 1}(),
                      kx::Int64=3, ky::Int64=3)

    eq_type = df[:, :eq_type][1]
    xcol = xvar
    xcol = .&(eq_type == :fi, xvar == :q) ? :rq : xcol
    if isempty(xvals)
        xvals = unique(df[:, xcol])
    end
    
    if isempty(yvals)
        yvals = unique(df[:, yvar])
    end
    
    X, Y, Z = form_2pm_mesh_grid(df, xcol, yvar, 
                                 Symbol(:s_, eq_type, :_bpr); 
                                 xvals=xvals, yvals=yvals)
        
    sd = Dict{Symbol, Any}(:df => df,
                           :eq_type => eq_type,
                           :xvar => xvar,
                           :xvals => xvals,
                           :X => X,
                           :yvar => yvar,
                           :yvals => yvals,
                           :Y => Y,
                           :interp => Dict{Symbol, Any}())
    

    zcond(x::String) = any([occursin(String(y), x) for y in p2m.zfuns])
    zvars = [Symbol(x) for x in names(df) if .&(zcond(x), Symbol(x) != :eq_type)]
    for zvar in zvars
        _, _, Z = form_2pm_mesh_grid(df, xcol, yvar, zvar;
                                     xvals=xvals, yvals=yvals)
        sd[zvar] = Dict(:zvals => Z')
        
        sd[:interp][zvar] = Dierckx.Spline2D(sd[:xvals], sd[:yvals],
                                             sd[zvar][:zvals]; kx=kx, ky=ky)
    end
    
    return sd
end

function get_2pm_contour_plot_path_name(zvar::Symbol;
                                        plot_folder::String="/home/artur/BondPricing/Julia/Plots/M2P",
                                        contour_fname_ext::String="eps")
    
    return string(plot_folder, "/", zvar, ".", contour_fname_ext)
end

# Payoff Functions {{{2
function get_contour_equilibria_funs(fi_funs, mp_funs, pool_funs, sep_funs,
                                     fi_fv::Float64, fi_fv_fun, ilc_otc::Float64)
    # k_otc::Float64)

    jeq_ind = (x, y) -> mp_funs[:r_mp_mbr](x, y) >= fi_funs[:r_fi_mbr](x, y)
    fi_ind = (x, y) -> jeq_ind(x,y) == false
    pool_ind = (x, y) -> jeq_ind(x,y) ? pool_funs[:s_pool_fv](x, y) >= sep_funs[:s_sep_fv](x, y) : false
    sep_ind = (x, y) -> jeq_ind(x,y) ? !pool_ind(x, y) : false

    fun_dict = Dict{Symbol, Any}(:jeq_ind => jeq_ind,
                                 :fi_ind => fi_ind,
                                 :sep_ind => sep_ind,
                                 :pool_ind => pool_ind,
                                 :mbr => Dict{Symbol, Any}())
    for zvar in [:fv, :mbr, :lev]
        fun_dict[:mbr][:fi] = (x, y) -> fi_ind(x, y) ? fi_funs[Symbol(:r_fi_, zvar)](x, y) : .0
        fun_dict[:mbr][:pool] = (x, y) -> pool_ind(x, y) ? pool_funs[Symbol(:r_pool_, zvar)](x, y) : .0
        fun_dict[:mbr][:sep] = (x, y) -> sep_ind(x, y) ? sep_funs[Symbol(:r_sep_, zvar)](x, y) : .0
    end

    fun_dict[:eq_bool] = (x, y) -> (eq_cat_dict[:fi][1] * fun_dict[:fi_ind](x, y) +
                                    eq_cat_dict[:sep][1] * fun_dict[:sep_ind](x, y) +
                                    eq_cat_dict[:pool][1] * fun_dict[:pool_ind](x, y))
    fun_dict[:r_mbr] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_funs[:r_fi_mbr](x, y) +
                                  fun_dict[:sep_ind](x,y) * sep_funs[:r_sep_mbr](x, y) +
                                  fun_dict[:pool_ind](x,y) * pool_funs[:r_pool_mbr](x, y))
    fun_dict[:s_fv] = (x, y) -> (fun_dict[:fi_ind](x,y) * fi_fv + #fi_funs[:fv](x, y) +
                                 fun_dict[:sep_ind](x,y) * sep_funs[:s_sep_fv](x, y) +
                                 fun_dict[:pool_ind](x,y) * pool_funs[:s_pool_fv](x, y))
     fun_dict[:bool_otc] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(ilc_otc)) ? 1 : 0
     fun_dict[:bool_otc_ep] = (x, y) -> (fun_dict[:s_fv](x, y) < fi_fv_fun(ilc_otc)) ? 1 : fun_dict[:eq_bool](x, y)

     catd = Dict(zip([eq_cat_dict[x][1] for x in keys(eq_cat_dict)],
                     [eq_cat_dict[x][2] for x in keys(eq_cat_dict)]))
     fun_dict[:cat_otc_ep] = (x, y) ->  catd[fun_dict[:bool_otc_ep](x, y)]

    return fun_dict
end

 
# Full Information v.s. Misrepresentation {{{2
function plot_fi_misrep(df::DataFrame; 
                        figaspect::Float64=.8)
    
    #= q_vals = vcat(.0, df[:, :q]) =#
    #= fi_fv_vals = vcat(df[1, :s_fi_fv], df[:, :r_fi_fv]) =#
    #= misrep_fv_vals = vcat(df[1, :s_fi_fv], df[:, :r_misrep_fv]) =#
    #= fi_mbr_vals = vcat(df[1, :s_fi_mbr], df[:, :r_fi_mbr]) =#
    #= misrep_mbr_vals = vcat(df[1, :s_fi_mbr], df[:, :r_misrep_mbr]) =#

    q_vals = df[:, :q]
    fi_fv_vals = df[:, :r_fi_fv]
    misrep_fv_vals = df[:, :r_mp_fv]
    fi_mbr_vals = df[:, :r_fi_mbr]
    misrep_mbr_vals = df[:, :r_mp_mbr]

    Seaborn.set(style="darkgrid")
    fig = PyPlot.figure(figsize=(8,6))#figsize=Tuple(PyPlot.figaspect(my_figaspect)))
    ax1 = fig.add_subplot(211)
    ax1.plot(q_vals, fi_fv_vals, "--")
    ax1.plot(q_vals, misrep_fv_vals, "r-")
    ax1.set_title(" ")

    # Axes Limits ##################################################
    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()
    ax1.set_xlim([minimum(q_vals), maximum(q_vals)])
    ax1.set_ylim([ymin, ymax])
    ax1.set_ylabel("Firm Value", labelpad=10)
    # ##############################################################

    ax2 = fig.add_subplot(212)
    ax2.plot(q_vals, fi_mbr_vals, "--")
    ax2.plot(q_vals, misrep_mbr_vals, "r-")
    ax2.set_ylabel("Market-to-Book-Ratio", labelpad=10)
    # Axes Limits ##################################################
    xmin, xmax = ax2.get_xlim()
    ymin, ymax = ax2.get_ylim()
    ax2.set_xlim([minimum(q_vals), maximum(q_vals)])
    ax2.set_ylim([ymin, ymax])
    ax2.set_xlabel("Shock Probability (\$q\$)", labelpad=10)
    # ##############################################################

    fig.suptitle("Full Information v.s. Misrepresentation Payoffs")
    fig.tight_layout(pad=1.8)
    
    return fig
end

# Contour Plots {{{2
function get_2pm_contour_plot_title(fr, eq_type::Symbol,
                                    zvar::Symbol;
                                    zdf::Symbol=Symbol(""),
                                    rbdisc_otc::Float64=NaN,
                                    ft::Symbol=Symbol(""),
                                    mu_s::Float64=NaN,
                                    params_list::Array{Symbol,1}=[:V0, :D, :sig],
                                    s_fi_fv::Float64=NaN, s_otc_fv::Float64=NaN)

    title_eq_type = [eq_type_2pm_title[k][2] for k in keys(eq_type_2pm_title) 
                     if (eq_type_2pm_title[k][1] == eq_type)][1]

    if !isnan(mu_s)
        params_list = vcat(:mu_s, params_list)
    end

    title_params = join([[string("\$", contour_2pm_tlabels[x][1], "= \$ ",
                                 p2m.str_format_fun(contour_2pm_tlabels[x][2],
                                                    getfield(fr, x)))
                          for x in params_list]..., 
                         string("\$", contour_2pm_tlabels[:rbdisc][1], "= \$ ",
                                str_format_fun(contour_2pm_tlabels[:rbdisc][2], fr.rt.rbdisc))],  ", ")


    firm_type_title = ""
    if ft == :st
        firm_type_title = "Safe Type's "
    else
        firm_type_title = "Risky Type's "
    end

    if occursin("diff", String(zdf))
        differential = occursin("perc", String(zdf)) ? " Differential (%) " : " Differential "        
        plot_title = latexstring(firm_type_title, 
                                 eq_type_2pm_title[eq_type][2],
                                 " v.s. Full Information Eq. ",
                                 contour_2pm_tlabels[zvar][1], differential,
                                 "\n for ", title_params)
    else
        if eq_type == :mp
            plot_title = latexstring(firm_type_title,
                                     contour_2pm_tlabels[zvar][1],
                                     " in case of Misrepresentation",
                                     "\n for ", title_params)
        elseif eq_type in [:pool, :sep]
            plot_title = latexstring(firm_type_title, "Optimal ",
                                     contour_2pm_tlabels[zvar][1], " in a ",
                                     eq_type_2pm_title[eq_type][2], " Equilibrium ",
                                     "\n for ", title_params)
        else
            otc_title = ""
            if !isnan(rbdisc_otc)
                otc_title = string(", \$", contour_2pm_tlabels[:rbdisc_otc][1], "= \$ ",
                                     str_format_fun(contour_2pm_tlabels[:rbdisc_otc][2], rbdisc_otc))
            end

            plot_title = latexstring(firm_type_title, "Optimal ",
                                     contour_2pm_tlabels[zvar][1], " in the ",
                                     title_eq_type,
                                     "\n for ", title_params, otc_title)
        end
    end


    if .&(!isnan(s_fi_fv), !isnan(s_otc_fv))
        plot_title = latexstring(plot_title,
                                 "\n Safe Type's Full Information Firm Value = ",
                                 str_format_fun("%.2f", s_fi_fv),
                                 ", Safe Type's OTC Firm Value = ",
                                 str_format_fun("%.2f", s_otc_fv), ") \n")

    elseif !isnan(s_fi_fv)
        plot_title = latexstring(plot_title,
                                 "\n (Safe Type's Full Information Firm Value = ",
                                 str_format_fun("%.2f", s_fi_fv), ") \n")
    elseif !isnan(s_otc_fv)
        plot_title = latexstring(plot_title,
                                 "\n (Safe Type's OTC Firm Value = ",
                                 str_format_fun("%.2f", s_otc_fv), ") \n")
    end    

    return plot_title
end


function plot_2pm_iso_curves(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2};
                             iso_xlabel::Symbol=:q,
                             iso_ylabel::Symbol=:mu_s,
                             seaborn_style=iso_plt_inputs[:seaborn_style],
                             iso_levels=iso_plt_inputs[:iso_levels],
                             heat_levels=iso_plt_inputs[:heat_levels],
                             iso_cmap=iso_cmaps["full_info"],
                             heat_cmap::String="",
                             fig_aspect=iso_plt_inputs[:fig_aspect],
                             iso_fontsize=iso_plt_inputs[:iso_fontsize],
                             use_subgrid=iso_plt_inputs[:use_subgrid],
                             subgrid_rows=iso_plt_inputs[:subgrid_rows],
                             iso_cols=iso_plt_inputs[:iso_cols],
                             heat_cols=iso_plt_inputs[:heat_cols],
                             cat_Z=[],
                             cat_cmap="GnBu",
                             cat_alpha::Float64=.25,
                             rmp_cmap=rmp_cmaps["misrep"],
                             rmp_alpha::Float64=.15)

    if isempty(heat_cmap)
        heat_cmap = iso_cmap
    end

    if !isempty(seaborn_style)
        Seaborn.set(style=seaborn_style)
    end

    w, h = figaspect(fig_aspect)
    fig = PyPlot.figure(figsize=(w, h))

    # Choose between subgrids or subplots ##################################
    if use_subgrid
        fig = PyPlot.figure(figsize=(w, h))
        ax1 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols),
                                  (0, 0), rowspan=subgrid_rows, colspan=iso_cols)
        ax2 = PyPlot.subplot2grid((subgrid_rows, iso_cols + heat_cols),
                                  (0, iso_cols), rowspan=subgrid_rows, colspan=heat_cols)
    else
        fig, axs = PyPlot.subplots(1, 2, figsize=(w, h), sharey=true)
        ax1 = axs[1] # fig.add_subplot(121)
        ax2 = axs[2] # fig.add_subplot(122)
    end
    # ######################################################################

    CS = ax1.contour(X, Y, Z, levels=iso_levels, cmap=iso_cmap)
    ax1.clabel(CS, inline=5, fontsize=iso_fontsize)
    ax1.set_xlabel(latexstring("\$", xylabels[iso_xlabel][1], "\$"), labelpad=10)
    ax1.set_ylabel(latexstring("\$", xylabels[iso_ylabel][1], "\$"), labelpad=10)

    CS2 = ax2.contourf(X, Y, Z, levels=heat_levels, cmap=heat_cmap)
    if use_subgrid
        ax2.tick_params(
            axis="y",          # changes apply to the x-axis
            which="both",      # both major and minor ticks are affected
            bottom=false,      # ticks along the bottom edge are off
            top=false,         # ticks along the top edge are off
            left=false,
            right=false,
            labelleft=false,
            labelbottom=false)
    end
    ax2.set_xlabel(latexstring("\$", xylabels[iso_xlabel][1], "\$"), labelpad=10)

    # Add Colorbar
    cb2 = fig.colorbar(CS2)

    if !isempty(cat_Z)
        cats = sort(unique(cat_Z))
        cat_tick_labels = [eq_cat_dict[x][2] for x in [:fi, :sep, :pool, :otc]
                           if eq_cat_dict[x][1] in cats]

        if size(cats, 1) < size([x for x in keys(eq_cat_dict)], 1)
            cat_Z = cat_Z .- 1
            cats = cats .- 1
        end

        cat_levels = size(cats, 1) - 1
        CS1 = ax1.contourf(X, Y, cat_Z,
                           cmap=cat_cmap, levels=cat_levels) #, alpha=cat_alpha)
        cb1 = fig.colorbar(CS1, ax=ax1, ticks=reverse(cats))#, orientation="horizontal")
        cb1.set_ticklabels(cat_tick_labels)
        cb1.set_clim(1, cat_levels + 1)
    end
    
    return fig, ax1, ax2
end


function plot_2pm_iso_contour_curves(X, Y, Z;
                                     rmp_diff_fun=nothing,
                                     fig_title::LaTeXString=LaTeXString(""),
                                     file_path_name::String="",
                                     seaborn_style=iso_plt_inputs[:seaborn_style],
                                     iso_levels=iso_plt_inputs[:iso_levels],
                                     heat_levels=iso_plt_inputs[:heat_levels],
                                     iso_cmap=iso_cmaps["full_info"],
                                     heat_cmap::String="",
                                     fig_aspect=iso_plt_inputs[:fig_aspect],
                                     iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                     use_subgrid=iso_plt_inputs[:use_subgrid],
                                     subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                     iso_cols=iso_plt_inputs[:iso_cols],
                                     heat_cols=iso_plt_inputs[:heat_cols],
                                     title_font_size=iso_plt_inputs[:title_font_size],
                                     fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                     tight_pad=iso_plt_inputs[:tight_pad],
                                     h_pad=iso_plt_inputs[:h_pad],
                                     w_pad=iso_plt_inputs[:w_pad],
                                     rmp_cmap=rmp_cmaps["misrep"],
                                     rmp_alpha::Float64=.15)

    fig, ax1, ax2 = plot_2pm_iso_curves(X, Y, Z;
                                        seaborn_style=seaborn_style,
                                        iso_levels=iso_levels,
                                        heat_levels=heat_levels,
                                        iso_cmap=iso_cmap,
                                        heat_cmap=heat_cmap,
                                        fig_aspect=fig_aspect,
                                        iso_fontsize=iso_fontsize,
                                        use_subgrid=use_subgrid,
                                        subgrid_rows=subgrid_rows,
                                        iso_cols=iso_cols,
                                        heat_cols=heat_cols,
                                        rmp_cmap=rmp_cmap,
                                        rmp_alpha=rmp_alpha)
    
    if !isempty(fig_title)
        fig.suptitle(fig_title, fontsize=title_font_size)
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig
end


# JEQ Contour Plots {{{2

function plot_2pm_equilibria_iso_contour_curves(X, Y, Z, eq_type_Z;
                                                ptd::Dict{Symbol,Dict{Symbol,Any}}=Dict{Symbol,Dict{Symbol,Any}}(),
                                                sup_title_sep::Bool=false,
                                                fig_title::LaTeXString=LaTeXString(""),
                                                file_path_name::String="",
                                                iso_xlabel::Symbol=:q,
                                                iso_ylabel::Symbol=:mu_s,
                                                seaborn_style=iso_plt_inputs[:seaborn_style],
                                                iso_levels=iso_plt_inputs[:iso_levels],
                                                heat_levels=iso_plt_inputs[:heat_levels],
                                                iso_cmap=iso_cmaps[:fi],
                                                heat_cmap::String="",
                                                fig_aspect=iso_plt_inputs[:fig_aspect],
                                                iso_fontsize=iso_plt_inputs[:iso_fontsize],
                                                use_subgrid=iso_plt_inputs[:use_subgrid],
                                                subgrid_rows=iso_plt_inputs[:subgrid_rows],
                                                iso_cols=iso_plt_inputs[:iso_cols],
                                                heat_cols=iso_plt_inputs[:heat_cols],
                                                title_font_size=iso_plt_inputs[:title_font_size],
                                                fig_dpi::Int64=iso_plt_inputs[:fig_dpi],
                                                tight_pad=iso_plt_inputs[:tight_pad],
                                                h_pad=iso_plt_inputs[:h_pad],
                                                w_pad=iso_plt_inputs[:w_pad],
                                                cat_cmap="GnBu",
                                                rmp_cmap=rmp_cmaps["misrep"])

    fig, ax1, ax2 = plot_2pm_iso_curves(X, Y, Z;
                                        iso_xlabel=iso_xlabel,
                                        iso_ylabel=iso_ylabel,
                                        iso_cmap=iso_cmap,
                                        iso_levels=15,
                                        cat_Z=eq_type_Z,
                                        cat_cmap=cat_cmap,
                                        rmp_cmap=rmp_cmap)
#    ax1.contourf(X, Y, eq_type_Z, cmap="GnBu_r", alpha=.4)

    if !isempty(ptd)
        dq = 25e-3 * ptd[:p1][:q]
        dmu_s = 25e-3 * ptd[:p1][:mu_s]
        for pt in keys(ptd)
            ax1.scatter(ptd[pt][:q], ptd[pt][:mu_s], 
                        marker=ptd[pt][:marker], 
                        c=ptd[pt][:color],
                        label=string(pt))
            ax1.annotate(string(pt, "  = (", ptd[pt][:q], ", ", ptd[pt][:mu_s], ")"), 
                         (ptd[pt][:q] + dq, ptd[pt][:mu_s] + dmu_s))
            # ax2.scatter(ptd[pt][:q], ptd[pt][:mu_s],
            #             marker=ptd[pt][:marker], 
            #             c=ptd[pt][:color])
        end
    end

    if !isempty(fig_title)
        ttl = fig.suptitle(fig_title, fontsize=title_font_size)

        if sup_title_sep
            ttl.set_position([.5, 1.05])
        end
    end
    PyPlot.tight_layout(pad=tight_pad, h_pad=h_pad, w_pad=w_pad)

    if !isempty(file_path_name)
        PyPlot.savefig(file_path_name, dpi=fig_dpi, bbox_inches="tight")
    end

    return fig, ax1, ax2
end


# Numeric Functions {{{1
function expected_payoff(fr, x0, f; xf::Float64=10., xN::Int64=10^3) 
    # muxh = x0 + .5 * fr.sig^2 
    muxh = fr.rt.rf - .5 * fr.sig^2 

    du = Distributions.Normal(muxh, fr.sig)
    x1grid = range(x0 - xf  * fr.sig, stop=muxh + xf * fr.sig, length=xN)
    d1x = x1grid[2] - x1grid[1]
    yu = sum([f(x1) * pdf(du, x1) for x1 in x1grid]) * d1x
    
    yl = 0.
    if fr.q > 0
        muxl = muxh - 1. * fr.sig
        dl = Distributions.Normal(muxl, fr.sig)
        x1gridl = x1grid .- (muxh - muxl)
        
        yl = sum([f(x1) * pdf(dl, x1) for x1 in x1gridl]) * d1x
    end
    
    return (1 - fr.q) * yu + fr.q * yl # sum([f(x1) * pdf(d, x1) for x1 in x1grid]) * d1x
end

# Asset Value
AV(fr, x::Float64) = fr.V0 * exp(x)

# Tax Shields
TS(fr) = fr.pi * fr.mu_b * fr.D

# Bankrupcty Condition
BKR(fr, x::Float64) = AV(fr, x) + TS(fr) < fr.mu_b*fr.D

function num_bondpr(fr, x1; xf::Float64=10., xN::Int64=10^3) 
    bpr = .0
    if fr.mu_b > 0
        dp(x2) = BKR(fr, x2) ? fr.alpha * AV(fr, x2) : fr.mu_b*fr.D 
        bpr = exp(-fr.rt.rbdisc) * expected_payoff(fr, x1, dp; xf=xf, xN=xN)/fr.mu_b
    end
    
    return bpr
end

function num_debtpr(fr, x1::Float64; xf::Float64=10., xN::Int64=10^3) 
    return num_bondpr(fr, x1::Float64; xf=xf, xN=xN) * fr.mu_b
end

function num_eqpr(fr, x1::Float64; xf::Float64=10., xN::Int64=10^3) 
    ep(x) = max(AV(fr, x) + TS(fr) - fr.mu_b * fr.D, .0) 
    return exp(-fr.rt.rf) * expected_payoff(fr, x1, ep; xf=xf, xN=xN)
end

function num_mbr(fr, x1::Float64; xf::Float64=10., xN::Int64=10^3) 
    return num_eqpr(fr, x1; xf=xf, xN=xN) / (fr.V0 - num_debtpr(fr, x1; xf=xf, xN=xN))
end

# Numeric Get Methods {{{2
function get_num_dpr(fr; x1::Float64=0., q::Float64=NaN, 
                  mu_b::Float64=NaN, xf::Float64=10., xN::Int64=10^3) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return num_debtpr(fr, x1; xf=xf, xN=xN)
end

function get_num_epr(fr; x1::Float64=0., q::Float64=NaN, 
                  mu_b::Float64=NaN, xf::Float64=10., xN::Int64=10^3) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return num_eqpr(fr, x1; xf=xf, xN=xN)
end

function get_num_fv(fr; x1::Float64=0., q::Float64=NaN, 
                 mu_b::Float64=NaN, xf::Float64=10., xN::Int64=10^3) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return (num_debtpr(fr, x1; xf=xf, xN=xN) 
            + num_eqpr(fr, x1; xf=xf, xN=xN))
end

function get_num_mbr(fr; x1::Float64=0., q::Float64=NaN, 
                  mu_b::Float64=NaN, xf::Float64=10., xN::Int64=10^3) 
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    return num_mbr(fr, x1; xf=xf, xN=xN)
end

function get_num_misrep_mbr(sf, rf; x1::Float64=0., 
                            sq::Float64=NaN, rq::Float64=NaN, 
                            mu_b::Float64=NaN, xf::Float64=10., xN::Int64=10^3) 
    # Safe Firm
    sf = setq(sf; q=sq)
    sf = set_mu_b(sf; mu_b=mu_b)
    
    # Risky Firm
    rf = setq(rf; q=rq)
    rf = set_mu_b(rf; mu_b=mu_b)
    
    
    rE1 = get_num_epr(rf)
    sD1 = get_num_dpr(sf)
    
    return  rE1 / (rf.V0 - sD1)
end

function get_num_prices(sf, rf; 
                        mu_b_min::Float64=.0, 
                        mu_b_max::Float64=NaN, mubN::Int64=30)
    mu_b_max = isnan(mu_b_max) ? max(sf.V0/sf.D, rf.V0/rf.D) : mu_b_max
    mu_b_grid = range(mu_b_min, stop=mu_b_max, length=mubN)
    
    sdpr_vec = [get_num_dpr(sf; mu_b=mu_b) for mu_b in mu_b_grid]
    sepr_vec = [get_num_epr(sf; mu_b=mu_b) for mu_b in mu_b_grid]
    rdpr_vec = [get_num_dpr(rf; mu_b=mu_b) for mu_b in mu_b_grid]
    repr_vec = [get_num_epr(rf; mu_b=mu_b) for mu_b in mu_b_grid]
    sfv_vec = [get_num_fv(sf; mu_b=mu_b) for mu_b in mu_b_grid]
    rfv_vec = [get_num_fv(rf; mu_b=mu_b) for mu_b in mu_b_grid]
    smbr_vec = [get_num_mbr(sf; mu_b=mu_b) for mu_b in mu_b_grid]
    rmbr_vec = [get_num_mbr(rf; mu_b=mu_b) for mu_b in mu_b_grid]
    r_misrep_mbr_vec = [get_num_misrep_mbr(sf, rf; mu_b=mu_b) for mu_b in mu_b_grid]

    df = DataFrame(:qs => sf.q, :qr => rf.q,
                   :mu_b => mu_b_grid, 
                   :s_debt => sdpr_vec, :s_equity => sepr_vec, 
                   :s_fv => sfv_vec, :s_mbr => smbr_vec,
                   :r_debt => rdpr_vec, :r_equity => repr_vec, 
                   :r_fv => rfv_vec, :r_mbr => rmbr_vec, 
                   :r_misrep_mbr => r_misrep_mbr_vec)
    
    return df
end

# 3-Period Model Functions {{{1
# Debt Rollover Cost
num_drc(fr, x1::Float64) = fr.mu_b*(num_bondpr(fr, x1) - fr.D)

# Equity Dilution
num_eq_diff(fr, x::Float64, x1::Float64) = FV(fr, x) + num_drc(fr, x1)

function get_x1(fr, x::Float64; 
                q::Float64=NaN, mu_b::Float64=NaN, 
                xf::Float64=3., xN::Int64=10^5) 
    
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    # DELTA x
    diff(x1::Float64) = num_eq_diff(fr, x, x1) - (1 + fr.rt.rf) * fr.V0 * exp(x1)
    
    x1grid = range(x - xf  * fr.sig, stop=x + xf * fr.sig, length=xN)
    
    x1diff = [diff(x1) for x1 in x1grid]
    
    x1star = x1grid[argmin(abs.(x1diff))]
    
    println(diff(x1star))
    
    return x1star#, eqdiff(x, x1star)
end

function num_default_cond(fr, x1::Float64; q::Float64=NaN, mu_b::Float64=NaN)
    
    # NET CASH FLOW:
    ncf = fr.mu_b*(num_bondpr(fr, x1) - (1-fr.pi) * fr.D)
    
    # Determine Default
    E1 = num_eqpr(fr, x1)
    
    # Default Condition
    return .&(ncf < 0, abs.(ncf) > E1)
end

# Bond Payoff At time t=1        
function t1_payoffs(fr, x1::Float64; 
                    q::Float64=NaN, mu_b::Float64=NaN)
    
    fr = setq(fr; q=q)
    fr = set_mu_b(fr; mu_b=mu_b)
    
    default = false
    bp = .0
    ep = NaN
    x1star = NaN
    if fr.mu_b > 0.
#         default = default_cond(fr, x1star)
        
        x1star = get_x1(fr, x1)
        
        # NET CASH FLOW:
        ncf = fr.mu_b*(num_bondpr(fr, x1star) - (1-fr.pi) * fr.D)
        
        # Determine Default
        E1 = num_eqpr(fr, x1star)
        
        # Default Condition
        default = .&(ncf < 0, abs.(ncf) > E1)
 
        bp = default ? fr.alpha * FV(fr, x1) / fr.mu_b : fr.D
        ep = default ? 0. : E1
    else
        ep = num_eqpr(fr, x1)
    end
        
    return Dict{Symbol, Any}(:default => default,
                             :x1 => x1,
                             :x1star => x1star,
                             :bond => bp, :equity => ep)            
end    

# function t0_prices(fr;
#                    q::Float64=NaN, mu_b::Float64=NaN, 
#                    xf::Float64=3., xN::Int64=10^3)
#     x1grid = range(- xf  * fr.sig, stop= xf * fr.sig, length=xN)
    

# End Module {{{1
end

