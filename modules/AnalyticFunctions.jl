module_path = "/home/artur/BondPricing/bond-model/modules/"
push!(LOAD_PATH, module_path)
# modlpath = joinpath(module_path, "ModelObj")
# if !(modlpath in LOAD_PATH)
#     push!(LOAD_PATH, modlpath)
# end

module AnalyticFunctions

using Distributions
using Interpolations

# User-defined package:
using ModelObj: extract_param


# * CVM Functions ########################################
# include("cvm_funs.jl")
# Analytic Functions for the CVM Model

# Rates
function rgrow(r, gross_delta)
    return r - gross_delta
end

function rdisc(r, xi, k)
    return r + xi * k
end

# Auxiliary Functions
function cvm_a(sigma, r, gross_delta)
   return  (rgrow(r, gross_delta) - .5 * sigma ^ 2) / sigma ^ 2
end


function cvm_z(sigma, r, gross_delta)
    return (cvm_a(sigma, r, gross_delta)^2 * sigma^4 + 2 * r * sigma^2)^.5 / sigma^2
end


function cvm_zhat(sigma, r, gross_delta, xi, k)
    return (cvm_a(sigma, r, gross_delta)^2 * sigma^4 +
             2 * rdisc(r, xi, k) * sigma^2)^.5 / sigma^2
end


function cvm_eta(sigma, r, gross_delta)
    return cvm_z(sigma, r, gross_delta) - cvm_a(sigma, r, gross_delta)
end
  
  
# ######### Bond Price-Specific Auxiliary Functions #########
function cvm_h1(v, ttm, sigma, r, gross_delta)
    return (-v - cvm_a(sigma, r, gross_delta) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end


function cvm_h2(v, ttm, sigma, r, gross_delta)
    return (-v + cvm_a(sigma, r, gross_delta) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end


function cvm_q1(v, ttm, sigma, r, gross_delta, xi, k)
    return (-v - cvm_zhat(sigma, r, gross_delta, xi, k) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end


function cvm_q2(v, ttm, sigma, r, gross_delta, xi, k)
    return (-v + cvm_zhat(sigma, r, gross_delta, xi, k) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end
  

function cvm_F(v, ttm, sigma, r, gross_delta)
    return Distributions.cdf(Normal(), cvm_h1(v, ttm, sigma, r, gross_delta)) +
            exp(-2 * cvm_a(sigma, r, gross_delta) * v) *
           Distributions.cdf(Normal(), cvm_h2(v, ttm, sigma, r, gross_delta))
end


function cvm_G(v, ttm, sigma, r, gross_delta, xi, k)
    return exp((- cvm_a(sigma, r, gross_delta) + cvm_zhat(sigma, r, gross_delta, xi, k)) * v) *
             Distributions.cdf(Normal(), cvm_q1(v, ttm, sigma, r, gross_delta, xi, k)) +
           exp((- cvm_a(sigma, r, gross_delta) - cvm_zhat(sigma, r, gross_delta, xi, k)) * v) *
             Distributions.cdf(Normal(), cvm_q2(v, ttm, sigma, r, gross_delta, xi, k))
end
  
# Bond Price
function cvm_bond_price(vb, v, ttm, mu_b, m, c, p, sigma, r, gross_delta, xi, k, alpha, pi)
    return (c / rdisc(r, xi, k)) +
           exp(-rdisc(r, xi, k) * ttm) * (p - c / rdisc(r, xi, k)) * (1 - cvm_F(v, ttm, sigma, r, gross_delta)) +
           (alpha * vb / (mu_b * m) - c / rdisc(r, xi, k)) * cvm_G(v, ttm, sigma, r, gross_delta, xi, k)
end
  
  
# Credit-Risk-Free Bond Price
function rfbond_price(ttm, c, p, r, xi, k)
    return ((c / rdisc(r, xi, k)) +
            (p - c / rdisc(r, xi, k)) * exp(- rdisc(r, xi, k) * ttm))
end
  
  
######### vb-Specific Auxiliary Functions #########
function cvm_b(x, m, sigma, r, gross_delta, xi, k)
    return (1 / (cvm_z(sigma, r, gross_delta) + x)) *
          exp(- rdisc(r, xi, k) * m) *
          (Distributions.cdf(Normal(), x * sigma * sqrt(m)) -
           exp(r * m) * Distributions.cdf(Normal(), - cvm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
end


function cvm_B(x, m, sigma, r, gross_delta)
    return (1 / (cvm_z(sigma, r, gross_delta) + x)) *
            (Distributions.cdf(Normal(), x * sigma * sqrt(m)) -
              exp(.5 * (cvm_z(sigma, r, gross_delta)^2 - x^2) * sigma^2 * m) *
              Distributions.cdf(Normal(), - cvm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
end
  

function cvm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi)
    return ((1 - pi) * (c * m) +
             (1 - exp(- rdisc(r, xi, k) * m)) *
             (p - c / rdisc(r, xi, k))) / cvm_eta(sigma, r, gross_delta)
end


function cvm_numerator2(m, c, p, sigma, r, gross_delta, xi, k)
    return (p - c / rdisc(r, xi, k)) *
            (cvm_b(- cvm_a(sigma, r, gross_delta), m, sigma, r, gross_delta, xi, k) +
              cvm_b(cvm_a(sigma, r, gross_delta), m, sigma, r, gross_delta, xi, k)) +
             (c / rdisc(r, xi, k)) *
             (cvm_B(-cvm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta) +
              cvm_B(cvm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta))
end


function cvm_denominator(m, sigma, r, gross_delta, iota, xi, k, alpha)
    return (gross_delta - iota) / (cvm_eta(sigma, r, gross_delta) - 1) + (alpha / m) *
            (cvm_B(- cvm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta) +
             cvm_B(cvm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta))
end


# Default Boundary
function cvm_vb(mu_b, m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
    value = mu_b * (cvm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi) +
                    cvm_numerator2(m, c, p, sigma, r, gross_delta, xi, k)) /
                    cvm_denominator(m, sigma, r, gross_delta, iota, xi, k, alpha)
  
    return max(convert(Float64, value), 1e-4)
end


######### Equity-Specific Auxiliary Functions #########
function cvm_gamma(sigma, r, gross_delta)
    return cvm_a(sigma, r, gross_delta) + cvm_z(sigma, r, gross_delta)
end
  

function cvm_k_fun(v, x, w, u, sigma, m)
    return exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
           exp(-u * v) * Distributions.cdf(Normal(),(-v + (u - x) * (sigma^2) * m) / (sigma * sqrt(m))) -
           exp(-(x + w) * v) * Distributions.cdf(Normal(), (-v + w * (sigma^2) * m) / (sigma * sqrt(m)))
end


function cvm_K(v, x, w, u, sigma, m)
    return (Distributions.cdf(Normal(), w * sigma * sqrt(m)) -
              exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
              Distributions.cdf(Normal(),(u - x) * sigma * sqrt(m))) * exp(- u * v) +
             cvm_k_fun(v, x, w, u, sigma, m)
end

  
function cvm_A(v, y, sigma, r, gross_delta, m)
    return 1 / (cvm_z(sigma, r, gross_delta) - y) *
        (cvm_K(v, cvm_a(sigma, r, gross_delta), y, cvm_gamma(sigma, r, gross_delta), sigma, m) +
         cvm_k_fun(v, cvm_a(sigma, r, gross_delta), - y, - cvm_eta(sigma, r, gross_delta), sigma, m)) +
        1 / (cvm_z(sigma, r, gross_delta) + y) *
        (cvm_K(v, cvm_a(sigma, r, gross_delta), - y, cvm_gamma(sigma, r, gross_delta), sigma, m) +
         cvm_k_fun(v,cvm_a(sigma, r, gross_delta), y, - cvm_eta(sigma, r, gross_delta), sigma, m))
end


# Equity Value
function cvm_eq(v, mu_b, m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
    vbl = cvm_vb(mu_b, m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
    net_delta = (gross_delta - iota)

    cf1 = 1 / (r - rgrow(r, gross_delta))
    cf20 = 1 / (sigma^2 * cvm_z(sigma, r, gross_delta))
    cf21 = 1 / (cvm_gamma(sigma, r, gross_delta) + 1) 
    cf2 = - cf20 * cf21
    cf30 = (1 - pi) * (c * m)
    cf31 = (1 - exp(- rdisc(r, xi, k) * m))
    cf32 = (p - c / rdisc(r, xi, k))
    cf3 = - (cf30 +  cf31 * cf32) * cf20
    cf4 = 1 / cvm_eta(sigma, r, gross_delta)
    cf5 = 1 / cvm_gamma(sigma, r, gross_delta)
    cf6 = cf20
    cf7 = exp(- rdisc(r, xi, k) * m) * cf32
    cf8 = alpha/m
    cf9 = - c / rdisc(r, xi, k)
    
    vf1 = exp(v)
    vf2 = exp(-cvm_gamma(sigma, r, gross_delta) * v)
    vf3 = (1 - exp(- cvm_gamma(sigma, r, gross_delta) * v))
    vf4 = cvm_A(v, cvm_a(sigma, r, gross_delta), sigma, r, gross_delta, m) 
    vf5 = cvm_A(v, cvm_zhat(sigma, r, gross_delta, xi, k), sigma, r, gross_delta, m)
    
    return  mu_b * (net_delta * (vbl/mu_b) * (cf1 * vf1 +  cf2 * vf2) +
                    cf3 * (cf4 + cf5 * vf3) +
                    cf6 * (cf7 * vf4 - (cf8 * (vbl/mu_b) + cf9) * vf5))
end


# * SVM Analytic Functions ###############################
# include("svm_analytic_funs.jl")
# SVM Functions

# ######## Pre-Volatility Shock Bond Pricing Functions #########
# Pre-Volatility Shock Auxiliary Functions
function rdisc_pvs(r, xi, k, _lambda)
    return r + xi * k + _lambda
end

function zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)
    return sqrt(cvm_a(sigma, r, gross_delta)^2 * sigma^4 +
                2 * rdisc_pvs(r, xi, k, _lambda) * sigma^2) / sigma^2
end

function q1_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
    return (-v - zhat_pvs(sigma, r, gross_delta, xi, k, _lambda) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end

function q2_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
    return (-v + zhat_pvs(sigma, r, gross_delta, xi, k, _lambda) * sigma^2 * ttm) / (sigma * sqrt(ttm))
end

function G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
    return exp((- cvm_a(sigma, r, gross_delta) + zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)) * v) *
        Distributions.cdf(Normal(), q1_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)) +
        exp((- cvm_a(sigma, r, gross_delta) - zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)) * v) *
        Distributions.cdf(Normal(), q2_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda))
end

function psi_v_td(vt, v, ttm, sigma, r, gross_delta)
    return (Distributions.cdf(Normal(), (-v + vt) / (sigma * sqrt(ttm)) + cvm_a(sigma, r, gross_delta) * sigma * sqrt(ttm))
            - exp(-2 * cvm_a(sigma, r, gross_delta) * vt) *
            Distributions.cdf(Normal(), (-v - vt) / (sigma * sqrt(ttm)) + cvm_a(sigma, r, gross_delta) * sigma * sqrt(ttm)))
end

function dv_psi_v_td(vt, v, ttm, sigma, r, gross_delta)
    return (-1.0 ./(sigma * sqrt(ttm))) *
        (Distributions.pdf(Normal(), (-v + vt) / (sigma * sqrt(ttm)) + cvm_a(sigma, r, gross_delta) * sigma * sqrt(ttm))
         - exp(-2 * cvm_a(sigma, r, gross_delta) * vt) *
         Distributions.pdf(Normal(), (-v - vt) / (sigma * sqrt(ttm)) + cvm_a(sigma, r, gross_delta) * sigma * sqrt(ttm)))
end

# Market Value of (Credit-Risk-Free) Debt at Upper Barrier:
function rf_debt(mu_b, m, c, p, r, xi, k)
    return mu_b * (((c / rdisc(r, xi, k)) * (m - (1 - exp(-(rdisc(r, xi, k)) * m)) / (rdisc(r, xi, k))) +
                    p * (1 - exp(-(rdisc(r, xi, k)) * m)) / (rdisc(r, xi, k))))
end

# Notice Psi = 1- F and does not depend on lambda
# Prior to the volatility shock, G_pvs takes rdisc = r + xi*k + lambda!!!!
function default_vb_payoff(vb, mu_b, m, alpha)
    return (alpha * vb) / (mu_b * m)
end

function on_default_payoff(v, vb, ttm, mu_b, m, c, p, r, xi, k, alpha)
    dvbpayoff = default_vb_payoff(vb, mu_b, m, alpha) * exp(v)
    rfbond = rfbond_price(ttm, c, p, r, xi, k)
    rfdebt = rf_debt(mu_b, m, c, p, r, xi, k)

    return dvbpayoff * (rfdebt > dvbpayoff) + rfbond * (rfdebt <= dvbpayoff)
end

function coupon_pv_fun(v, ttm, c, sigma, r, gross_delta, xi, k, _lambda)
    return (c / rdisc_pvs(r, xi, k, _lambda)) * ((1 - exp(- rdisc_pvs(r, xi, k, _lambda) * ttm) *
           (1 - cvm_F(v, ttm, sigma, r, gross_delta))) -
           G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda))
end

function default_payoff_pv_fun(v, vb, ttm, mu_b, m, sigma, r, gross_delta, xi, k, alpha, _lambda)
    return default_vb_payoff(vb, mu_b, m, alpha) *
           G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
end

function maturity_payoff_pv_fun(v, ttm, p, sigma, r, gross_delta, xi, k, _lambda)
    return p * exp(- rdisc_pvs(r, xi, k, _lambda) * ttm) * (1 - cvm_F(v, ttm, sigma, r, gross_delta))
end

function no_vol_shock_cf_pv(v, vb, ttm, mu_b, m, c, p, sigma, r, gross_delta, xi, k, alpha, _lambda)
    return coupon_pv_fun(v, ttm, c, sigma, r, gross_delta, xi, k, _lambda) +
           default_payoff_pv_fun(v, vb, ttm, mu_b, m, sigma, r, gross_delta, xi, k, alpha, _lambda) +
           maturity_payoff_pv_fun(v, ttm, p, sigma, r, gross_delta, xi, k, _lambda)
end

# Notice that time-to-maturity becomes ttm-u and maturity becomes ttm:
function vol_shock_cf_integrand(vt, v_prior, v_post, vb_prior, vb_post, u, ttm, mu_b, c, p,
                                sigma_prior, sigma_post, r, gross_delta, xi, k, alpha, pi, _lambda)
    return _lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
           cvm_bond_price(vb_post, v_post, ttm - u, mu_b, ttm, c, p, sigma_post, r, gross_delta, xi, k, alpha, pi) *
           (-dv_psi_v_td(vt, v_prior, u, sigma_prior, r, gross_delta) / (vb_prior * exp(v_prior)))
end


# * Analytic Functions Get Methods #######################
# include("af_get_methods.jl")
function get_k_struct(svm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN, p::Float64=NaN)
    if isnan(mu_b)
        mu_b = svm.mu_b
    end

    if isnan(m)
        m=svm.m
    end
    
    if isnan(c)
        c=svm.c
    end

    if isnan(p)
        p=svm.p
    end
    
    return mu_b, m, c, p
end


function get_rgrow(svm)
   return rgrow(svm.pm.r, svm.pm.gross_delta) 
end


function get_rdisc(svm)
   return rdisc(svm.pm.r, svm.pm.xi, svm.pm.kappa) 
end


function get_cvm_vb(svm, sigma;
                    mu_b::Float64=NaN, m::Float64=NaN,
                    c::Float64=NaN, p::Float64=NaN)
    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)

    return cvm_vb(mu_b, m, c, p, 
                  sigma, svm.pm.r, svm.pm.gross_delta, 
                  svm.pm.iota, svm.pm.xi, svm.pm.kappa,
                  svm.pm.alpha, svm.pm.pi)
end


function get_agg_c(svm; mu_b::Float64=NaN, m::Float64=NaN, c::Float64=NaN)
    mu_b, m, c, _ = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=NaN)
   
    return mu_b * c * m
end


function get_agg_p(svm; mu_b::Float64=NaN, m::Float64=NaN, p::Float64=NaN)
    mu_b, m, _, p = get_k_struct(svm; mu_b=mu_b, m=m, c=NaN, p=p)
   
    return mu_b * p * m
end


function get_param(svm, pname::Symbol)
    val = NaN
    # svm_pos = findin([string(x) for x in fieldnames(svm)], [pname])
    # svm_pm_pos = findin([string(x) for x in fieldnames(svm.pm)], [pname])
    
    svm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm))])
    svm_pm_pos = findall(x -> x==pname, [x for x in fieldnames(typeof(svm.pm))])
    
    if pname == :C
        return get_agg_c(svm; mu_b=svm.mu_b, m=svm.m, c=svm.c)
    elseif pname == :P
        return get_agg_p(svm; mu_b=svm.mu_b, m=svm.m, p=svm.p)
    elseif pname == :delta
        return svm.pm.gross_delta - svm.pm.iota
    else
        return extract_param(svm, pname)
    end
    
    # elseif length(svm_pos) > 0
    #     # return getfield(svm, fieldnames(svm)[svm_pos[1]])
    #     return getfield(svm, fieldnames(typeof(svm))[svm_pos[1]])
    # else
    #     # return getfield(svm.pm, fieldnames(svm.pm)[svm_pm_pos[1]])
    #     return getfield(svm.pm, fieldnames(typeof(svm.pm))[svm_pm_pos[1]])
    # end
end


function get_leverage(debt::Union{Float64, Array{Float64, 1}},
                      equity::Union{Float64, Array{Float64, 1}})
    return (debt ./ (debt .+ equity)) .* 100
end

function get_mbr(V0::Float64,
                 debt::Union{Float64, Array{Float64, 1}},
                 equity::Union{Float64, Array{Float64, 1}})
    return (equity ./ (V0 .- debt) .- 1.) .* 100
end


# * END MODULE
end
