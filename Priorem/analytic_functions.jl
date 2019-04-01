# analytic functions
module AnalyticFunctions

    using Distributions
    using Interpolations
    using ProgressMeter

    function rgrow(r, gross_delta)
        r - gross_delta
    end

    function rdisc(r, xi, k)
        r + xi * k
    end

    # Auxiliary Functions
    function bsm_a(sigma, r, gross_delta)
        (rgrow(r, gross_delta) - .5 * sigma ^ 2) / sigma ^ 2
    end

    function bsm_z(sigma, r, gross_delta)
        return (bsm_a(sigma, r, gross_delta)^2 * sigma^4 + 2 * r * sigma^2)^.5 / sigma^2
    end

    function bsm_zhat(sigma, r, gross_delta, xi, k)
        (bsm_a(sigma, r, gross_delta)^2 * sigma^4 +
                          2 * rdisc(r, xi, k) * sigma^2)^.5 / sigma^2
    end

    function bsm_eta(sigma, r, gross_delta)
        bsm_z(sigma, r, gross_delta) - bsm_a(sigma, r, gross_delta)
    end


    # ######### Bond Price-Specific Auxiliary Functions #########
    function bsm_h1(v, ttm, sigma, r, gross_delta)
        return (-v - bsm_a(sigma, r, gross_delta) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function bsm_h2(v, ttm, sigma, r, gross_delta)
        return (-v + bsm_a(sigma, r, gross_delta) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function bsm_q1(v, ttm, sigma, r, gross_delta, xi, k)
        return (-v - bsm_zhat(sigma, r, gross_delta, xi, k) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function bsm_q2(v, ttm, sigma, r, gross_delta, xi, k)
        return (-v + bsm_zhat(sigma, r, gross_delta, xi, k) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function bsm_F(v, ttm, sigma, r, gross_delta)
        return cdf(Normal(), bsm_h1(v, ttm, sigma, r, gross_delta)) +
               exp(-2 * bsm_a(sigma, r, gross_delta) * v) *
                cdf(Normal(), bsm_h2(v, ttm, sigma, r, gross_delta))
    end

    function bsm_G(v, ttm, sigma, r, gross_delta, xi, k)
        return exp((- bsm_a(sigma, r, gross_delta) + bsm_zhat(sigma, r, gross_delta, xi, k)) * v) *
               cdf(Normal(), bsm_q1(v, ttm, sigma, r, gross_delta, xi, k)) +
               exp((- bsm_a(sigma, r, gross_delta) - bsm_zhat(sigma, r, gross_delta, xi, k)) * v) *
               cdf(Normal(), bsm_q2(v, ttm, sigma, r, gross_delta, xi, k))
    end

    # Bond Price
    function zhi_bond_price(vb, v, ttm, m, c, p, sigma, r, gross_delta, xi, k, alpha, pi)
        return (c / rdisc(r, xi, k)) + exp(-rdisc(r, xi, k) * ttm) * (p - c / rdisc(r, xi, k)) *
               (1 - bsm_F(v, ttm, sigma, r, gross_delta)) + (alpha * vb / m - c / rdisc(r, xi, k)) * bsm_G(v, ttm, sigma, r, gross_delta, xi, k)
    end


    # Credit-Risk-Free Bond Price
    function rfbond_price(ttm, c, p, r, xi, k)
        return ((c / rdisc(r, xi, k)) +
                (p - c / rdisc(r, xi, k)) * exp(- rdisc(r, xi, k) * ttm))
    end


    ######### vb-Specific Auxiliary Functions #########
    function bsm_b(x, m, sigma, r, gross_delta, xi, k)
        return (1 / (bsm_z(sigma, r, gross_delta) + x)) *
            exp(- rdisc(r, xi, k) * m) *
            (cdf(Normal(), x * sigma * sqrt(m)) -
             exp(r * m) * cdf(Normal(), - bsm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
    end

    function bsm_B(x, m, sigma, r, gross_delta)
        return (1 / (bsm_z(sigma, r, gross_delta) + x)) *
            (cdf(Normal(), x * sigma * sqrt(m)) -
                exp(.5 * (bsm_z(sigma, r, gross_delta)^2 - x^2) * sigma^2 * m) *
                cdf(Normal(), - bsm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
    end

    function bsm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi)
        return ((1 - pi) * (m * c) +
                (1 - exp(- rdisc(r, xi, k) * m)) *
                (p - c / rdisc(r, xi, k))) / bsm_eta(sigma, r, gross_delta)
    end

    function bsm_numerator2(m, c, p, sigma, r, gross_delta, xi, k)
        return (p - c / rdisc(r, xi, k)) *
               (bsm_b(- bsm_a(sigma, r, gross_delta), m, sigma, r, gross_delta, xi, k) +
                bsm_b(bsm_a(sigma, r, gross_delta), m, sigma, r, gross_delta, xi, k)) +
               (c / rdisc(r, xi, k)) *
               (bsm_B(-bsm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta) +
                bsm_B(bsm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta))
    end

    function bsm_denominator(m, sigma, r, gross_delta, iota, xi, k, alpha)
        return (gross_delta - iota) / (bsm_eta(sigma, r, gross_delta) - 1) + (alpha / m) *
            (bsm_B(- bsm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta) +
             bsm_B(bsm_zhat(sigma, r, gross_delta, xi, k), m, sigma, r, gross_delta))
    end

    # Default Boundary
    function zhi_vb(m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        value = (bsm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi) +
                 bsm_numerator2(m, c, p, sigma, r, gross_delta, xi, k)) /
                 bsm_denominator(m, sigma, r, gross_delta, iota, xi, k, alpha)

        return max(convert(AbstractFloat, value), 1e-4)
    end

    ######### Equity-Specific Auxiliary Functions #########
    function bsm_gamma(sigma, r, gross_delta)
        return bsm_a(sigma, r, gross_delta) + bsm_z(sigma, r, gross_delta)
    end

    function bsm_k_fun(v, x, w, u, sigma, m)
        return exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
               exp(-u * v) * cdf(Normal(),(-v + (u - x) * (sigma^2) * m) /
                               (sigma * sqrt(m))) -
               exp(-(x + w) * v) * cdf(Normal(), (-v + w * (sigma^2) * m) /
                                     (sigma * sqrt(m)))
    end

    function bsm_K(v, x, w, u, sigma, m)
        return (cdf(Normal(), w * sigma * sqrt(m)) -
                exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
                cdf(Normal(),(u - x) * sigma * sqrt(m))) * exp(- u * v) +
               bsm_k_fun(v, x, w, u, sigma, m)
    end

    function bsm_A(v, y, sigma, r, gross_delta, m)
        return 1 / (bsm_z(sigma, r, gross_delta) - y) *
            (bsm_K(v, bsm_a(sigma, r, gross_delta), y, bsm_gamma(sigma, r, gross_delta), sigma, m) +
             bsm_k_fun(v, bsm_a(sigma, r, gross_delta), - y, - bsm_eta(sigma, r, gross_delta), sigma, m)) +
            1 / (bsm_z(sigma, r, gross_delta) + y) *
            (bsm_K(v, bsm_a(sigma, r, gross_delta), - y, bsm_gamma(sigma, r, gross_delta), sigma, m) +
             bsm_k_fun(v,bsm_a(sigma, r, gross_delta), y, - bsm_eta(sigma, r, gross_delta), sigma, m))
    end

    # Equity Value
    function zhi_eq(v, m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        zhivb = zhi_vb(m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        net_delta = (gross_delta - iota)

        return (net_delta / (r - rgrow(r, gross_delta))) * zhivb * exp(v) -
               (net_delta * zhivb / (bsm_z(sigma, r, gross_delta) * sigma^2)) *
            (exp(-bsm_gamma(sigma, r, gross_delta) * v) /
             (bsm_gamma(sigma, r, gross_delta) + 1)) -
            (((1 - pi) * (c * m) + (1 - exp(- rdisc(r, xi, k) * m)) *
              (p - c / rdisc(r, xi, k))) / (bsm_z(sigma, r, gross_delta) * sigma^2)) *
            (1 / bsm_eta(sigma, r, gross_delta) +
                (1 - exp(- bsm_gamma(sigma, r, gross_delta) * v)) /
             bsm_gamma(sigma, r, gross_delta)) +
            1 / (bsm_z(sigma, r, gross_delta) * sigma^2) *
            (exp(- rdisc(r, xi, k) * m) * (p - c / rdisc(r, xi, k)) *
             bsm_A(v, bsm_a(sigma, r, gross_delta), sigma, r, gross_delta, m) -
             (alpha * zhivb / m -
              c / rdisc(r, xi, k)) *
             bsm_A(v, bsm_zhat(sigma, r, gross_delta, xi, k), sigma, r, gross_delta, m))
    end


    # ######## Pre-Volatility Shock Bond Pricing Functions #########
    # Pre-Volatility Shock Auxiliary Functions
    function rdisc_pvs(r, xi, k, _lambda)
        return r + xi * k + _lambda
    end

    function zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)
        return sqrt(bsm_a(sigma, r, gross_delta)^2 * sigma^4 +
                          2 * rdisc_pvs(r, xi, k, _lambda) * sigma^2) / sigma^2
    end

    function q1_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
        return (-v - zhat_pvs(sigma, r, gross_delta, xi, k, _lambda) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function q2_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
        return (-v + zhat_pvs(sigma, r, gross_delta, xi, k, _lambda) * sigma^2 * ttm) / (sigma * sqrt(ttm))
    end

    function G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
        return exp((- bsm_a(sigma, r, gross_delta) + zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)) * v) *
                                                cdf(Normal(), q1_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)) +
               exp((- bsm_a(sigma, r, gross_delta) - zhat_pvs(sigma, r, gross_delta, xi, k, _lambda)) * v) *
                                                cdf(Normal(), q2_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda))
    end

    function psi_v_td(vt, v, ttm, sigma, r, gross_delta)
        return (cdf(Normal(), (-v + vt) / (sigma * sqrt(ttm)) + bsm_a(sigma, r, gross_delta) * sigma * sqrt(ttm))
                - exp(-2 * bsm_a(sigma, r, gross_delta) * vt) *
                cdf(Normal(), (-v - vt) / (sigma * sqrt(ttm)) + bsm_a(sigma, r, gross_delta) * sigma * sqrt(ttm)))
    end

    function dv_psi_v_td(vt, v, ttm, sigma, r, gross_delta)
        return (-1./(sigma * sqrt(ttm))) *
               (pdf(Normal(), (-v + vt) / (sigma * sqrt(ttm)) + bsm_a(sigma, r, gross_delta) * sigma * sqrt(ttm))
                - exp(-2 * bsm_a(sigma, r, gross_delta) * vt) *
                pdf(Normal(), (-v - vt) / (sigma * sqrt(ttm)) + bsm_a(sigma, r, gross_delta) * sigma * sqrt(ttm)))
    end

    # Market Value of (Credit-Risk-Free) Debt at Upper Barrier:
    function rf_debt(m, c, p, r, xi, k)
        return ((c / rdisc(r, xi, k)) * (m - (1 - exp(-(rdisc(r, xi, k)) * m)) / (rdisc(r, xi, k))) +
                p * (1 - exp(-(rdisc(r, xi, k)) * m)) / (rdisc(r, xi, k)))
    end

    # Notice Psi = 1- F and does not depend on lambda
    # Prior to the volatility shock, G_pvs takes rdisc = r + xi*k + lambda!!!!
    function on_default_payoff(v, vb, ttm, m, c, p, r, xi, k, alpha)
        return (alpha * vb * exp(v) / m) * (rf_debt(m, c, p, r, xi, k) > (alpha * vb * exp(v) / m)) +
                (c / rdisc(r, xi, k) + exp(-rdisc(r, xi, k) * ttm) * (p - c / rdisc(r, xi, k))) *
                (rf_debt(m, c, p, r, xi, k) <= (alpha * vb * exp(v) / m))
    end

    function coupon_pv_fun(v, ttm, c, sigma, r, gross_delta, xi, k, _lambda)
        return (c / rdisc_pvs(r, xi, k, _lambda)) * ((1 - exp(- rdisc_pvs(r, xi, k, _lambda) * ttm) *
                                                     (1 - bsm_F(v, ttm, sigma, r, gross_delta))) -
                                                     G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda))
    end

    function default_payoff_pv_fun(v, vb, ttm, m, sigma, r, gross_delta, xi, k, alpha, _lambda)
        return (alpha * vb / m) * G_pvs(v, ttm, sigma, r, gross_delta, xi, k, _lambda)
    end

    function maturity_payoff_pv_fun(v, ttm, p, sigma, r, gross_delta, xi, k, _lambda)
        return p * exp(- rdisc_pvs(r, xi, k, _lambda) * ttm) * (1 - bsm_F(v, ttm, sigma, r, gross_delta))
    end

    function no_vol_shock_cf_pv(v, vb, ttm, m, c, p, sigma, r, gross_delta, xi, k, alpha, _lambda)
        return coupon_pv_fun(v, ttm, c, sigma, r, gross_delta, xi, k, _lambda) +
               default_payoff_pv_fun(v, vb, ttm, m, sigma, r, gross_delta, xi, k, alpha, _lambda) +
               maturity_payoff_pv_fun(v, ttm, p, sigma, r, gross_delta, xi, k, _lambda)
    end

    # Notice that time-to-maturity becomes ttm-u and maturity becomes ttm:
    function vol_shock_cf_integrand(vt, v_prior, v_post, vb_prior, vb_post, u, ttm, c, p,
                                   sigma_prior, sigma_post, r, gross_delta, xi, k, alpha, pi, _lambda)
        return _lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
               zhi_bond_price(vb_post, v_post, ttm - u, ttm, c, p, sigma_post, r, gross_delta, xi, k, alpha, pi) *
               (-dv_psi_v_td(vt, v_prior, u, sigma_prior, r, gross_delta) / (vb_prior * exp(v_prior)))
    end

# #########

    function approx_bond_int(vt, v_upper_thresh, ttm, c, p, r, xi, k, _lambda, sigmal, gross_delta)
        return _lambda * exp(-(rdisc(r, xi, k) + _lambda) * ttm) *
               rfbond_price(ttm, c, p, r, xi, k) *
               psi_v_td(vt, v_upper_thresh, ttm, sigmal, r, gross_delta)
    end

    function bond_pr_vol_cf(vt, ttm, vmax, vb, vbh,
                            c, p, r, gross_delta,
                            xi, k, alpha, pi,
                            _lambda, sigmal, sigmah)
        tauN = 10^3
        vN = 10^3
        dt1 = (ttm - 1e-4)/tauN
        dv = vmax/vN
        cf1_mat = @spawn [vol_shock_cf_integrand(vt, v, v - log(vbh / vb),
                                                 vb, vbh, u, ttm,
                                                 c, p, sigmal, sigmah, r,
                                                 gross_delta, xi, k,
                                                 alpha, pi, _lambda) *
                          vb * exp(v) * (v > 0) for u=linspace(1e-4, ttm, tauN),
                                                    v=linspace(0.0, vmax, vN)]

        N = 10^4
        dt2 = ttm/N
        cf2_vec = @spawn [approx_bond_int(vt, vmax, tau,  c, p, r,
                                           xi, k, _lambda, sigmal, gross_delta)
                          for tau=linspace(dt2, ttm, N)]
        return sum(fetch(cf1_mat)) * dt1 * dv + sum(fetch(cf2_vec)) * dt2
    end

    # Objective is to find V such that the value of the
    # newly-issued (tau = m) risky bond price when sigma = sigmah
    # is sufficiently close to the credit-risk-free bond price.
    function get_bond_vmax(V0, m, c, p, sigmah, r,
                            gross_delta, iota, xi, k,
                            alpha, pi, print=false)
        bondVmax = 1.25*V0
        vb = zhi_vb(m, c, p, sigmah, r, gross_delta, iota, xi, k, alpha, pi)

        vmax = log(bondVmax/vb)
        rfbond = rfbond_price(m, c, p, r, xi, k)
        bondpr = zhi_bond_price(vb, vmax, m, m, c, p, sigmah, r,
                                    gross_delta, xi, k, alpha, pi)

        per_diff = (rfbond - bondpr) / rfbond

        cond = per_diff > 1e-3
        while cond
            bondVmax = 1.025 * bondVmax
            vmax = log(bondVmax / vb)
            bondpr = zhi_bond_price(vb, vmax, m, m, c, p, sigmah, r,
                                gross_delta, xi, k, alpha, pi)
            per_diff = (rfbond - bondpr) / rfbond
            cond = per_diff > 1e-3
        end

        if print
            println(string("Bond Percentage Difference: ", per_diff))
        end

        return bondVmax
    end


    function bond_pr_main(vt, ttm, vb, vbh,
                          m, c, p,
                          r, gross_delta, iota,
                          xi, k, alpha, pi,
                          _lambda, sigmal, sigmah, vmax=NaN)

        if vt <= 0
            return on_default_payoff(vt, vb, ttm,
                                     m, c, p, r, xi,
                                     k, alpha)
        else
            Vt = vb*exp(vt)
            # Compute vmax
            if isnan(vmax)
                bondVmax = get_bond_vmax(Vt, m, c, p, sigmah, r,
                                         gross_delta, iota, xi,
                                         k, alpha, pi)
                vmax = log(bondVmax/vb)
            end

            if vt > vmax
                return rfbond_price(ttm, c, p, r, xi, k)
            else
                cf0 = no_vol_shock_cf_pv(vt, vb, ttm,
                                         m, c, p, sigmal,
                                         r, gross_delta,
                                         xi, k, alpha, _lambda)

                cf1 = bond_pr_vol_cf(vt, ttm, vmax, vb, vbh,
                                     c, p, r, gross_delta,
                                     xi, k, alpha, pi,
                                     _lambda, sigmal, sigmah)

                return min(cf0 + cf1, rfbond_price(ttm, c, p, r, xi, k))
            end
        end
    end

    function get_debt_price(vt, vb, vbh, ttm_grid_size,
                            m, c, p,
                            r, gross_delta, iota,
                            xi, k, alpha, pi,
                            _lambda, sigmal, sigmah)
        t0 = .0001

        # Grid for Interpolation
        ttm_grid = linspace(t0, m, ttm_grid_size)

        # Since vmax does not change with maturity, compute
        Vt = vb*exp(vt)
        # it only once:
        bondVmax = get_bond_vmax(Vt, m, c, p, sigmah, r,
                                 gross_delta, iota, xi,
                                 k, alpha, pi, true)
        vmax = log(bondVmax/vb)

        # Compute Bond Prices for Different Maturities
        tmp = @spawn [bond_pr_main(vt, ttm, vb, vbh,
                                    m, c, p,
                                    r, gross_delta, iota,
                                    xi, k, alpha, pi,
                                    _lambda, sigmal, sigmah, vmax) for ttm=ttm_grid]
        ttmvec = fetch(tmp)

        # Interpolate
        itp = interpolate(ttmvec, BSpline(Cubic(Line())), OnGrid())
        sitp = Interpolations.scale(itp, ttm_grid) # Scale

        # Refined Grid
        ttm_grid_refined_size = 10^4
        ttm_grid_refined = linspace(t0, m, ttm_grid_refined_size)
        dt = (m - t0)/ttm_grid_refined_size

        # Integrate
        return (1/m) * sum([sitp[x] for x=ttm_grid_refined]) * dt
    end

    function pvvb_bond_surf(V0, pmin, pmax, pN, vN,
                            vbmin, vbmax, vbN,
                            m, c, p,
                            r, gross_delta, iota,
                            xi, k,
                            alpha, pi,
                            _lambda, sigmal, sigmah)

        bondVmax = get_bond_vmax(V0, m, c, p, sigmah, r,
                                 gross_delta, iota, xi, k,
                                 alpha, pi)

        pgrid = linspace(pmin, pmax, pN)
        vgrid = linspace(0, log(bondVmax/vbmin), vN)
        vbgrid = linspace(vbmin, vbmax, vbN)

        cube_future = @spawn [bond_pr_main(v, m, vb, vb, m, c, p,
                                 r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                                 sigmal, sigmah, log(bondVmax/vb)) for p=pgrid, v=vgrid, vb=vbgrid]

        cube = fetch(cube_future)
        save("data.jld", "data", cube)
        return cube
    end

    function tmp()
        return 3
    end

# End Module:
end
