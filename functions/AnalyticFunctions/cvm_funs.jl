# Analytic Functions for the CVM Model

    function rgrow(r, gross_delta)
        r - gross_delta
    end

    function rdisc(r, xi, k)
        r + xi * k
    end

    # Auxiliary Functions
    function cvm_a(sigma, r, gross_delta)
        (rgrow(r, gross_delta) - .5 * sigma ^ 2) / sigma ^ 2
    end

    function cvm_z(sigma, r, gross_delta)
        return (cvm_a(sigma, r, gross_delta)^2 * sigma^4 + 2 * r * sigma^2)^.5 / sigma^2
    end

    function cvm_zhat(sigma, r, gross_delta, xi, k)
        (cvm_a(sigma, r, gross_delta)^2 * sigma^4 +
                          2 * rdisc(r, xi, k) * sigma^2)^.5 / sigma^2
    end

    function cvm_eta(sigma, r, gross_delta)
        cvm_z(sigma, r, gross_delta) - cvm_a(sigma, r, gross_delta)
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
        return cdf(Normal(), cvm_h1(v, ttm, sigma, r, gross_delta)) +
               exp(-2 * cvm_a(sigma, r, gross_delta) * v) *
                cdf(Normal(), cvm_h2(v, ttm, sigma, r, gross_delta))
    end

    function cvm_G(v, ttm, sigma, r, gross_delta, xi, k)
        return exp((- cvm_a(sigma, r, gross_delta) + cvm_zhat(sigma, r, gross_delta, xi, k)) * v) *
               cdf(Normal(), cvm_q1(v, ttm, sigma, r, gross_delta, xi, k)) +
               exp((- cvm_a(sigma, r, gross_delta) - cvm_zhat(sigma, r, gross_delta, xi, k)) * v) *
               cdf(Normal(), cvm_q2(v, ttm, sigma, r, gross_delta, xi, k))
    end

    # Bond Price
    function cvm_bond_price(vb, v, ttm, m, c, p, sigma, r, gross_delta, xi, k, alpha, pi)
        return (c / rdisc(r, xi, k)) + exp(-rdisc(r, xi, k) * ttm) * (p - c / rdisc(r, xi, k)) *
               (1 - cvm_F(v, ttm, sigma, r, gross_delta)) + (alpha * vb / m - c / rdisc(r, xi, k)) * cvm_G(v, ttm, sigma, r, gross_delta, xi, k)
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
            (cdf(Normal(), x * sigma * sqrt(m)) -
             exp(r * m) * cdf(Normal(), - cvm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
    end

    function cvm_B(x, m, sigma, r, gross_delta)
        return (1 / (cvm_z(sigma, r, gross_delta) + x)) *
            (cdf(Normal(), x * sigma * sqrt(m)) -
                exp(.5 * (cvm_z(sigma, r, gross_delta)^2 - x^2) * sigma^2 * m) *
                cdf(Normal(), - cvm_z(sigma, r, gross_delta) * sigma * sqrt(m)))
    end

    function cvm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi)
        return ((1 - pi) * (m * c) +
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
    function zhi_vb(m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        value = (cvm_numerator1(m, c, p, sigma, r, gross_delta, xi, k, pi) +
                 cvm_numerator2(m, c, p, sigma, r, gross_delta, xi, k)) /
                 cvm_denominator(m, sigma, r, gross_delta, iota, xi, k, alpha)

        return max(convert(AbstractFloat, value), 1e-4)
    end

    ######### Equity-Specific Auxiliary Functions #########
    function cvm_gamma(sigma, r, gross_delta)
        return cvm_a(sigma, r, gross_delta) + cvm_z(sigma, r, gross_delta)
    end

    function cvm_k_fun(v, x, w, u, sigma, m)
        return exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
               exp(-u * v) * cdf(Normal(),(-v + (u - x) * (sigma^2) * m) /
                               (sigma * sqrt(m))) -
               exp(-(x + w) * v) * cdf(Normal(), (-v + w * (sigma^2) * m) /
                                     (sigma * sqrt(m)))
    end

    function cvm_K(v, x, w, u, sigma, m)
        return (cdf(Normal(), w * sigma * sqrt(m)) -
                exp(.5 * ((u - x)^2 - w^2) * (sigma^2) * m) *
                cdf(Normal(),(u - x) * sigma * sqrt(m))) * exp(- u * v) +
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
    function zhi_eq(v, m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        zhivb = zhi_vb(m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
        net_delta = (gross_delta - iota)

        return (net_delta / (r - rgrow(r, gross_delta))) * zhivb * exp(v) -
               (net_delta * zhivb / (cvm_z(sigma, r, gross_delta) * sigma^2)) *
            (exp(-cvm_gamma(sigma, r, gross_delta) * v) /
             (cvm_gamma(sigma, r, gross_delta) + 1)) -
            (((1 - pi) * (c * m) + (1 - exp(- rdisc(r, xi, k) * m)) *
              (p - c / rdisc(r, xi, k))) / (cvm_z(sigma, r, gross_delta) * sigma^2)) *
            (1 / cvm_eta(sigma, r, gross_delta) +
                (1 - exp(- cvm_gamma(sigma, r, gross_delta) * v)) /
             cvm_gamma(sigma, r, gross_delta)) +
            1 / (cvm_z(sigma, r, gross_delta) * sigma^2) *
            (exp(- rdisc(r, xi, k) * m) * (p - c / rdisc(r, xi, k)) *
             cvm_A(v, cvm_a(sigma, r, gross_delta), sigma, r, gross_delta, m) -
             (alpha * zhivb / m -
              c / rdisc(r, xi, k)) *
             cvm_A(v, cvm_zhat(sigma, r, gross_delta, xi, k), sigma, r, gross_delta, m))
    end


