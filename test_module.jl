module test

    function bond_pr_progress(v, vb, m, c, p,
                             r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                             sigmal, sigmah, bondVmax; progress=nothing)
        if progress != nothing
            next!(progress)
        end
        return bond_pr_main(v, m, vb, vb, m, c, p,
                                 r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                                 sigmal, sigmah, log(bondVmax/vb))
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

        # cube_future = @spawn [bond_pr_main(v, m, vb, vb, m, c, p,
        #                          r, gross_delta, iota,xi, k, alpha, pi, _lambda,
        #                          sigmal, sigmah, log(bondVmax/vb)) for p=pgrid, v=vgrid, vb=vbgrid]

        p = Progress(length(pgrid) * length(vgrid) * length(vbgrid), 1)
        cube_future = @spawn [bond_pr_progress(v, vb, m, c, p,
                                 r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                                 sigmal, sigmah, bondVmax; progress=p) for p=pgrid, v=vgrid, vb=vbgrid]
       return fetch(cube_future)

    end

    # ########################################################################

    # Notice that time-to-maturity becomes ttm-u and maturity becomes ttm:
    function vol_shock_cf_integrand_Psi(vt, v_prior, vb_prior, u, ttm,
                                    sigma_prior, sigma_post,
                                    r, gross_delta, xi, k, _lambda)
        return _lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                   (-dv_psi_v_td(vt, v_prior, u, sigma_prior, r, gross_delta))
    end

    # Notice that time-to-maturity becomes ttm-u and maturity becomes ttm:
    function vol_shock_cf_integrand_FPsi(vt, v_prior, vb_prior, vb_post, u, ttm,
                                     sigma_prior, sigma_post, r, gross_delta,
                                     xi, k, _lambda)
        return _lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                   (exp(-rdisc(r, xi, k) * (ttm - u)) * (1 - bsm_F(v_prior - log(vb_post/vb_prior),
                                                                   ttm - u, sigma_post, r, gross_delta))) *
                   (-dv_psi_v_td(vt, v_prior, u, sigma_prior, r, gross_delta)) *
                   (v_prior - log(vb_post/vb_prior) > 0)
    end

    function vol_shock_cf_integrand_GPsi(vt, v_prior, vb_prior, vb_post, u, ttm,
                                     sigma_prior, sigma_post, r, gross_delta,
                                     xi, k, _lambda)
        return _lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                   bsm_G(v_prior - log(vb_post/vb_prior), ttm - u, sigma_post, r, gross_delta, xi, k) *
                   (-dv_psi_v_td(vt, v_prior, u, sigma_prior, r, gross_delta)) *
                   (v_prior - log(vb_post/vb_prior) > 0)
    end

    function bond_int_terms(vt, v, vbl, vbh,
                        u, ttm, vmax,
                        sigmal, sigmah,
                        r, gross_delta,
                        xi, k, _lambda)

        tauN = 10^3
        vN = 10^3
        dt = (ttm - 1e-4)/tauN
        dv = vmax/vN

        ugrid = linspace(1e-4, ttm, tauN)
        vgrid = linspace(0, vmax, vN)


        cube_Psi_future = @spawn [vol_shock_cf_integrand_Psi(vt, v, vbl,
                                                             u, ttm,
                                                             sigmal, sigmah,
                                                             r, gross_delta,
                                                             xi, k, _lambda)
                                  for u=ugrid, v=vgrid]
        cube_Psi_Int = sum(fetch(cube_Psi_future)) * dt * dv
        # println(string("Cube_Psi_Int: ", cube_Psi_Int))

        cube_FPsi_future = @spawn [vol_shock_cf_integrand_FPsi(vt, v, vbl, vbh,
                                                               u, ttm,
                                                               sigmal, sigmah,
                                                               r, gross_delta,
                                                               xi, k, _lambda)
                                   for u=ugrid, v=vgrid]
        cube_FPsi_Int = sum(fetch(cube_FPsi_future)) * dt * dv
        # println(string("Cube_FPsi_Int: ", cube_FPsi_Int))


        cube_GPsi_future = @spawn [vol_shock_cf_integrand_GPsi(vt, v, vbl, vbh,
                                                               u, ttm,
                                                               sigmal, sigmah,
                                                               r, gross_delta,
                                                               xi, k, _lambda)
                                   for u=ugrid, v=vgrid]
        cube_GPsi_Int = sum(fetch(cube_GPsi_future)) * dt * dv
        # println(string("Cube_GPsi_Int: ", cube_GPsi_Int))

        return (cube_Psi_Int, cube_FPsi_Int, cube_GPsi_Int)
    end

end
