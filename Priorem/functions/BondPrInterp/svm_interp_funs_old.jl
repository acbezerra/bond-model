# SVM Interpolation Functions

    function bond_res_Psi(vt, m, ugrid, vmax, _lambda, sigma, 
			  c, p, r, gross_delta, xi, kappa)
			 
	resPsi_future = @spawn [_lambda * exp(-rdisc_pvs(r, xi, kappa, _lambda) * u) *
                                rfbond_price(m - u, c, p, r, xi, kappa) *
				psi_v_td(vt, v_upper_thresh, u,
                                         sigma, r, gross_delta) for u=ugrid]
        return sum(fetch(resPsi_future))
    end


    function bond_int_Psi(vt, # ttm,  # vmax,
                          sigmal,
                          r, gross_delta,
                          xi, k, _lambda, ugrid, vgrid)

        # tauN = 10^3
        # vN = 10^3
        # dt = (ttm - 1e-4)/tauN
        # dv = vmax/vN
        #
        # ugrid = linspace(1e-4, ttm, tauN)
        # vgrid = linspace(0, vmax, vN)

        cube_Psi_future = @spawn [_lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                                  (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta))
                                  for u=ugrid, v=vgrid]
        return sum(fetch(cube_Psi_future)) # * dt * dv
    end

    function bond_int_FPsi(vt, vbhl, # ttm, # vmax,
                            sigmal, sigmah,
                            r, gross_delta,
                            xi, k, _lambda, ugrid, vgrid)

        # vbhl is the ratio vbh/vbl

        ttm = ugrid[end]
        cube_FPsi_future = @spawn [_lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                                   (exp(-rdisc(r, xi, k) * (ttm - u)) * (1 - cvm_F(v - log(vbhl),
                                                               ttm - u, sigmah, r, gross_delta))) *
                                   (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                                   (v - log(vbhl) > 0)
                                  for u=ugrid, v=vgrid]
        return sum(fetch(cube_FPsi_future)) # * dt * dv
    end

    function bond_int_GPsi(vt, vbhl, #ttm, #  vmax,
                            sigmal, sigmah,
                            r, gross_delta,
                            xi, k, _lambda, ugrid, vgrid)

        # vbhl is the ratio vbh/vbl

        ttm = ugrid[end]
        cube_GPsi_future = @spawn [_lambda * exp(- rdisc_pvs(r, xi, k, _lambda) * u) *
                                   cvm_G(v - log(vbhl), ttm - u,
                                         sigmah, r, gross_delta, xi, k) *
                                   (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                                   (v - log(vbhl) > 0)
                                  for u=ugrid, v=vgrid]
        return sum(fetch(cube_GPsi_future)) # * dt * dv
    end


    function interpd_bondpr_funs(vtgrid, vbhlgrid, ugrid, vgrid,
                                 sigmal, sigmah,
                                 r, gross_delta,
                                 xi, kappa, _lambda)

        # ttm = ugrid[end]
        dt = (ugrid[end] - ugrid[1])/length(ugrid)
        dv = (vgrid[end] - vgrid[1])/length(vgrid)


        v_upper_thresh = vgrid[end]
        ttm = ugrid[end]
        N = 10^4
        dt2 = ttm/N
        ugrid2=linspace(dt2, ttm, N)
        resPsi_Future = @spawn [bond_res_Psi(vt, v_upper_thresh,
                                             sigmal,
                                             r, gross_delta,
                                             xi, kappa, _lambda,
                                             ugrid2) * dt2
                                    for vt=vtgrid]
        resPsi = fetch(resPsi_Future)

        cubePsi_Future = @spawn [bond_int_Psi(vt, sigmal,
                                            r, gross_delta,
                                            xi, kappa, _lambda,
                                            ugrid, vgrid) * dt * dv
                                 for vt=vtgrid]
        cubePsi = fetch(cubePsi_Future)

        cubeFPsi_Future = @spawn [bond_int_FPsi(vt, vbhl,
                                        sigmal, sigmah,
                                        r, gross_delta,
                                        xi, kappa, _lambda,
                                        ugrid, vgrid) * dt * dv
                                  for vbhl=vbhlgrid, vt=vtgrid]
        cubeFPsi = fetch(cubeFPsi_Future)

        cubeGPsi_Future = @spawn [bond_int_GPsi(vt, vbhl,
                                        sigmal, sigmah,
                                        r, gross_delta,
                                        xi, kappa, _lambda,
                                        ugrid, vgrid) * dt * dv
                                  for vbhl=vbhlgrid, vt=vtgrid]
        cubeGPsi = fetch(cubeGPsi_Future)

        return resPsi, cubePsi, cubeFPsi, cubeGPsi
    end

    function bond_pr_interp(st, Vt, vbl, ttm, vmax,
                        resPsi, iPsi, iFPsi, iGPsi)

        # ####################################
        # ######## Extract Parameters ########
        # ####################################
        # Capital Structure
        m = st.pm.m
        c = st.c
        p = st.p

        # Default
        alpha = st.pm.alpha
        pi = st.pm.pi

        # Dynamics
        r = st.pm.r
        gross_delta = st.pm.gross_delta
        iota = st.pm.iota

        # Liquidity
        xi = st.pm.xi
        k = st.pm.kappa

        # Volatility
        _lambda = st.pm.lambda
        sigmal = st.pm.sigmal
        sigmah = st.pm.sigmah
        # ####################################

        vt = log(Vt/vbl)

        # Default Barrier
        vbh = zhi_vb(m, c, p, sigmah, r,
                        gross_delta, iota, xi, k, alpha, pi)

        if vt < 0
                return on_default_payoff(vt, vbl, ttm,
                                         m, c, p, r, xi,
                                         k, alpha)
        else
            if vt > vmax
                return rfbond_price(ttm, c, p, r, xi, k)
            else
                # Maturity or Default prior to Volatility Shock:
                cf0 = no_vol_shock_cf_pv(vt, vbl, ttm,
                                         m, c, p, sigmal,
                                         r, gross_delta,
                                         xi, k, alpha, _lambda)

                # Volatility Shock Prior to Maturity:
                cf1 = c/rdisc(r, xi, k) * iPsi[vt] +
                      (p - c/rdisc(r, xi, k)) * iFPsi[vbh/vbl, vt] +
                      (alpha * vbh/m - c/rdisc(r, xi, k)) * iGPsi[vbh/vbl, vt]
                # println(cf1)
                cf2 = rfbond_price(ttm, c, p, r, xi, k) * resPsi[vt]
                # println(cf2)

                return min(cf0 + cf1 + cf2, rfbond_price(ttm, c, p, r, xi, k))
            end
        end
    end

