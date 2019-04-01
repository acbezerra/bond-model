# SVM Interpolation Functions

# Functions from the Analytic Functions Module:
# rdisc
# rdisc_pvs
# rfbond_price
# psi_v_td
# dv_psi_v_td
# cvm_F
# cvm_G
# zhi_vb
# on_default_payoff
# no_vol_shock_cf_pv



# ##################################################################
# ######################## J2: f21, f22 ############################
# ################################################################## 
function f21_integrand(vt, ttm, u, vmax, _lambda, sigma,
		      r, gross_delta, xi, kappa)

    return _lambda * exp(-rdisc_pvs(r, xi, kappa, _lambda) * u) *
            psi_v_td(vt, vmax, u, sigma, r, gross_delta)
end

function f22_integrand(vt, ttm, u, vmax, _lambda, sigma,
		      r, gross_delta, xi, kappa)

    return  _lambda * exp(-rdisc(r, xi, kappa) * ttm + _lambda * u) *
            psi_v_td(vt, vmax, u, sigma, r, gross_delta)
end

function f21_int(vt, vmax, _lambda, sigma, 
                r, gross_delta, xi, kappa; 
		ttm=1.0, N=10^2, ugrid=nothing)

    if ugrid == nothing
	du, ugrid = grid_creator(1e-4, ttm, N)
    else
        ttm = ugrid[end]
	du = ugrid[2] - ugrid[1]
    end

    f21_integrand_vec  = @spawn [f21_integrand(vt, ttm, u, vmax, _lambda, sigma,
				      r, gross_delta, xi, kappa)
				for u in ugrid]

    itp = interpolate(fetch(f21_integrand_vec), BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f21_integrand_sitp = Interpolations.scale(itp, ugrid)

    du2, ugrid2 = grid_creator(ugrid[1], ugrid[end], 10^4) # grid_creator(1e-4, ugrid[end], 10^4) 

    return sum([f21_integrand_sitp(x) for x in ugrid2]) * du2
end

function f22_int(vt, vmax, _lambda, sigma, 
                r, gross_delta, xi, kappa; 
		ttm=1.0, N=10^2, ugrid=nothing)

    if ugrid == nothing
	du, ugrid = grid_creator(1e-4, ttm, N)
    else
        ttm = ugrid[end]
	du = ugrid[2] - ugrid[1]
    end

    f22_integrand_vec  = @spawn [f22_integrand(vt, ttm, u, vmax, _lambda, sigma,
				      r, gross_delta, xi, kappa)
				for u=ugrid]

    itp = interpolate(fetch(f22_integrand_vec), BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f22_integrand_sitp = Interpolations.scale(itp, ugrid)

    du2, ugrid2 = grid_creator(ugrid[1], ugrid[end], 10^4) #  grid_creator(1e-4, ugrid[end], 10^4) 

    return sum([f22_integrand_sitp(x) for x in ugrid2]) * du2
end


# ##################################################################
# ########################### f1_int ###############################
# ################################################################## 

function g11_int(vt, vgrid, u, _lambda, sigmal,
          	r, gross_delta, xi, kappa)
     dv = vgrid[2] - vgrid[1] 
     f11_integrand_v = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) * 
			 (- dv_psi_v_td(vt, v, u, sigmal, r, gross_delta))
			 for v in vgrid]
     return sum(fetch(f11_integrand_v)) * dv
end


function f11v_int(vt, ugrid, vgrid, _lambda,
                  sigmal, r, gross_delta, xi, kappa)

    du = ugrid[2] - ugrid[1]
    f11v_integrand = @spawn [g11_int(vt, vgrid, u, _lambda, sigmal,
			          r, gross_delta, xi, kappa) for u in ugrid]

    return sum(fetch(f11v_integrand)) * du
end


function integral_u(f, u; zN=10^4)

	dz = (u-1e-5)/zN
	# zgrid = linspace(dz, u, zN) 
	zgrid = range(dz, stop=u, length=zN)
	return sum(fetch(@spawn [f[z] for z in zgrid])) * dz
end


# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 

function f12_int(vt, vbhl, ugrid, vgrid,
                _lambda, sigmal, sigmah,
    	    r, gross_delta, xi, kappa)

    # vbhl is the ratio vbh/vbl

    ttm = ugrid[end]
    du = ugrid[2] - ugrid[1] 
    dv = vgrid[2] - vgrid[1] 

    f12_integrand = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) *
                               (exp(-rdisc(r, xi, kappa) * (ttm - u)) * 
    			    (1 - cvm_F(v - log(vbhl), ttm - u, sigmah, r, gross_delta))) *
                               (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                               (v - log(vbhl) > 0)
                              for u in ugrid, v in vgrid]
   
    return sum(fetch(f12_integrand)) * du * dv
	
end


# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 

function f13_int(vt, vbhl, ugrid, vgrid, 
                _lambda, sigmal, sigmah,
                r, gross_delta,
                xi, kappa)

    # vbhl is the ratio vbh/vbl

    ttm = ugrid[end]
    du = ugrid[2] - ugrid[1] 
    dv = vgrid[2] - vgrid[1] 

    f13_integrand = @spawn [_lambda * exp(- rdisc_pvs(r, xi, kappa, _lambda) * u) *
                               cvm_G(v - log(vbhl), ttm - u,
                                     sigmah, r, gross_delta, xi, kappa) *
                               (-dv_psi_v_td(vt, v, u, sigmal, r, gross_delta)) *
                               (v - log(vbhl) > 0)
                              for u in ugrid, v in vgrid]
    return sum(fetch(f13_integrand)) * du * dv
end

