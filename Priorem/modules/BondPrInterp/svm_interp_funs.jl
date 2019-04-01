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


function grid_creator(z0, z1, N)
	dz = (z1 - z0) / float(N)
	if z0^2 < 1e-6
		# return dz, linspace(z0+dz, z1, N)
		return dz, range(z0+dz, stop=z1, length=N)
	else
		# return dz, linspace(z0, z1, N)
		return dz, range(z0, stop=z1, length=N) 
	end
end


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

    tmp = fetch(f21_integrand_vec)
    itp = interpolate(tmp, BSpline(Cubic(Line(Interpolations.OnGrid()))))
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

    tmp = interpolate(fetch(f22_integrand_vec), BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f22_integrand_sitp = Interpolations.scale(tmp, ugrid)

    du2, ugrid2 = grid_creator(1e-4, ugrid[end], 10^4) #  grid_creator(1e-4, ugrid[end], 10^4) 

    return sum([f22_integrand_sitp(x) for x in ugrid2]) * du2
end



# function f0_integrand(vt, ttm, u, vmax, _lambda, sigma,
# 		      c, p, r, gross_delta, xi, kappa)
# 
# 	return _lambda * exp(-rdisc_pvs(r, xi, kappa, _lambda) * u) *
# 	       rfbond_price(ttm - u, c, p, r, xi, kappa) *
#                psi_v_td(vt, vmax, u, sigma, r, gross_delta) 
# end


# function f0_int(vt, vmax, _lambda, sigma, 
#                 c, p, r, gross_delta, xi, kappa; 
# 		ttm=1.0, N=10^2, ugrid=nothing)
# 
#     if ugrid == nothing
# 	du, ugrid = grid_creator(1e-4, ttm, N)
#     else
#         ttm = ugrid[end]
# 	du = ugrid[2] - ugrid[1]
#     end
# 
#     f0_integrand_vec  = @spawn [f0_integrand(vt, ttm, u, vmax, _lambda, sigma,
# 				      c, p, r, gross_delta, xi, kappa)
# 				for u=ugrid]
# 
#     tmp = interpolate(fetch(f0_integrand_vec), BSpline(Cubic(Line(OnGrid()))))
#     f0_integrand_sitp = Interpolations.scale(tmp, ugrid)
# 
#     du2, ugrid2 = grid_creator(1e-4, ugrid[end], 10^4) 
# 
#     return sum([f0_integrand_sitp[x] for x=ugrid2]) * du2
# end


# function f0_interp(vmax, _lambda, sigma, 
#                    c, p, r, gross_delta, xi, kappa, 
# 		   vtgrid, tau_grid; uN=1e3)
# 
# 	if length(vtgrid) > 2 && length(tau_grid) > 2
# 		f0_int_vec = fetch(@spawn [f0_int(vt, vmax, _lambda, sigma, 
#                                     c, p, r, gross_delta, xi, kappa; 
# 				    ttm=tau, N=uN) for vt=vtgrid, tau=tau_grid])
# 		f0_int_itp = interpolate(f0_int_vec, BSpline(Cubic(Line(OnGrid()))))
# 		return Interpolations.scale(f0_int_itp, vtgrid, tau_grid)
# 	elseif length(tau_grid) > 2 
# 		f0_int_vec = fetch(@spawn [f0_int(vtgrid[1], vmax, _lambda, sigma, 
#                                     c, p, r, gross_delta, xi, kappa; 
# 				    ttm=tau, N=uN) for tau=tau_grid])
# 		f0_int_itp = interpolate(f0_int_vec, BSpline(Cubic(Line(OnGrid()))))
# 		return Interpolations.scale(f0_int_itp, tau_grid)
# 	else	
# 		f0_int_vec = fetch(@spawn [f0_int(vt, vmax, _lambda, sigma, 
#                                     c, p, r, gross_delta, xi, kappa; 
# 				    ttm=tau_grid[1], N=uN) for vt=vtgrid])
# 		f0_int_itp = interpolate(f0_int_vec, BSpline(Cubic(Line(OnGrid()))))
# 		return Interpolations.scale(f0_int_itp, vtgrid)
# 	end
# end


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


function f11v_int(vt, _lambda, sigmal,
          	r, gross_delta, xi, kappa;
 		vmax=1.2, vN=10^3, vgrid=nothing,
    	        ttm=1.0, uN=10^3, ugrid=nothing)

    if vgrid == nothing
        dv, vgrid = grid_creator(1e-4, vmax, vN)
    end

    if ugrid == nothing
         du, ugrid = grid_creator(1e-4, ttm, uN) 
    end

    return ugrid, fetch(@spawn [g11_int(vt, vgrid, u, _lambda, sigmal,
			r, gross_delta, xi, kappa) for u in ugrid]) 
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
                              for u=ugrid, v in vgrid]
   
    return sum(fetch(f12_integrand)) * du * dv
	
#     f2_itp = interpolate(fetch(f2_integrand), BSpline(Cubic(Line(OnGrid()))))
#     f2_sitp = Interpolations.scale(f2_itp, ugrid, vgrid)
# 
#     ugrid2 = linspace(ugrid[1], ugrid[end], 10^3)
#     vgrid2 = linspace(vgrid[1], vgrid[end], 10^3)
#     du2 = ugrid2[2] - ugrid2[1]
#     dv2 = vgrid2[2] - vgrid2[1]
# 
#     tmp = fetch(@spawn [f2_sitp[u, v] for u=ugrid2, v=vgrid2])
# 
    return sum(tmp) * du2 * dv2

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

