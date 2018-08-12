
# ##################################################################
# ########################### f0_int ###############################
# ################################################################## 

function f0_fd(vtgrid, m, vmax, _lambda, sigma,
		     c, p, r, gross_delta, xi, kappa;
		     uN=1e3)
	
	du, ugrid = grid_creator(1e-4, m, uN)

	f0_fd_fut = fetch(@spawn [f0_int(vt, vmax, _lambda, sigma, 
                			c, p, r, gross_delta, xi, kappa; 
					ugrid=ugrid)
					for vt=vtgrid])

	f0_fd_itp = interpolate(f0_fd_fut, BSpline(Cubic(Line())), OnGrid())
	return Interpolations.scale(f0_fd_itp, vtgrid)
end


# ##################################################################
# ########################### f1_int ###############################
# ################################################################## 

function f1_fd(vtgrid, _lambda, sigmal,
		     r, gross_delta, xi, kappa;
	             vmax=1.2, vN=1e3,
		     ttm=1.0, uN=1e3)

        dv, vgrid = grid_creator(1e-4, vmax, vN)
        du, ugrid = grid_creator(1e-4, ttm, uN) 

	f1_vec_future = @spawn [g1_int(vt, vgrid, u, 
				       _lambda, sigmal,
				       r, gross_delta, 
				       xi, kappa) for vt=vtgrid, u=ugrid]
	f1_vec = sum(fetch(f1_vec_future), 2) * du
	
	f1_itp = interpolate(f1_vec[1:end], BSpline(Cubic(Line())), OnGrid())
	return Interpolations.scale(f1_itp, vtgrid)

end


# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 

function f2_fd(vtgrid, vbhlgrid, 
		     _lambda, sigmal, sigmah,
		     r, gross_delta, xi, kappa;
		     vmax=1.2, vN=1e3,ttm=1.0, uN=1e3)

	du, ugrid = grid_creator(1e-4, ttm, uN)
	dv, vgrid = grid_creator(1e-4, vmax, vN)

	f2_surf = @spawn [f2_int(vt, vbhl, ugrid,  vgrid,  
				       _lambda, sigmal, sigmah,
				       r, gross_delta, 
				       xi, kappa) for vt=vtgrid, vbhl=vbhlgrid]
	
	f2_itp = interpolate(fetch(f2_surf), BSpline(Cubic(Line())), OnGrid())
	return Interpolations.scale(f2_itp, vtgrid, vbhlgrid)

end


# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 

function f3_fd(vtgrid, vbhlgrid, 
		     _lambda, sigmal, sigmah,
		     r, gross_delta, xi, kappa;
		     vmax=1.2, vN=1e3,ttm=1.0, uN=1e3)

	du, ugrid = grid_creator(1e-4, ttm, uN)
	dv, vgrid = grid_creator(1e-4, vmax, vN)

	f3_surf = @spawn [f3_int(vt, vbhl, ugrid,  vgrid,  
			      _lambda, sigmal, sigmah,
			      r, gross_delta, 
			      xi, kappa) for vt=vtgrid, vbhl=vbhlgrid]
	
	f3_itp = interpolate(fetch(f3_surf), BSpline(Cubic(Line())), OnGrid())
	return Interpolations.scale(f3_itp, vtgrid, vbhlgrid)

end

# ##################################################################
# ########################### Interp ###############################
# ################################################################## 

function fd_interp(svm;
        vtmax=.8, vtN=20, vtgrid=nothing,
        vbhlmin=.75, vbhlmax=1.25, vbhlN=5, vbhlgrid=nothing,
        vmax=1.2, vN=10^3, uN=10^3)
	
        # ####################################
        # ######## Extract Parameters ########
        # ####################################
        # Capital Structure
        m = svm.pm.m
        c = svm.c
        p = svm.p
    
        # Taxes & Default
        
        # Dynamics
        r = svm.pm.r
        gross_delta = svm.pm.gross_delta
        # iota?
    
        # Liquidity
        xi = svm.pm.xi
        kappa = svm.pm.kappa
        
        # Volatility
        _lambda = svm.pm.lambda
        sigmal = svm.pm.sigmal
        sigmah = svm.pm.sigmah 
        # ####################################
        
        # ####################################
        # ##### Form Interpolation Grids #####
        # ####################################
        if vtgrid==nothing
           _, vtgrid=grid_creator(0.0, vtmax, vtN)
        end
        
        if vbhlgrid==nothing
           _, vbhlgrid=grid_creator(vbhlmin, vbhlmax, vbhlN)
        end
        # ####################################

	f0 = @spawn f0_fd(vtgrid, m, vmax, 
                          _lambda, sigmal, 
                          c, p, r, 
                          gross_delta, xi, kappa, uN=uN)
        
        f1 = @spawn f1_fd(vtgrid, _lambda, sigmal, 
                          r, gross_delta, xi, kappa;
                          vmax=vmax, vN=uN,
                          ttm=m, uN=uN)
    
        f2 = @spawn f2_fd(vtgrid, vbhlgrid,
                          _lambda, sigmal, sigmah,
                          r, gross_delta, xi, kappa;
                          vmax=vmax, vN=uN,
                          ttm=m, uN=uN)
        
        f3 = @spawn f3_fd(vtgrid, vbhlgrid,
                          _lambda, sigmal, sigmah,
                          r, gross_delta, xi, kappa;
                          vmax=vmax, vN=uN,
                          ttm=m, uN=uN)
        
        return fetch(f0), fetch(f1), fetch(f2), fetch(f3)
end
