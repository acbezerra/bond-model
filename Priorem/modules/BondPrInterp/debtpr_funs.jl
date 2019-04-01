
function get_pv_rfdebt(svm)
   return rf_debt(svm.m, svm.c, svm.p, 
                  svm.pm.r, svm.pm.xi, svm.pm.kappa)
end
# function get_pv_rfdebt(svm)
# 	
# 	m = svm.m
#  	c = svm.c
# 	p = svm.p
# 	rdisc = rdisc(svm.pm.r, svm.pm.xi, svm.pm.kappa) 
# 
#         return ((c/ rdisc) * (m- (1 - exp(- rdisc * m)) / rdisc) +
# 		 p * (1 -  exp(- rdisc * m)) / rdisc) 
# end


function get_cvm_debt_price(svm, sigma; c=nothing, p=nothing,
                            Vt=nothing, ttmMax=nothing, 
			    vmax=nothing, N1=50, N2=10^4)
	# ###############################
	# Capital Structure	
	m = svm.m
    
        if c==nothing
            c = svm.c
        end
        if p==nothing
            p = svm.p
        end
    
 	# Dynamics
	r = svm.pm.r	
	gross_delta = svm.pm.gross_delta

	# Default and Taxes
	alpha = svm.pm.alpha
	pii = svm.pm.pi

	# Liquidity
	xi = svm.pm.xi
	kappa = svm.pm.kappa
	# ###############################

	if Vt==nothing
		Vt = svm.pm.V0
	end

	if ttmMax == nothing
		ttmMax = m
	end

	vb = get_cvm_vb(svm, sigma)
	v = log(Vt/vb)
	
	if vmax==nothing
		Vmax = get_bond_Vmax(svm)
		vmax = log(Vmax/float(vb))
	end

	
	# Create time-to-maturity grids
	_, ttm_grid1 = grid_creator(0.0, ttmMax, N1)
	dt2, ttm_grid2 = grid_creator(0.0, ttmMax, N2)


	# Compute Bond Prices
	bond_vec = @spawn [cvm_bond_price(vb, v, ttm, m, c, p, 
					  sigma, r, gross_delta, 
					  xi, kappa, alpha, pii)
			   for ttm=ttm_grid1]


	# Interpolate	
    
	bond_interp_itp = interpolate(fetch(bond_vec), BSpline(Cubic(Line(OnGrid()))))
	bond_interp_sitp =  Interpolations.scale(bond_interp_itp, ttm_grid1)

	# Refine Bond Price Grid
	bond_vec2 = @spawn [bond_interp_sitp(ttm) for ttm=ttm_grid2]

	# Integrate
	return sum(fetch(bond_vec2)) * dt2
end


function get_svm_debt_price(svm, vbl, f11, f12, f13, f21, f22; 
			    c=nothing, p=nothing,
                            Vt=nothing, ttmMax=nothing, vmax=nothing, N=10^4) 
		
        if c==nothing
            c = svm.c
        end
        if p==nothing
            p = svm.p
        end
 
	# Here ttmMax, Vt and vbl should default to the values
	# used to compute the f functions!
	# vmax can be either the value used to compute the f functions 
	# or the value given by get_bond_vmax function, if smaller.
	
	if Vt==nothing
		Vt = svm.pm.V0
	end

	if ttmMax == nothing
		ttmMax = svm.m
	end

	if vmax==nothing
		Vmax = get_bond_Vmax(svm)
		vmax = log(Vmax/float(vbl))
	end

	
	# Create  time-to-maturity grid
	dt, ttm_grid = grid_creator(0.0, ttmMax, N)

	# Get Bond Prices at Different Maturities
	bondpr_vec = @spawn [bondpr(svm, Vt, vbl, ttm, vmax, f11, f12, f13,
				    f21, f22, c=c, p=p, fd_method=false) for
			     ttm=ttm_grid]

	return sum(fetch(bondpr_vec)) * dt
end

