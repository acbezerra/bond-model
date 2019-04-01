
# function deriv_calculator(f, x, h)
#     return (f(x+h) - f(x))/h
# end

function eq_fd_newly_issued_bonds(svm, vbl, mu_b, c, p, vgrid;
                                  vtN::Int64=1000, ftype::String="bf")

    rfbond = rfbond_price(svm.m, c, p, svm.pm.r, svm.pm.xi, svm.pm.kappa)
    dpayoff = on_default_payoff(0., vbl, svm.m,
                                mu_b, svm.m, c, p,
                                svm.pm.r, svm.pm.xi,
                                svm.pm.kappa, svm.pm.alpha)

    _, v_subgrid = grid_creator((1 + 1e-4) * minimum(svm.bs.vtgrid), maximum(vgrid), vtN)

    bpr_vec = fetch(@spawn [get_svm_bond_price(svm, vbl, svm.m;
                                               Vt=vbl*exp(v), mu_b=mu_b,
                                               c=c, p=p, ftype=ftype) for v in v_subgrid])

    # bpr = sp_interpolate.interp1d(vcat(.0, v_subgrid), vcat(dpayoff, bpr_vec), 
    #                              kind="cubic", fill_value="extrapolate")
    bpr = Dierckx.Spline1D(vcat(.0, v_subgrid), vcat(dpayoff, bpr_vec); k=3, bc="extrapolate")

    return Array([minimum([maximum([bpr(v)[1], dpayoff]), rfbond]) for v in vgrid])[2:end-1]
end


function eq_fd_core_cvmh_eq_values(svm, vbl, mu_b, c, p, vgrid)
    tic  = time()
    println("Computing Constant Volatility Equity Values")

    # vbh = get_cvm_vb(svm, svm.pm.sigmah; mu_b=mu_b, c=c, p=p)
    cvm_eqh_all = fetch(@spawn [get_cvm_eq(svm, vbl * exp(v), svm.pm.sigmah;
                                           mu_b=mu_b, c=c, p=p) for v in vgrid])
    
    println("Finished computing Constant Volatility Equity Values")
    println(string("Time to compute Constant Volatility Equity Values: ",  time() - tic))
    println(" ")

    return Array{Float64}(cvm_eqh_all)[2:end-1]
end


function eq_fd_core_coeffs(svm, vgrid)
    deltav = vgrid[1] - vgrid[2] 
    nu = get_rgrow(svm) - .5 * svm.pm.sigmal^2

    qu = .5 * (nu / deltav + svm.pm.sigmal^2 / (deltav^2))
    qd = .5 * (-nu / deltav + svm.pm.sigmal^2 / (deltav^2))
    qm = - (svm.pm.r + svm.pm.lambda + svm.pm.sigmal^2 / (deltav^2))

    return Dict{Symbol,Float64}(:deltav => deltav,
                                :nu => nu,
                                :qu => qu,
                                :qd => qd,
                                :qm => qm)
end


function eq_fd_core_matrices(svm, vbl, mu_b, c, p,
                             vgrid, eq_vbl, eq_max,
                             bond_prices, cvm_eqh,
                             coeffs)

    # Gamma Vector:
    Gamma = (get_param(svm, :delta) * vbl * exp.(vgrid[2:end-1]) .+ 
             mu_b .* (-(1 - get_param(svm, :pi)) .* (svm.m * c) .+ 
             bond_prices .- p) .+
             get_param(svm, :lambda) * cvm_eqh)

    
    println(string("Shape of Gamma matrix: ", size(Gamma)))
    Gamma[1] += coeffs[:qu] * eq_max    
    Gamma[end] += coeffs[:qd] * eq_vbl

    # A Matrix:
    A = (coeffs[:qm] * Array(Diagonal(ones(size(Gamma, 1)))) +
         coeffs[:qu] * [1. *(y==x-1) for x in 1:size(Gamma, 1), y in 1:size(Gamma, 1)] +
         coeffs[:qd] * [1. *(y==x+1) for x in 1:size(Gamma, 1), y in 1:size(Gamma, 1)])

    return Gamma, A
end


function eq_fd_core_eq_values(svm, vbl, vgrid, eq_vbl, eq_max, Gamma, A)
    # ###### Compute Pre-Volatility Shock Equity Function: ######
    # Form Function and add lim_v->infty E(v)
    eq_vals = vcat(eq_max, -\(A, Gamma), eq_vbl)

    # Interpolate to back-out equity value at VB:
    eq_spl = Dierckx.Spline1D([vbl * exp(v) for v in reverse(vgrid)], reverse(eq_vals); k=3, bc="extrapolate")

    # eq_spl = sp_interpolate.interp1d((vbl * exp.(reverse(vgrid))),
    # reverse(eq_vals), kind="cubic", fill_value="extrapolate")
    # eq_spl = sp_interpolate.PchipInterpolator((vbl * exp.(reverse(vgrid))),
    #                                   reverse(eq_vals))

    println(string("equity: ", eq_spl(get_param(svm, :V0))))
    
    # Compute Derivative at Default Barrier:
    eq_spl_deriv = fetch(@spawn [Dierckx.derivative(eq_spl, vbl * exp(v)) for v in vgrid])
    eq_deriv = eq_spl_deriv[end]

    # Compute Derivative at Default Barrier:
    # eq_spl_deriv = eq_spl.derivative()(vbl .* exp(reverse(vgrid)))
    # eq_deriv = eq_pchip_deriv[1]
    # println(string("eq_deriv: ", eq_deriv))
    
    # # #######################################
    # h = 1e-5
    # eq_spl_deriv = fetch(@spawn [deriv_calculator(eq_spl, vbl * exp(v), h) 
    #                               for v=reverse(vgrid)])
    # eq_deriv = eq_spl_deriv[1]
    # # #######################################

    # Equity Values
    # Equity set as function of V:
    e0 = eq_spl(get_param(svm, :V0))
    # e0 = eq_sitp(log(get_param(svm, :V0)/float(vbl)))

    # Look for negative values
    eq_min_val = minimum(reverse(eq_vals)[2:end])  # interpolate in a neighborhood
    eq_negative = eq_min_val .< -0.05 
    eq_deriv_min_val = minimum(eq_spl_deriv) 

    eq_dict = Dict(:e0 => e0,
                   :eq_max => eq_max,
                   :eq_vals=>  eq_vals,
                   :eq_deriv => eq_deriv,
                   :eq_spl_deriv => eq_spl_deriv, 
                   :eq_vb => eq_vbl,
                   :eq_min_val => eq_min_val,
                   :eq_negative => eq_negative,
                   :eq_deriv_min_val => eq_deriv_min_val)
end


function eq_fd_core(svm, vbl, mu_b, c, p, vgrid, bond_prices)
    core_tic = time()

    # #################################
    # ######## Boundary Values ########
    # #################################
    # Upper Barrier: Value of Equity
    eq_max = get_cvm_eq(svm, vbl * exp(vgrid[1]),
                        svm.pm.sigmal; mu_b=mu_b, c=c, p=p)
    println(string("eq_max: ", eq_max))
    
    # Lower Barrier:
    eq_vbl = maximum([0., get_param(svm, :alpha) * vbl - get_pv_rfdebt(svm, mu_b=mu_b, c=c, p=p)])
    println(string("eq_vbl: ", eq_vbl))
    
    # ##################################
    # ######## CVM Equity Values #######
    # ##################################
    cvm_eqh = eq_fd_core_cvmh_eq_values(svm, vbl, mu_b, c, p, vgrid)

    # #################################
    # ######### Coefficients: #########
    # #################################
    coeffs = eq_fd_core_coeffs(svm, vgrid)

    # #################################
    # ########### Matrices: ###########
    # #################################
    Gamma, A = eq_fd_core_matrices(svm, vbl, mu_b, c, p,
                                   vgrid, eq_vbl, eq_max,
                                   bond_prices, cvm_eqh,
                                   coeffs)  

    # #################################
    # ######### Equity Values #########
    # #################################
    println("Computing equity values... ")
    core_eq_tic = time()
    eq_dict = eq_fd_core_eq_values(svm, vbl, vgrid,
                                   eq_vbl, eq_max,
                                   Gamma, A)
    
    println(string("Equity Core Function Computation Time: ", time() - core_eq_tic))
    
    println(string("Total Equity FD Core Function Computation Time: ", time() - core_tic))

    # return Gamma, A, eq_dict
    return eq_dict
end


function eq_fd_export_results(svm, vbl::Float64,
                              mu_b::Float64,
                              c::Float64,
                              p::Float64,
                              eq_dict,
                              debt::Float64=NaN)
    if isnan(debt)
        debt = get_svm_debt_price(svm, vbl; mu_b=mu_b, c=c, p=p)
    end
        
    # Get Firm Value:
    firm_value = debt + eq_dict[:e0]

    # Compute Leverage:
    lev = (debt / firm_value) * 100

    # Compute ROE:
    roe = (eq_dict[:e0] / (get_param(svm, :V0) - debt) - 1.) * 100

    results = Dict(:V0 =>  get_param(svm, :V0),
                   :r =>  svm.pm.r,
                   :gross_delta =>  get_param(svm, :gross_delta),
                   :iota=> svm.pm.iota,
 		   :delta=> get_param(svm, :delta),
                   :alpha=> svm.pm.alpha,
                   :pi=> svm.pm.pi,
                   :xi=> svm.pm.xi,
                   :kappa=> svm.pm.kappa,
                   :lambda=> svm.pm.lambda,
                   :sigmal => svm.pm.sigmal,
                   :sigmah => svm.pm.sigmah,
                   :mu_b => mu_b,
                   :m =>  svm.m,
                   :c =>  c,
                   :p =>  p,
                   :vb =>  vbl,
                   :debt =>  debt,
                   :equity =>  eq_dict[:e0],
                   :eq_deriv =>  eq_dict[:eq_deriv],
                   :firm_value =>  firm_value,
                   :eq_min_val =>  eq_dict[:eq_min_val],
                   :eq_vb =>  eq_dict[:eq_vb],
                   :eq_negative =>  eq_dict[:eq_negative],
                   :eq_deriv_min_val =>  eq_dict[:eq_deriv_min_val],
                   :leverage =>  lev,
                   :ROE =>  roe)

    df = DataFrame(results)

    return results, df
end

 
function eq_fd(svm, vbl::Float64;
    mu_b::Float64=NaN,
    c::Float64=NaN,
    p::Float64=NaN,
    debt::Float64=NaN,
    ftype::String="bf")

    tic = time()

    if isnan(mu_b)
        mu_b=svm.mu_b
    end
    
    if isnan(c)
        c = svm.c
    end
    
    if isnan(p)
        p = svm.p
    end

    # V MAX:
    println("Computing Equity Vmax")
    eq_Vmax = get_eq_Vmax(svm; mu_b=mu_b, c=c, p=p)
    println(string("Equity Vmax: ", eq_Vmax))
    println(" ")

    vtN = 1500
#     vtgrid = reverse(grid_creator(0.0, log(eq_Vmax/float(vbl)), vtN)[2])
    vtgrid = reverse(range(0.0, stop=log(eq_Vmax/float(vbl)), length=vtN))

    # #################################
    # ###### Newly-Issued Bonds #######
    # #################################
    # Newly-Issued Bond Prices
    bond_prices = eq_fd_newly_issued_bonds(svm, vbl, mu_b, c, p, vtgrid; ftype=ftype)
    
    eq_dict = eq_fd_core(svm, vbl, mu_b, c, p, vtgrid, bond_prices)
    _, df = eq_fd_export_results(svm, vbl, mu_b, c, p, eq_dict)

    println(string("Total computation time: ", time() - tic))

    return df
end





# function eq_fin_diff_core(svm, vbl, v_grid, bondpr)
#         # #################################
#         # ######### Equity Values #########
#         # #################################
#         # Upper Barrier: Value of Equity
#         eq_max = get_cvm_eq(vbl * exp(v_grid[1]), svm.pm.sigmal, svm)

#         # Note: I can use the CVM value here because both 
#         # the cvm and svm equity values will have converged 
#         # to the credit-risk-free equity value in the upper
#         # barrier.
        
#         # ##################################
#         # ##### Baseline Model Equity Values #####
#         # ##################################
#         println("Computing Constant Volatility Equity Values")
#         bsm_eqh_vals_Future = @spawn [get_cvm_eq(vbl * exp(v), svm.pm.sigmah, svm) for v=v_grid]
#         bsm_eqh_vals = fetch(bsm_eqh_vals_Future)
    
#         # Interpolate
# #         itp = interpolate(reverse(bsm_eqh_sub_vals), BSpline(Cubic(Line())), OnGrid())
# #         bsm_interp_eqh = Interpolations.scale(itp, reverse(v_sub_grid))  # Scale

#         # Equity Values
# #         bsm_eqh = bsm_interp_eqh(v_grid[2:end-1])
#         bsm_eqh = bsm_eqh_vals[2:end-1]
#         println("Finished computing Constant Volatility Equity Values")
    
#         # #################################
#         # ######### Coefficients: #########
#         # #################################
#         deltav = v_grid[1] - v_grid[2]
#         nu = get_rgrow(svm) - .5 * get_param(svm, :sigmal)^2

#         qu = .5 * (nu / deltav + get_param(svm, :sigmal)^2 / (deltav^2))
#         qd = .5 * (-nu / deltav + get_param(svm, :sigmal)^2 / (deltav^2))
#         qm = -(get_param(svm, :r) + get_param(svm, :lambda) + 
#                                     get_param(svm, :sigmal)^2 / (deltav^2))

#         println("deltav: ", string(deltav))
#         println("nu: ", string(nu))
#         println("qu: ", string(qu))
#         println("qm: ", string(qm))
#         println("qd: ", string(qd))
    
#         # Present Value of Debt:
#         pv_debt = get_pv_rfdebt(svm)
    
#         println("Finished computing pv rfdebt")
    
#         # Gamma Vector:
#         Gamma = (get_param(svm, :delta) * vbl * exp.(v_grid[2:end-1]) .- 
#                     (1 - get_param(svm, :pi)) .* get_param(svm, :c) .+ 
#                 bondpr[2:end-1] .- get_param(svm, :p) .+
#                     get_param(svm, :lambda) * bsm_eqh)

#         println("Shape of Gamma matrix: ", string(size(Gamma)))
#         Gamma[1] += qu * eq_max

#         eq_vbl = max(0., get_param(svm, :alpha) * vbl - pv_debt)
#         println("eq_vbl: ", string(eq_vbl))
#         Gamma[end] += qd * eq_vbl

#         # A Matrix:
#         A = (qm * Array(Diagonal(ones(length(Gamma)))) +
#              qu * [1. *(y==x-1) for x in 1:length(Gamma), y in 1:length(Gamma)] +
#              qd * [1. *(y==x+1) for x in 1:length(Gamma), y in 1:length(Gamma)])
        
#         # ###### Compute Pre-Volatility Shock Equity Function: ######
#         # Form Function and add lim_v->infty E(v)
#         eq_vals = vcat(eq_max, -\(A, Gamma), eq_vbl)
    
#         # Interpolate to back-out equity value at VB:
#         eq_interp_itp = interpolate(reverse(eq_vals), BSpline(Cubic(Line())), OnGrid())
#         eq_sitp = Interpolations.scale(eq_interp_itp, reverse(v_grid))  # Scale

#         # Compute Derivative at Default Barrier:
#         h = 1e-5
#         inv_v_grid = reverse(v_grid)
#         eq_sitp_deriv_fut = @spawn [deriv(eq_sitp, v, h) * (1/(vbl * exp(v))) 
#                                         for v=inv_v_grid]
#         eq_sitp_deriv = fetch(eq_sitp_deriv_fut)
#         eq_deriv_vb = eq_sitp_deriv[1]

#         # Equity Values

#         # Equity set as function of V:
#         e0 = float(eq_sitp[log(get_param(svm, :V0)/vbl)])
    
#         # Look for negative values
#         eq_min_val = minimum(reverse(eq_vals)[2:end])  # interpolate in a neighborhood
#         eq_negative = eq_min_val .< -0.05 
#         eq_deriv_min_val = minimum(eq_sitp_deriv) 

#         eq_dict = Dict(:e0 => e0,
#                        :eq_max => eq_max,
#                        :eq_vals=>  eq_vals,
#                        "eq_sitp" =>  eq_sitp,
#                        :eq_deriv => eq_deriv_vb,
#                        "eq_sitp_deriv" => eq_sitp_deriv,
#                        :eq_vb => eq_vbl,
#                        :eq_min_val => eq_min_val,
#                        :eq_negative => eq_negative,
#                        :eq_deriv_min_val => eq_deriv_min_val, 
#                        "invA_Gamma" => - inv(A) * Gamma)
                                                                
#         return Gamma, A, eq_dict
# end


# function eq_export_fin_diff_res(svm, vbl, eq_dict; debt=nothing, bprfs=nothing)
#      # Note: bprfs is a container for the interpolated functions 
#      #       used in the computation of bond prices at different
#      #       maturities.
#      if bprfs!=nothing
# 	debt = get_svm_debt_price(svm, vbl, bprfs["f0mat"],
#                                             bprfs["f1mat"],
# 	                                    bprfs["f2mat"],
# 	                                    bprfs["f3mat"])
#      elseif debt==nothing
#          println("WARNING: Assuming debt is issued at par")
#          println("   ===>  Setting debt = P")
#          debt = get_param(svm, :p)
#      end

#     # Get Firm Value:
#     firm_value = debt + eq_dict[:e0]

#     # Compute Leverage:
#     lev = (debt / firm_value) * 100

#     # Compute ROE:
#     roe = (eq_dict[:e0] / (svm.pm.V0 - debt) - 1) * 100

#     results = Dict(:V0=> svm.pm.V0,
#                    :r=> svm.pm.r,
#                    :gross_delta=> svm.pm.gross_delta,
#                    :iota=> svm.pm.iota,
# 		   :delta=> get_param(svm, :delta),
#                    :alpha=> svm.pm.alpha,
#                    :pi=> svm.pm.pi,
#                    :xi=> svm.pm.xi,
#                    :kappa=> svm.pm.kappa,
#                    :lambda=> svm.pm.lambda,
#                    :sigmal=> svm.pm.sigmal,
#                    :sigmah=> svm.pm.sigmah,
#                    :m=> svm.m,
#                    :c=> svm.c,
#                    :p=> svm.p,
#                    "vb"=> vbl,
# 		   "debt"=> debt,
#                    "equity"=> eq_dict[:e0],
#                    :eq_deriv=> eq_dict[:eq_deriv],
#                    "firm_value"=> firm_value,
#                    :eq_min_val=> eq_dict[:eq_min_val],
#                    :eq_vb=> eq_dict[:eq_vb],
#                    :eq_negative=> eq_dict[:eq_negative],
#                    :eq_deriv_min_val=> eq_dict[:eq_deriv_min_val],
#                    "leverage"=> lev,
#                    "ROE"=> roe)

#     df = DataFrame(results)

#     return results, df
# end


# function equity_fin_diff(svm, vbl, f11fd, f12fd, f13fd, f21fd, f22fd;
#                          p=nothing, debt=nothing, bprfs=nothing, 
# 			 full_output=false)
	
# 	if p==nothing
# 	    p = svm.p
# 	end
	
# 	# V MAX:
# 	println("Computing Equity Vmax")
# 	eqVmax = get_eq_Vmax(svm)
# 	println("Equity Vmax: ", string(eqVmax))
# 	println(" ")
	
# 	println("Computing Bond Vmax")
# 	bondVmax = get_bond_Vmax(svm)
# 	println("Bond Vmax: ", string(bondVmax))
# 	println(" ")
	
# 	# Ensure eq_vmax > bond_vmax
# 	eqVmax = max(1.25 * bondVmax, eqVmax)
# 	println("Adjusted Equity Vmax: ", string(eqVmax))
	
# 	# v max
# 	eq_vmax = log(eqVmax / vbl)
# 	bond_vmax = log(bondVmax / vbl)
	
# 	# Ensure eq_vmax > bond_vmax
# 	if bond_vmax >= eq_vmax
# 	    eq_vmax = 1.1 * bond_vmax
# 	end
# 	println("eq_vmax: ", string(eq_vmax))
# 	# #######################################
# 	# ####### Form grid of v values #########
# 	# #######################################
# 	lenv = 1500
# 	v_grid = range(eq_vmax, stop=0.0, length=lenv)
	
# 	bondpr_grid_fut = @spawn [bondpr(svm, vbl * exp(v), vbl,
# 	                                 svm.m, bond_vmax,
#                                          f11fd, f12fd, f13fd,
#                                          f21fd, f22fd) for v=v_grid]
# 	bondpr_grid = fetch(bondpr_grid_fut)
	
# 	# refine grids
# 	#  v_grid_ref = range(v_grid[1], stop=v_grid[end], length=10^4)
# 	#  bondpr_itp = interpolate(reverse(bondpr_grid), BSpline(Cubic(Line())), OnGrid())
# 	#  bondpr_sitp = Interpolations.scale(bondpr_itp, reverse(v_grid))  # Scale
# 	#  bondpr_ref = fetch(@spawn [bondpr_sitp[v] for v=v_grid_ref])
	
# 	# #################################
# 	# ######### Equity Values #########
# 	# #################################
# 	println("Finished computing bondpr_grid")
# 	Gamma, A, eq_dict = eq_fin_diff_core(svm, vbl, v_grid, bondpr_grid)
	
# 	println("Exporting results... ")
# 	_, df = eq_export_fin_diff_res(svm, vbl, eq_dict, debt=debt, bprfs=bprfs)
	
# 	return v_grid, eq_dict, Gamma, A, df
# end

