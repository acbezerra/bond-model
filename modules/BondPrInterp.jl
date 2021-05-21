module_path = "/home/artur/BondPricing/bond-model/modules/"
push!(LOAD_PATH, module_path)
# modnames = ["ModelObj", "AnalyticFunctions"]
# for modl in modnames
#     if !(joinpath(module_path, modl) in LOAD_PATH)
#         push!(LOAD_PATH, joinpath(module_path, modl))
#     end
# end


# Interpolate Primitive SVM Bond Pricing Functions
module BondPrInterp

using Distributed  # @spawn macro
using Distributions
using Interpolations
# using PyCall
# @pyimport scipy.interpolate as sp_interpolate
# @pyimport numpy as np
# using LinearAlgebra
using Dierckx

# User-defined packages:
using ModelObj: grid_creator, set_bpr_grids
using AnalyticFunctions: rdisc, rdisc_pvs, rfbond_price,
                         cvm_bond_price, rf_debt,
                         get_rdisc, get_rgrow, get_k_struct,
                         get_cvm_vb, get_param,
                         psi_v_td, dv_psi_v_td,
                         cvm_F, cvm_G, cvm_vb, get_cvm_vb,
                         on_default_payoff, no_vol_shock_cf_pv


# * SVM Primitive Bond Pricing Functions #####################
# include("svm_primitive_bpr_funs.jl")
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


# ** J2: f21 and f22
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


# ** f1_int
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

# ** f2_int
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

# ** f3_int
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


# * SVM Interpolate Bond Pricing Functions ###################
# include("svm_interp_bpr_funs.jl")
# ** J2: f21, f22
# ##################################################################
# ######################### J2: f21, f22 ###########################
# ################################################################## 

# function f21_inputs(svm, vtgrid, ttmgrid; vmax=1.2, uN=10^3)
function f21_inputs(svm)   
    f21_fd_input = @spawn [f21_int(vt,
                                   svm.bi.vmax,
                                   svm.pm.lambda,
                                   svm.pm.sigmal,
                                   svm.pm.r,
                                   svm.pm.gross_delta,
                                   svm.pm.xi,
                                   svm.pm.kappa;
                                   ttm=tau, N=svm.bi.uN) for vt in svm.bs.vtgrid,
                                                            tau in svm.bs.ttmgrid]

    return fetch(f21_fd_input)
end


# function f22_inputs(svm, vtgrid, ttmgrid; vmax=1.2, uN=10^3)
function f22_inputs(svm)
    f22_fd_input = @spawn [f22_int(vt,
                                   svm.bi.vmax,
                                   svm.pm.lambda,
                                   svm.pm.sigmal,
                                   svm.pm.r,
                                   svm.pm.gross_delta,
                                   svm.pm.xi,
                                   svm.pm.kappa; 
				   ttm=tau, N=svm.bi.uN) for vt in svm.bs.vtgrid,
                                                            tau in svm.bs.ttmgrid]

    return fetch(f22_fd_input)
end


# ** f11_int
# ##################################################################
# ########################### f11_int ##############################
# ################################################################## 
# function f11_inputs(svm, vtgrid, ttmgrid;
#                        vmax=1.2, vN=10^3, uN=10^3)
function f11_inputs(svm)
    dv, vgrid = grid_creator(1e-4, svm.bi.vmax, svm.bi.vN)

    f11_vec = @spawn [f11v_int(vt,
                               grid_creator(0.0, ttm, svm.bi.uN)[2],
                               vgrid,
                               svm.pm.lambda,
                               svm.pm.sigmal,
                               svm.pm.r,
                               svm.pm.gross_delta,
                               svm.pm.xi,
                               svm.pm.kappa) for vt in svm.bs.vtgrid,
                                                ttm in svm.bs.ttmgrid]

    return fetch(f11_vec)
end


# ** f2_int
# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 
# function f12_inputs(svm, vtgrid, ttmgrid, vbhlgrid;
# 		       vmax=1.2, vN=10^3, uN=10^3)
function f12_inputs(svm)
    _, vgrid = grid_creator(1e-4, svm.bi.vmax, svm.bi.vN)

    f12_surf = @spawn [f12_int(vt,
                               vbhl,
                               grid_creator(0.0, ttm, svm.bi.uN)[2],
                               vgrid,  
			       svm.pm.lambda,
                               svm.pm.sigmal,
                               svm.pm.sigmah,
                               svm.pm.r,
                               svm.pm.gross_delta, 
			       svm.pm.xi,
                               svm.pm.kappa) for vt in svm.bs.vtgrid,
                                                ttm in svm.bs.ttmgrid,
                                               vbhl in svm.bs.vbhlgrid]

    return fetch(f12_surf)
end


# ** f3_int 
# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 
# function f13_inputs(svm, vtgrid, ttmgrid, vbhlgrid;
# 		       vmax=1.2, vN=10^3, uN=10^3)
function f13_inputs(svm)
    _, vgrid = grid_creator(1e-4, svm.bi.vmax, svm.bi.vN)

    f13_surf = @spawn [f13_int(vt, vbhl,
                               grid_creator(0.0, ttm, svm.bi.uN)[2],
                               vgrid,  
			       svm.pm.lambda,
                               svm.pm.sigmal,
                               svm.pm.sigmah,
                               svm.pm.r,
                               svm.pm.gross_delta, 
			       svm.pm.xi,
                               svm.pm.kappa) for vt in svm.bs.vtgrid,
                                                ttm in svm.bs.ttmgrid,
                                               vbhl in svm.bs.vbhlgrid]

    return fetch(f13_surf)
end	

# ** Collect Surfaces to Form Bond Pricing Function
# ##################################################################
# ########################### Interp ###############################
# ################################################################## 
function bpr_surfs(svm)
    # ####################################
    # ##### Form Interpolation Grids #####
    # ####################################
    svm = set_bpr_grids(svm)
    # ####################################
    
    svm.bs.f11_surf = fetch(@spawn f11_inputs(svm))
    svm.bs.f12_surf = fetch(@spawn f12_inputs(svm))
    svm.bs.f13_surf = fetch(@spawn f13_inputs(svm))
    svm.bs.f21_surf = fetch(@spawn f21_inputs(svm))
    svm.bs.f22_surf = fetch(@spawn f22_inputs(svm))

    return svm
end


function gen_ref_surfs(svm)
    return Dict("ttmgrid_ref" => range(minimum(svm.bs.ttmgrid),
                                       stop=svm.bi.ttm_max,
                                       length=svm.bi.ttmN_ref),
                "vtgrid_ref" => range(minimum(svm.bs.vtgrid),
                                      stop=svm.bi.vtmax,
                                      length=svm.bi.vtN_ref),
                "vbhlgrid_ref" => range(minimum(svm.bs.vbhlgrid),
                                        stop=svm.bi.vbhlmax,
                                        length=svm.bi.vbhlN_ref))
end


function bpr_interp_f1f(svm, mat, ref_surfs)
    xN = size(svm.bs.vtgrid, 1)
    # spl = Dierckx.Spline2D(svm.bs.vtgrid, svm.bs.ttmgrid, mat;
    #                        kx=3, ky=3, s=0.0)
    spl = Dierckx.Spline2D(svm.bs.vtgrid, svm.bs.ttmgrid, mat;
                           kx=3, ky=3, s=xN)


    tmp_surf = fetch(@spawn [Dierckx.evaluate(spl, vt, ttm) for
                             vt in ref_surfs["vtgrid_ref"],
                             ttm in ref_surfs["ttmgrid_ref"]])
    f_itp = Interpolations.interpolate(tmp_surf,
                                       BSpline(Cubic(Line(Interpolations.OnGrid()))))

    return Interpolations.scale(f_itp, ref_surfs["vtgrid_ref"],
                                ref_surfs["ttmgrid_ref"])
end


function bpr_interp_f2f(svm, mat)
    f_itp = Interpolations.interpolate(mat, BSpline(Cubic(Line(Interpolations.OnGrid()))))
    return Interpolations.scale(f_itp, svm.bs.vtgrid,
                                svm.bs.ttmgrid, svm.bs.vbhlgrid)
end


function bpr_interp(svm)
    ref_surfs = gen_ref_surfs(svm)

    svm.bf.f11 = bpr_interp_f1f(svm, svm.bs.f11_surf, ref_surfs)
    svm.bf.f12 = bpr_interp_f2f(svm, svm.bs.f12_surf)
    svm.bf.f13 = bpr_interp_f2f(svm, svm.bs.f13_surf)
    svm.bf.f21 = bpr_interp_f1f(svm, svm.bs.f21_surf, ref_surfs)
    svm.bf.f22 = bpr_interp_f1f(svm, svm.bs.f22_surf, ref_surfs)

    return svm
end


# * SVM Bond Pricing Fixed TTM Functions #####################
# include("svm_bpr_fixed_ttm_funs.jl")
function interp_bpr_fixed_tau_surfs(svm, bse)
    xN = size(svm.bs.vtgrid, 1)
    svm.bft.f11 = Dierckx.Spline1D(Array(svm.bs.vtgrid), 
                                   bse[:f11_surf]; 
                                   k=3, bc="extrapolate")
    svm.bft.f12 = Dierckx.Spline2D(Array(svm.bs.vtgrid), 
                                   Array(svm.bs.vbhlgrid), 
                                   bse[:f12_surf]; 
                                   kx=3, ky=3, s=xN)
    svm.bft.f13 = Dierckx.Spline2D(Array(svm.bs.vtgrid), 
                                   Array(svm.bs.vbhlgrid), 
                                   bse[:f13_surf]; 
                                   kx=3, ky=3, s=xN)
    svm.bft.f21 = Dierckx.Spline1D(Array(svm.bs.vtgrid), 
                                   bse[:f21_surf]; 
                                   k=3, bc="extrapolate")
    svm.bft.f22 = Dierckx.Spline1D(Array(svm.bs.vtgrid), 
                                   bse[:f22_surf]; 
                                   k=3, bc="extrapolate")
    
    return svm
end


function bpr_interp_fixed_ttm(svm; 
                              ttm::Float64=NaN,
                              vtmax::Float64=NaN,
                              vtN::Int64=0,
                              vbhlmin::Float64=NaN,
                              vbhlmax::Float64=NaN,
                              vbhlN::Int64=0)
    # Adjust grids ################################
    svm_grids = Dict()

    # vtgrid
    svm_grids[:vtgrid] = svm.bs.vtgrid
    if any([!isnan(vtmax), vtN > 5])
        isnan(vtmax) ? vtmax = svm.bi.vtmax  : svm.bit.vtmax=vtmax
        vtN < 5 ? vtN = svm.bi.vtN  : svm.bit.vtN = vtN
     
        svm.bs.vtgrid = range(0.0, stop=vtmax, length=vtN)
    end

    # ttmgrid
    svm_grids[:ttmgrid] = svm.bs.ttmgrid
    isnan(ttm) ? svm.bit.ttm = svm.m : svm.bit.ttm = ttm 
    svm.bs.ttmgrid = range(svm.bit.ttm, stop=svm.bit.ttm, length=1)

    # vbhlgrid  
    svm_grids[:vbhlgrid] = svm.bs.vbhlgrid
    if any([!isnan(vbhlmin), !isnan(vbhlmax), vbhlN > 5])
        isnan(vbhlmin) ? vbhlmin = svm.bi.vbhlmin  : svm.bit.vbhlmin = vbhlmin
        isnan(vbhlmax) ? vbhlmax = svm.bi.vbhlmax  : svm.bit.vbhlmax = vbhlmax
        vbhlN < 5 ? vbhlN = svm.bi.vbhlN  : svm.bit.vbhlN = vbhlN
        
        svm.bs.vbhlgrid = range(vbhlmin, stop=vbhlmax, length=vbhlN)
    end
    # #############################################

    # Form Surfaces ###############################
    f11_surf = fetch(@spawn BondPrInterp.f11_inputs(svm))[:,1]
    f12_surf = fetch(@spawn BondPrInterp.f12_inputs(svm))[:, 1, :]
    f13_surf = fetch(@spawn BondPrInterp.f13_inputs(svm))[:, 1, :]
    f21_surf = fetch(@spawn BondPrInterp.f21_inputs(svm))[:,1]
    f22_surf = fetch(@spawn BondPrInterp.f22_inputs(svm))[:,1]

    bst = Dict(:f11_surf => f11_surf,
               :f12_surf => f12_surf,
               :f13_surf => f13_surf,
               :f21_surf => f21_surf,
               :f22_surf => f22_surf)
    # #############################################
    
    
    # Interpolate Surfaces
    svm = interp_bpr_fixed_tau_surfs(svm, bst)
    
    # Restore Grids
    svm.bs.ttmgrid = svm_grids[:ttmgrid]
    svm.bs.vtgrid = svm_grids[:vtgrid]
    svm.bs.vbhlgrid = svm_grids[:vbhlgrid]
    
    return svm
end


# * Bond Pricing Functions ###################################
# include("bondpr_funs.jl")
function get_cvm_bond_price(svm, ttm::Float64,
                            sigma::Float64;
                            Vt::Float64=NaN,
                            vb::Float64=NaN,
                            mu_b::Float64=NaN,
                            m::Float64=NaN,
                            c::Float64=NaN,
                            p::Float64=NaN)
    # ####################################
    # ######## Extract Parameters ########
    # ####################################
    if isnan(Vt)
        Vt=svm.pm.V0
    end

    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)
    
    if isnan(vb)
        vb = get_cvm_vb(svm, sigma; mu_b=mu_b, m=m, c=c, p=p)
    end
    # ####################################

    return cvm_bond_price(vb, log(Vt/float(vb)), ttm,
                          mu_b, m, c, p, sigma,
                          svm.pm.r, svm.pm.gross_delta,
                          svm.pm.xi, svm.pm.kappa,
                          svm.pm.alpha, svm.pm.pi)

end


# ##########################################################
# #################### Fin Diff Method #####################
# ##########################################################
# 1. Notice the value of ttm passed to
# the function call must match the 
# value of ttm used to calibrate the 
# functions f0, f1, f2 and f3.
# 2. Notice also the order of the 
# arguments in the f2 and f3 functions:
# first vt, then vbhl.
# 3. Finally, vbhl in the f2 and f3 
# functions below is the ratio of the 
# post- to -pre-volatility shock
# bankruptcy barriers.

# ##########################################################
# ############# Prices at Different Maturities #############
# ##########################################################


# ##########################################################
# ##########################################################
function compute_jbpr_inputs(svm, Vt::Float64,
                             vbl::Float64,
                             mu_b::Float64,
                             c::Float64,
                             p::Float64)
    vt = log(Vt/float(vbl))
    
    # Default Barrier Ratio
    vbh = get_cvm_vb(svm, svm.pm.sigmah; mu_b=mu_b, c=c, p=p)
    vbhl = vbh/float(vbl)

    return vt, vbh, vbhl 
end


function get_interp_values(svm, vt::Float64,
                           ttm::Float64,
                           vbhl::Float64;
                           ftype::String="bf")
    if ftype == "bf"
        return svm.bf.f11(vt, ttm)[1],
               svm.bf.f12(vt, ttm, vbhl), 
               svm.bf.f13(vt, ttm, vbhl),
               svm.bf.f21(vt, ttm)[1],
               svm.bf.f22(vt, ttm)[1]
    elseif .&(ftype == "bft", abs.(ttm - svm.bit.ttm) < 1e-6)
        return svm.bft.f11(vt),
               svm.bft.f12(vt, vbhl), 
               svm.bft.f13(vt, vbhl),
               svm.bft.f21(vt),
               svm.bft.f22(vt)
    else
        println("Unable to compute bond pricing function values.")
        println(string("Function type: ", ftype))
        return 
    end
end


#= function get_svm_bond_price(svm, vbl::Float64, =#
#=                             ttm::Float64; =#
#=                             Vt::Float64=NaN, =#
#=                             mu_b::Float64=NaN, =#
#=                             c::Float64=NaN, =#
#=                             p::Float64=NaN, =#
#=                             ftype::String="bf") =#
#=     # #################################### =#
#=     # ######## Extract Parameters ######## =#
#=     # #################################### =#
#=     if isnan(Vt) =#
#=         Vt = svm.pm.V0 =#
#=     end =#

#=     mu_b, _, c, p = get_k_struct(svm; mu_b=mu_b, m=NaN, c=c, p=p) =#
#=     # #################################### =#
    
#=     # Maximum vt =#
#=     vtmax = svm.bi.vtmax =# 

#=     vt, vbh, vbhl = compute_jbpr_inputs(svm, Vt, vbl, mu_b, c, p) =# 
#=     # #################################### =#
#=     if vt <= 0 =#
#=         return on_default_payoff(vt, vbl, ttm, =#
#=                                  mu_b, svm.m, c, p, =#
#=                                  svm.pm.r, svm.pm.xi, =#
#=                                  svm.pm.kappa, svm.pm.alpha) =#
#=     elseif vt > minimum(svm.bs.vtgrid) =#
#=         rfbond = rfbond_price(ttm, c, p, =#
#=                               svm.pm.r, svm.pm.xi, svm.pm.kappa) =#

#=         if vt > vtmax =#
#=                return rfbond =# 
#=         end =#

#=         # Maturity or Default prior to Volatility Shock: =#
#=         cf0 = no_vol_shock_cf_pv(vt, vbl, ttm, =#
#=                                  mu_b, svm.m, =#
#=                                  c, p, =#
#=                                  svm.pm.sigmal, =#
#=                                  svm.pm.r, svm.pm.gross_delta, =#
#=                                  svm.pm.xi, svm.pm.kappa, =#
#=                                  svm.pm.alpha, svm.pm.lambda) =#
       
#=         f11, f12, f13, f21, f22 = get_interp_values(svm, vt, ttm, vbhl; ftype=ftype) =#
        
#=         # Volatility Shock Prior to Maturity: =#
#=         cf1 = c/get_rdisc(svm) * f11 + =#
#=               (p - c/get_rdisc(svm)) * f12 + =#
#=               (svm.pm.alpha * vbh/float(mu_b * svm.m) - c/get_rdisc(svm)) * f13 =#
       
#=         cf2 = c/get_rdisc(svm) * f21 + (p - c/get_rdisc(svm)) * f22 =# 
       
#=         return min(cf0 + cf1 + cf2, rfbond) =#
#=     else =#
#=         # In this case, firm is so close to default that I just =#
#=         # return the CVM bond price for sigma = sigmal: =#
#=         return get_cvm_bond_price(svm, ttm, svm.pm.sigmal; =#
#=                             Vt=Vt, vb=vbl, mu_b=mu_b, c=c, p=p) =#
#=     end =#
#= end =#


# Objective is to find V such that the value of the
# newly-issued (tau = m) risky bond price when sigma = sigmah
# is sufficiently close to the credit-risk-free bond price.
function get_bond_Vmax(svm; mu_b::Float64=NaN,
                       c::Float64=NaN,
                       p::Float64=NaN,
                       initV::Float64=NaN,
                       tol::Float64=.5 * 1e-3,
                       print::Bool=false)
    if isnan(mu_b)
        mu_b = svm.mu_b
    end
    
    if isnan(c)
        c=svm.c
    end

    if isnan(p)
        p=svm.p
    end

    if isnan(initV)
	initV = svm.pm.V0
    end

    bondVmax = 1.25 * initV
    vb = get_cvm_vb(svm, svm.pm.sigmah; mu_b=mu_b, c=c, p=p)
    rfbond = rfbond_price(svm.m, c, p,
                          svm.pm.r, svm.pm.xi, svm.pm.kappa)
    bondpr = get_cvm_bond_price(svm, svm.m, svm.pm.sigmah;
                                Vt=bondVmax, vb=vb, mu_b=mu_b, c=c, p=p)
    per_diff = (rfbond - bondpr) / rfbond
    
    cond = per_diff > tol 
    while cond
	bondVmax = 1.025 * bondVmax
        bondpr = get_cvm_bond_price(svm, svm.m, svm.pm.sigmah;
                                    Vt=bondVmax, vb=vb,
                                    mu_b=mu_b, c=c, p=p)
	per_diff = (rfbond - bondpr) / rfbond
	cond = per_diff > tol
    end
    
    if print
	println(string("Bond Percentage Difference: ", per_diff))
    end
    
    return bondVmax
end


# * Debt Pricing Functions ###################################
# include("debtpr_funs.jl") 
function get_pv_rfdebt(svm; mu_b::Float64=NaN,
                       m::Float64=NaN,
                       c::Float64=NaN,
                       p::Float64=NaN)

    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)
    
    return rf_debt(mu_b, m, c, p, 
                   svm.pm.r, svm.pm.xi, svm.pm.kappa)
end


function get_cvm_debt_price(svm, vb::Float64,
                            sigma::Float64; 
                            Vt::Float64=NaN,
                            mu_b::Float64=NaN,
                            m::Float64=NaN,
                            c::Float64=NaN,
                            p::Float64=NaN,
                            N1::Int64=100, N2::Int64=10^4)

    # ####################################
    # ######## Extract Parameters ########
    # ####################################
    if isnan(Vt)
        Vt=svm.pm.V0
    end

    mu_b, m, c, p = get_k_struct(svm; mu_b=mu_b, m=m, c=c, p=p)

    println(string(" mu_b = ", mu_b))
    println(string(" m = ", m))
    println(string(" c = ", c))
    println(string(" p = ", p))
    # ####################################
    
    # vb = get_cvm_vb(svm, sigma; c=c, p=p)
    v = log(Vt/float(vb))
	
    # Create time-to-maturity grids
    _, ttm_grid1 = grid_creator(0.01, m, N1)
    dt2, ttm_grid2 = grid_creator(0.0, m, N2)


    # Compute Bond Prices
    bond_vec = @spawn [get_cvm_bond_price(svm, ttm, sigma;
                                          Vt=Vt, vb=vb,
                                          mu_b=mu_b, m=m,
                                          c=c, p=p)
                       for ttm in ttm_grid1]

    # Interpolate	
    bond_interp_sitp = Dierckx.Spline1D(ttm_grid1, fetch(bond_vec); k=3, bc="extrapolate")

    # Refine Bond Price Grid
    bond_vec2 = @spawn [bond_interp_sitp(ttm) for ttm in ttm_grid2]
    
    # Integrate
    return mu_b * sum(fetch(bond_vec2)) * dt2
end


function get_svm_debt_price(svm, vbl::Float64; 
                            Vt::Float64=NaN,
                            mu_b::Float64=NaN,
                            c::Float64=NaN,
                            p::Float64=NaN,
                            ttmN0::Int64=10^2,
                            ttmN::Int64=10^4)

    # ####################################
    # ######## Extract Parameters ########
    # ####################################
    if isnan(Vt)
        Vt=svm.pm.V0
    end

    mu_b, _, c, p = get_k_struct(svm; mu_b=mu_b, m=NaN, c=c, p=p)
    # ####################################
    
    # Create  time-to-maturity grid
    _, ttm_grid = grid_creator(minimum(svm.bs.ttmgrid), 
                               maximum((svm.bs.ttmgrid)), ttmN0)

    # Get Bond Prices at Different Maturities
    bondpr_vec = fetch(@spawn [get_svm_bond_price(svm, vbl, ttm;
                                                  Vt=Vt, mu_b=mu_b, c=c, p=p)
                               for ttm in ttm_grid])
    bondpr_sitp = Dierckx.Spline1D(ttm_grid, bondpr_vec; k=3, bc="extrapolate")
    
    # Refine time-to-maturity grid
    dt2, ttm_grid2 = grid_creator(0.0, svm.m, ttmN)

    # Bond Prices
    bondpr_vec_ref = fetch(@spawn [bondpr_sitp(ttm) for ttm in ttm_grid2])

    return mu_b * sum(bondpr_vec_ref) * dt2
end

# * Bond Yields 
function get_svm_bond_price(svm, vbl::Float64,
                            ttm::Float64;
                            Vt::Float64=NaN,
                            mu_b::Float64=NaN,
                            c::Float64=NaN,
                            p::Float64=NaN,
                            ftype::String="bf")
    # ####################################
    # ######## Extract Parameters ########
    # ####################################
    if isnan(Vt)
        Vt = svm.pm.V0
    end

    mu_b, _, c, p = get_k_struct(svm; mu_b=mu_b, m=NaN, c=c, p=p)
    # ####################################
    
    # Maximum vt
    vtmax = svm.bi.vtmax 

    vt, vbh, vbhl = compute_jbpr_inputs(svm, Vt, vbl, mu_b, c, p) 
    # ####################################
    if vt <= 0
        return on_default_payoff(vt, vbl, ttm,
                                 mu_b, svm.m, c, p,
                                 svm.pm.r, svm.pm.xi,
                                 svm.pm.kappa, svm.pm.alpha)
    elseif vt > minimum(svm.bs.vtgrid)
        rfbond = rfbond_price(ttm, c, p,
                              svm.pm.r, svm.pm.xi, svm.pm.kappa)

        if vt > vtmax
               return rfbond 
        end

        # Maturity or Default prior to Volatility Shock:
        cf0 = no_vol_shock_cf_pv(vt, vbl, ttm,
                                 mu_b, svm.m,
                                 c, p,
                                 svm.pm.sigmal,
                                 svm.pm.r, svm.pm.gross_delta,
                                 svm.pm.xi, svm.pm.kappa,
                                 svm.pm.alpha, svm.pm.lambda)
       
        f11, f12, f13, f21, f22 = get_interp_values(svm, vt, ttm, vbhl; ftype=ftype)
        
        # Volatility Shock Prior to Maturity:
        cf1 = c/get_rdisc(svm) * f11 +
              (p - c/get_rdisc(svm)) * f12 +
              (svm.pm.alpha * vbh/float(mu_b * svm.m) - c/get_rdisc(svm)) * f13
       
        cf2 = c/get_rdisc(svm) * f21 + (p - c/get_rdisc(svm)) * f22 
       
        return min(cf0 + cf1 + cf2, rfbond)
    else
        # In this case, firm is so close to default that I just
        # return the CVM bond price for sigma = sigmal:
        return get_cvm_bond_price(svm, ttm, svm.pm.sigmal;
                            Vt=Vt, vb=vbl, mu_b=mu_b, c=c, p=p)
    end
end


# * BOND YIELDS
function get_bond_yield(svm;
                        KS::Dict{Symbol, Float64}=Dict{Symbol, Float64}(),
                        min_yield=1e-3,
                        max_yield=5.,
                        N=10^5,
                        ftype::String="bf")
    
    # Capital Structure Variables
    ks_vec = svm.model == "cvm" ? [:mu_b, :m, :c, :p] : [:mu_b, :m, :c, :p, :vbl]
    
    # Check for missing entries and NaN values in the KS dictionary:
    ks_vec_entries = [x for x in ks_vec if (x in keys(KS))]
    missing = (ks_vec_entries == ks_vec) ? any([isnan(KS[x]) for x in ks_vec_entries]) : true
    
    if missing
        println("Using Model Object Optimal Capital Structure paramters")
        
        if any([getfield(svm.optKS, x)==nothing 
                        for x in ks_vec])
            println("Missing optimal capital structure parameters. Exiting...")
            return
        end
        
        for x in ks_vec
            KS[x] = getfield(svm.optKS, x)
        end
    end
    
    # Yield - Bond Price Function
    ybpr(y) = (KS[:c]/y) * (1 - exp(-y * KS[:m])) + KS[:p] * exp(-y * KS[:m])
    
    # Vector of Yield Candidates
    yield_vec = range(min_yield, stop=max_yield, length=N)
    
    # Vector of Bond Prices
    ybpr_vec = [ybpr(y) for y in yield_vec]
    
    
    # Newly Issued Bond Price
    nibpr = NaN
    if svm.model == "cvm"
        nibpr = get_cvm_bond_price(svm, KS[:m], svm.pm.sigmal; 
                                   mu_b=KS[:mu_b], m=KS[:m],
                                   c=KS[:c], p=KS[:p])
    elseif svm.model == "svm"
        nibpr = get_svm_bond_price(svm, KS[:vbl], KS[:m]; 
                                   mu_b=KS[:mu_b], c=KS[:c], 
                                   p=KS[:p],ftype=ftype)
    end
    
    # Bond Yield
    return yield_vec[argmin(abs.(ybpr_vec .- nibpr))]
end


function get_bond_spread(svm;
                        KS::Dict{Symbol, Float64}=Dict{Symbol, Float64}(),
                        min_yield=1e-3,
                        max_yield=5.,
                        N=10^5,
                        ftype::String="bf")

  byield = get_bond_yield(svm; KS=KS, min_yield=min_yield,
                          max_yield=max_yield, N=N, ftype=ftype)

  # Spreads in Basis Points
  return (byield - svm.pm.r) * 10^4
end

# * END MODULE
end
