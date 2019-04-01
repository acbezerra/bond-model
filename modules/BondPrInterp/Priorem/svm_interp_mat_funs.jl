
# ##################################################################
# ######################### J2: f21, f22 ###########################
# ################################################################## 
function f21_mat_inputs(svm, ttmgrid, vbhlgrid; vmax=1.2, uN=10^3)

    f21_int_vec = @spawn [f21_int(log(svm.pm.V0/float(svm.vbh)) + log(vbhl),
                                  vmax, svm.pm.lambda, svm.pm.sigmal, 
                                  svm.pm.r, svm.pm.gross_delta,
                                  svm.pm.xi, svm.pm.kappa; 
				  ttm=tau, N=uN) for tau in ttmgrid,
                                                    vbhl in vbhlgrid]

    return fetch(f21_int_vec)
end


function f22_mat_inputs(svm, ttmgrid, vbhlgrid; vmax=1.2, uN=10^3)

    f22_int_vec = @spawn [f22_int(log(svm.pm.V0/float(svm.vbh)) + log(vbhl),
                                  vmax, svm.pm.lambda, svm.pm.sigmal, 
                                  svm.pm.r, svm.pm.gross_delta,
                                  svm.pm.xi, svm.pm.kappa; 
				  ttm=tau, N=uN) for tau in ttmgrid,
                                                    vbhl in vbhlgrid]
    return fetch(f22_int_vec)
end


# ##################################################################
# ########################### f1_int ###############################
# ################################################################## 
function f11_mat_inputs(svm, ttmgrid, vbhlgrid;
                        vmax=1.2, vN=10^3, uN=10^3)		

    dv, vgrid = grid_creator(0.0, vmax, vN)

    f11v_vec = @spawn [f11v_int(log(svm.pm.V0/float(svm.vbh)) + log(vbhl),
                                grid_creator(0.0, u, uN)[2],
                                vgrid, svm.pm.lambda, svm.pm.sigmal,
                                svm.pm.r, svm.pm.gross_delta,
                                svm.pm.xi, svm.pm.kappa) for u in ttmgrid,
                                                          vbhl in vbhlgrid]

    return fetch(f11v_vec)
end


# ##################################################################
# ########################### f2_int ###############################
# ################################################################## 
function f12_mat_inputs(svm, ttmgrid, vbhlgrid;
		        vmax=1.2, vN=10^3, uN=10^3)
    
    dv, vgrid = grid_creator(0.0, vmax, vN)

    f12_surf = @spawn [f12_int(log(svm.pm.V0/float(svm.vbh)) + log(vbhl),
                               vbhl, grid_creator(0.0, u, uN)[2],
		               vgrid, svm.pm.lambda,
                               svm.pm.sigmal, svm.pm.sigmah,
                               svm.pm.r, svm.pm.gross_delta,
                               svm.pm.xi, svm.pm.kappa) for u in ttmgrid,
                                                         vbhl in vbhlgrid]

    return fetch(f12_surf)
end


# ##################################################################
# ########################### f3_int ###############################
# ################################################################## 
function f13_mat_inputs(svm, ttmgrid, vbhlgrid;
		        vmax=1.2, vN=10^3, uN=10^3)

    dv, vgrid = grid_creator(0.0, vmax, vN)

    f13_surf = @spawn [f13_int(log(svm.pm.V0/float(svm.vbh)) + log(vbhl),
                               vbhl, grid_creator(0.0, u, uN)[2],
                               vgrid, svm.pm.lambda,
                               svm.pm.sigmal, svm.pm.sigmah,
                               svm.pm.r, svm.pm.gross_delta,
                               svm.pm.xi, svm.pm.kappa) for u in ttmgrid,
                                                         vbhl in vbhlgrid]
     

    return fetch(f13_surf)
end


# ##################################################################
# ########################### Interp ###############################
# ################################################################## 
function mat_interp_inputs(svm;
                          ttmgrid=nothing, ttm_max=svm.m, ttmN=20,
                          vbhlgrid=nothing, vbhlmin=.75, vbhlmax=1.25, vbhlN=10,
                          vmax=1.2, vN=10^3, uN=10^3)
   
    # ####################################
    # ##### Form Interpolation Grids #####
    # ####################################
    if ttmgrid==nothing
        _, ttmgrid = grid_creator(0.0, ttm_max, ttmN)
    end

    if vbhlgrid==nothing
        _, vbhlgrid=grid_creator(vbhlmin, vbhlmax, vbhlN)
    end
    # ####################################

    # Set Post-Volatility Shock Default Barrier
    svm.vbh = get_cvm_vb(svm, svm.pm.sigmah)

    f11 = @spawn f11_mat_inputs(svm, ttmgrid, vbhlgrid;
                                vmax=vmax, vN=vN, uN=uN)

    f12 = @spawn f12_mat_inputs(svm, ttmgrid, vbhlgrid;
                                vmax=vmax, vN=vN, uN=uN)
    
    f13 = @spawn f13_mat_inputs(svm, ttmgrid, vbhlgrid;
                                vmax=vmax, vN=vN, uN=uN)

    f21 = @spawn f21_mat_inputs(svm, ttmgrid, vbhlgrid;
                                vmax=vmax, uN=uN)

    f22 = @spawn f22_mat_inputs(svm, ttmgrid, vbhlgrid;
                                vmax=vmax, uN=uN)

    f11mat = fetch(f11)
    f12mat = fetch(f12)
    f13mat = fetch(f13)
    f21mat = fetch(f21)
    f22mat = fetch(f22)
    
    mat_input = Dict("vmax" => vmax,
                     "ttm_max" => ttm_max,
                     "ttmgrid"=> ttmgrid,
                     "vbhlgrid"=> vbhlgrid,
                     "vbhlmin" => vbhlmin,
                     "vbhlmax" => vbhlmax, 
                     "f21_input"=> f21mat,
                     "f22_input"=> f22mat,
		     "f11_input"=> f11mat,
		     "f12_input"=> f12mat,
		     "f13_input"=> f13mat)

    return mat_input
end


function mat_input_interp(bprd)

    f11_itp = interpolate(bprd["f11_input"], BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f11_sitp = Interpolations.scale(f11_itp, bprd["ttmgrid"], bprd["vbhlgrid"])

    f12_itp = interpolate(bprd["f12_input"], BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f12_sitp = Interpolations.scale(f12_itp, bprd["ttmgrid"], bprd["vbhlgrid"])

    f13_itp = interpolate(bprd["f13_input"], BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f13_sitp = Interpolations.scale(f13_itp, bprd["ttmgrid"], bprd["vbhlgrid"])

    # if length(bprd["vbhlgrid"]) > 1
    #     f12_sitp = Interpolations.scale(f12_itp, bprd["ttmgrid"], bprd["vbhlgrid"])
    #     f13_sitp = Interpolations.scale(f13_itp, bprd["ttmgrid"], bprd["vbhlgrid"])
    # else
    #     f12_sitp = Interpolations.scale(f12_itp, bprd["ttmgrid"])
    #     f13_sitp = Interpolations.scale(f13_itp, bprd["ttmgrid"])
    # end
    
    f21_itp = interpolate(bprd["f21_input"], BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f21_sitp = Interpolations.scale(f21_itp, bprd["ttmgrid"], bprd["vbhlgrid"])
    
    f22_itp = interpolate(bprd["f22_input"], BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f22_sitp = Interpolations.scale(f22_itp, bprd["ttmgrid"], bprd["vbhlgrid"])

    return Dict("f11"=> f11_sitp,
                "f12"=> f12_sitp,
                "f13"=> f13_sitp,
                "f21"=> f21_sitp,
                "f22"=> f22_sitp)
end
