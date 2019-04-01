

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
