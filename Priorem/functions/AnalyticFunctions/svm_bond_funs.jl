# # SVM Bond and Debt Pricing Functions - Python
# 
#     function approx_bond_int(vt, v_upper_thresh, ttm, c, p, r, xi, k, _lambda, sigmal, gross_delta)
#         return _lambda * exp(-rdisc_pvs(r, xi, k, _lambda) * ttm) *
#                rfbond_price(ttm, c, p, r, xi, k) *
#                psi_v_td(vt, v_upper_thresh, ttm, sigmal, r, gross_delta)
#     end
# 
#     function bond_pr_vol_cf(vt, ttm, vmax, vb, vbh,
#                             c, p, r, gross_delta,
#                             xi, k, alpha, pi,
#                             _lambda, sigmal, sigmah)
#         tauN = 10^3
#         vN = 10^3
#         dt1 = (ttm - 1e-4)/tauN
#         dv = vmax/vN
#         cf1_mat = @spawn [vol_shock_cf_integrand(vt, v, v - log(vbh / vb),
#                                                  vb, vbh, u, ttm,
#                                                  c, p, sigmal, sigmah, r,
#                                                  gross_delta, xi, k,
#                                                  alpha, pi, _lambda) *
#                           vb * exp(v) * (v > 0) for u=linspace(1e-4, ttm, tauN),
#                                                     v=linspace(0.0, vmax, vN)]
# 
#         N = 10^4
#         dt2 = ttm/N
#         cf2_vec = @spawn [approx_bond_int(vt, vmax, tau,  c, p, r,
#                                            xi, k, _lambda, sigmal, gross_delta)
#                           for tau=linspace(dt2, ttm, N)]
#         return sum(fetch(cf1_mat)) * dt1 * dv + sum(fetch(cf2_vec)) * dt2
#     end
# 
#     # Objective is to find V such that the value of the
#     # newly-issued (tau = m) risky bond price when sigma = sigmah
#     # is sufficiently close to the credit-risk-free bond price.
#     function get_bond_vmax(V0, m, c, p, sigmah, r,
#                             gross_delta, iota, xi, k,
#                             alpha, pi, print=false)
#         bondVmax = 1.25*V0
#         vb = zhi_vb(m, c, p, sigmah, r, gross_delta, iota, xi, k, alpha, pi)
# 
#         vmax = log(bondVmax/vb)
#         rfbond = rfbond_price(m, c, p, r, xi, k)
#         bondpr = zhi_bond_price(vb, vmax, m, m, c, p, sigmah, r,
#                                     gross_delta, xi, k, alpha, pi)
# 
#         per_diff = (rfbond - bondpr) / rfbond
# 
#         cond = per_diff > 1e-3
#         while cond
#             bondVmax = 1.025 * bondVmax
#             vmax = log(bondVmax / vb)
#             bondpr = zhi_bond_price(vb, vmax, m, m, c, p, sigmah, r,
#                                 gross_delta, xi, k, alpha, pi)
#             per_diff = (rfbond - bondpr) / rfbond
#             cond = per_diff > 1e-3
#         end
# 
#         if print
#             println(string("Bond Percentage Difference: ", per_diff))
#         end
# 
#         return bondVmax
#     end
# 
# 
#     function bond_pr_main(vt, ttm, vb, vbh,
#                           m, c, p,
#                           r, gross_delta, iota,
#                           xi, k, alpha, pi,
#                           _lambda, sigmal, sigmah, vmax=NaN)
# 
#         if vt <= 0
#             return on_default_payoff(vt, vb, ttm,
#                                      m, c, p, r, xi,
#                                      k, alpha)
#         else
#             Vt = vb*exp(vt)
#             # Compute vmax
#             if isnan(vmax)
#                 bondVmax = get_bond_vmax(Vt, m, c, p, sigmah, r,
#                                          gross_delta, iota, xi,
#                                          k, alpha, pi)
#                 vmax = log(bondVmax/vb)
#             end
# 
#             if vt > vmax
#                 return rfbond_price(ttm, c, p, r, xi, k)
#             else
#                 cf0 = no_vol_shock_cf_pv(vt, vb, ttm,
#                                          m, c, p, sigmal,
#                                          r, gross_delta,
#                                          xi, k, alpha, _lambda)
#                 # println(string("cf0: ", cf0))
# 
#                 cf1 = bond_pr_vol_cf(vt, ttm, vmax, vb, vbh,
#                                      c, p, r, gross_delta,
#                                      xi, k, alpha, pi,
#                                      _lambda, sigmal, sigmah)
#                 # println(string("cf1: ", cf1))
# 
#                 return min(cf0 + cf1, rfbond_price(ttm, c, p, r, xi, k))
#             end
#         end
#     end
# 
#     function get_debt_price(vt, vb, vbh, ttm_grid_size,
#                             m, c, p,
#                             r, gross_delta, iota,
#                             xi, k, alpha, pi,
#                             _lambda, sigmal, sigmah)
#         t0 = .0001
# 
#         # Grid for Interpolation
#         ttm_grid = linspace(t0, m, ttm_grid_size)
# 
#         # Since vmax does not change with maturity, compute
#         Vt = vb*exp(vt)
#         # it only once:
#         bondVmax = get_bond_vmax(Vt, m, c, p, sigmah, r,
#                                  gross_delta, iota, xi,
#                                  k, alpha, pi, true)
#         vmax = log(bondVmax/vb)
# 
#         # Compute Bond Prices for Different Maturities
#         tmp = @spawn [bond_pr_main(vt, ttm, vb, vbh,
#                                     m, c, p,
#                                     r, gross_delta, iota,
#                                     xi, k, alpha, pi,
#                                     _lambda, sigmal, sigmah, vmax) for ttm=ttm_grid]
#         ttmvec = fetch(tmp)
# 
#         # Interpolate
#         itp = interpolate(ttmvec, BSpline(Cubic(Line())), OnGrid())
#         sitp = Interpolations.scale(itp, ttm_grid) # Scale
# 
#         # Refined Grid
#         ttm_grid_refined_size = 10^4
#         ttm_grid_refined = linspace(t0, m, ttm_grid_refined_size)
#         dt = (m - t0)/ttm_grid_refined_size
# 
#         # Integrate
#         return (1/m) * sum([sitp[x] for x=ttm_grid_refined]) * dt
#     end
# 
#     function pvvb_bond_surf(V0, pmin, pmax, pN, vN,
#                             vbmin, vbmax, vbN,
#                             m, c, p,
#                             r, gross_delta, iota,
#                             xi, k,
#                             alpha, pi,
#                             _lambda, sigmal, sigmah)
# 
#         bondVmax = get_bond_vmax(V0, m, c, p, sigmah, r,
#                                  gross_delta, iota, xi, k,
#                                  alpha, pi)
# 
#         vbh = zhi_vb(m, c, p, sigmah, r, gross_delta, iota,
#                      xi, k, alpha, pi)
# 
# 
#         pgrid = linspace(pmin, pmax, pN)
#         vgrid = linspace(0, log(bondVmax/vbmin), vN)
#         vbgrid = linspace(vbmin, vbmax, vbN)
# 
#         cube_future = @spawn [bond_pr_main(v, m, vb, vbh, m, c, p,
#                                  r, gross_delta, iota,xi, k, alpha, pi, _lambda,
#                                  sigmal, sigmah, log(bondVmax/vb)) for p=pgrid, v=vgrid, vb=vbgrid]
# 
#        return fetch(cube_future)
# 
#     end
# 
