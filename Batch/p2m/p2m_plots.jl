# vim: set fdm=marker :

# * SYS & Other Args
# println(string("ARGUMENTS: ", ARGS))
# # ################ SYS ARGUMENTS ################
# # Position in the Safe Firm Measure array:
# run_mod = parse(Bool, ARGS[1])
# plt_fi = parse(Bool, ARGS[2])
# plt_mp = parse(Bool, ARGS[3])
# plt_pool = parse(Bool, ARGS[4])
# plt_sep = parse(Bool, ARGS[5])
# plt_ep = parse(Bool, ARGS[6])
# plt_dual = parse(Bool, ARGS[7])

run_mod = false
plt_fi = true
plt_mp = true
plt_pool = true
plt_sep = true
plt_ep = true
plt_dual = true

# Modularity {{{1
main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "bond-model/modules/")
# }}}1
# * Load Modules and Prepare Results {{{1
if run_mod
  p2m_init_path = string(main_path, "/bond-model/Batch/p2m/")
  include(string(joinpath(p2m_init_path, "p2m_initializer"), ".jl"))
else
    modls = ["2period"] # "ModelPlots",
    for modl in modls
      include(string(joinpath(module_path, modl), ".jl"))
    end
end


sf = resd[:sf]
rf = resd[:rf]
s_otc_fv = resd[:s_otc_fv]
ptr = resd[:ptr]
fid = resd[:fid]
poold = resd[:poold]
sepd = resd[:sepd]
pltd = resd[:pltd]
ext_fid = resd[:ext_fid]
ext_poold = resd[:ext_poold]
ext_pltd = resd[:ext_pltd]
ilc_otc = resd[:ilc_otc]
theta= 1.
otcdf = res[:otcdf]
fidf = res[:ext_fid][:df]
# }}}1
# Initial Comparison: FV v.s. MBR {{{1
qr = .5
p2m.plot_fi_fv_vs_mbr(df, qr; )
# }}}1
# Full Information Payoffs {{{1
if plt_fi
    xvar = :q
    yvar = :mu_s
    fty = :r
    fidf[!, :q] .= fidf[:, :rq]
    # mu_s_vals = unique(pooldf[:, :mu_s])
    zvars = [:debt, :eq, :fv, :mbr]
    plt_eq_type = :fi

    for z in zvars
        zvar = Symbol(fty, :_fi_, z)
        file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z;
                                                    ft=Symbol(fty, :t))

        fi_fig = p2m.plot_2pm_iso_contour_curves(ext_fid[:X],
                                                 ext_fid[:Y],
                                                 ext_fid[zvar][:zvals]';
                                                 iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                                 fig_title=plot_title,
                                                 file_path_name=file_path_name)
    end
end
# }}}1
# Misrepresentation Payoffs  {{{1
if plt_mp
    xvar = :q
    yvar = :mu_s
    fty = :r
    plt_eq_type = :mp
    for z in [:fv, :mbr]
        zvar = Symbol(fty, :_, plt_eq_type, :_, z)
        file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z;
                                                    ft=Symbol(fty, :t))
        p2m.plot_2pm_iso_contour_curves(ext_fid[:X], ext_fid[:Y], ext_fid[zvar][:zvals]';
                                        iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                        fig_title=plot_title,
                                        file_path_name=file_path_name)
    end

    for z in [:fv, :mbr]
        for zdf in [:diff, :perc_diff]
            zvar = Symbol(fty, :_mp_fi_, z, :_, zdf)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z; zdf=zdf,
                                                        ft=Symbol(fty, :t))
            p2m.plot_2pm_iso_contour_curves(ext_fid[:X], ext_fid[:Y], ext_fid[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)

        end
    end
end
# }}}1
# OTC Covenant Misrepresentation {{{1
tmp = p2m.otc_misrep_mbr(sf, rf, fidf, otcdf, ilc_otc; theta = 1.)
p2m.plot_r_otc_misrep_mbr(rf, tmp)
# }}}1 
# Pooling Prelim Analysis {{{1
qgrid = [.5]
rf = p2m.setq(rf, q=qgrid[1])
mu_s_grid = .015:.015:.975
fidf = p2m.get_opt_lev(sf, rf; )
pdff = vcat([p2m.get_pool_res(sf, rf, qgrid, mu_s; fidf=fidf) for mu_s in mu_s_grid]...)
fg = p2m.plot_pool_mub_mbr(pdff)
# }}}1
# Pooling Payoffs  {{{1
if plt_pool
    xvar = :q
    yvar = :mu_s
    plt_eq_type = :pool
    for fty in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
            zvar = Symbol(fty, :_, plt_eq_type, :_, z)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z;
                                                    ft=Symbol(fty, :t))
            p2m.plot_2pm_iso_contour_curves(ext_poold[:X], ext_poold[:Y], ext_poold[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end

    # Differentials
    for fty in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
           for zdf in [:diff, :perc_diff]
                zvar = Symbol(fty, :_pool_fi_, z, :_, zdf)
                file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
                plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z; zdf=zdf,
                                                            ft=Symbol(fty, :t))
                p2m.plot_2pm_iso_contour_curves(ext_poold[:X], ext_poold[:Y], ext_poold[zvar][:zvals]';
                                                iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                                fig_title=plot_title,
                                                file_path_name=file_path_name)
            end
        end
    end
end
# }}}1
# Separating Payoffs {{{1
if plt_sep
    mu_s_vals = unique(pooldf[:, :mu_s])
    xvar = :q
    yvar = :mu_s
    plt_eq_type = :sep

    for fty in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
            zvar = Symbol(fty, :_, plt_eq_type, :_, z)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z;
                                                    ft=Symbol(fty, :t))
            p2m.plot_2pm_iso_contour_curves(sepd[:X], sepd[:Y], sepd[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end

    # Differentials
    for z in [x for x in p2m.zfuns if x != :bpr]
       for zdf in [:diff, :perc_diff]
            zvar = Symbol(:s_sep_fi_, z, :_, zdf)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, plt_eq_type, z; zdf=zdf,
                                                        ft=:st)
            p2m.plot_2pm_iso_contour_curves(sepd[:X], sepd[:Y], sepd[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[plt_eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end
end
# }}}1
# Data Points {{{1
# Dictionary of Data Points
ptd = Dict{Symbol, Dict{Symbol, Any}}(
        [pt => Dict{Symbol, Any}(:q => ptr[pt][:q],
                                 :mu_s => ptr[pt][:mu_s],
                                 :marker => "o",
                                 :color => "mediumspringgreen")
        for pt in keys(ptr)]
)
ptd_ep = Dict([pt => ptd[pt] for pt in [:p1, :p2]])

# JEQ Plots
cat_cmaps = ["gnuplot", "GnBu"]
zvars = [:s_fv, :r_mbr]

iso_cmap = "cool"
cat_cmap = "bone_r"

# EP Equilibria {{{2
if plt_ep
    for i in 1:2
        zz = zvars[i]

        z_fty = Symbol(Symbol(split(String(zz), "_")[1]), :t)
        zvar = Symbol(split(String(zz), "_")[2])
        plot_zvar = Symbol(Symbol(split(String(zz), "_")[1]), :_ep_, zvar)
        file_path_name = p2m.get_2pm_contour_plot_path_name(plot_zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, :ep, zvar;
                                                    ft=z_fty)
        ep_fig = p2m.plot_2pm_equilibria_iso_contour_curves(ext_pltd[:X], ext_pltd[:Y],
                                                            ext_pltd[zz],
                                                            ext_pltd[:bool_Z];
                                                            fig_title=plot_title,
                                                            file_path_name=file_path_name,
                                                            iso_cmap=iso_cmap,
    #                                                       iso_levels=15,
                                                            cat_cmap=cat_cmap)

        file_path_name_pt = string(split(file_path_name, ".")[1], "_pt.",
                                   split(file_path_name, ".")[2])
        ep_fig = p2m.plot_2pm_equilibria_iso_contour_curves(ext_pltd[:X], ext_pltd[:Y],
                                                            ext_pltd[zz],
                                                            ext_pltd[:bool_Z];
                                                            ptd=ptd_ep,
                                                            fig_title=plot_title,
                                                            file_path_name=file_path_name_pt,
                                                            iso_cmap=iso_cmap,
    #                                                       iso_levels=15,
                                                            cat_cmap=cat_cmap)
    end
end
# }}}2
# EP v.s. OTC Equilibria {{{2
if plt_dual
    zvars = [:s_fv, :r_mbr]
    iso_cmap = "cool"
    cat_cmap = "bone_r"

    rbdisc_otc = sf.rt.rf + ilc_otc
    for i in 1:2
        zz = zvars[i]

        z_fty = Symbol(Symbol(split(String(zz), "_")[1]), :t)
        zvar = Symbol(split(String(zz), "_")[2])
        plot_zvar = Symbol(Symbol(split(String(zz), "_")[1]), :_dual_, zvar)
        file_path_name = p2m.get_2pm_contour_plot_path_name(plot_zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, :dual, zvar;
                                                    ft=z_fty, rbdisc_otc=rbdisc_otc,
                                                    s_otc_fv=s_otc_fv)


        p2m.plot_2pm_equilibria_iso_contour_curves(ext_pltd[:X], ext_pltd[:Y],
                                                   ext_pltd[zz],
                                                   ext_pltd[:bool_OTC_EP];
                                                   fig_title=plot_title,
                                                   file_path_name=file_path_name,
                                                   iso_cmap=iso_cmap,
                                                   iso_levels=15,
                                                   cat_cmap=cat_cmap)

        file_path_name_pt = string(split(file_path_name, ".")[1], "_pt.",
                                   split(file_path_name, ".")[2])
        p2m.plot_2pm_equilibria_iso_contour_curves(ext_pltd[:X], ext_pltd[:Y],
                                                   ext_pltd[zz],
                                                   ext_pltd[:bool_OTC_EP];
                                                   ptd=ptd,
                                                   fig_title=plot_title,
                                                   file_path_name=file_path_name_pt,
                                                   iso_cmap=iso_cmap,
                                                   iso_levels=15,
                                                   cat_cmap=cat_cmap)
    end
end
# }}}2
# }}}1
