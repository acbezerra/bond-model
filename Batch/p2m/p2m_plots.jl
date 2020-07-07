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
plt_fi = false
plt_mp = false
plt_pool = false
plt_sep = false
plt_ep = true
plt_dual = true

# Modularity {{{1
main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")

# * Load Modules and Prepare Results {{{1
if run_mod
  p2m_init_path = string(main_path, "/Julia/Batch/p2m/") 
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


# Full Information Payoffs {{{1
if plt_fi
    xvar = :q
    yvar = :mu_s
    ft = :r
    fidf[!, :q] .= fidf[:, :rq]
    # mu_s_vals = unique(pooldf[:, :mu_s])
    zvars = [:debt, :eq, :fv, :mbr]

    for z in zvars
        eq_type = :fi
        zvar = Symbol(ft, :_fi_, z)
        file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z;
                                                    ft=Symbol(ft, :t))

        fi_fig = p2m.plot_2pm_iso_contour_curves(ext_fid[:X], 
                                                 ext_fid[:Y],
                                                 ext_fid[zvar][:zvals]';
                                                 iso_cmap=p2m.iso_cmaps[eq_type],
                                                 fig_title=plot_title,
                                                 file_path_name=file_path_name)
    end
end

# Misrepresentation Payoffs  {{{1
if plt_mp
    xvar = :q
    yvar = :mu_s
    ft = :r
    eq_type = :mp
    for z in [:fv, :mbr]
        zvar = Symbol(ft, :_, eq_type, :_, z)
        file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z;
                                                    ft=Symbol(ft, :t))
        p2m.plot_2pm_iso_contour_curves(ext_fid[:X], ext_fid[:Y], ext_fid[zvar][:zvals]';
                                        iso_cmap=p2m.iso_cmaps[eq_type],
                                        fig_title=plot_title,
                                        file_path_name=file_path_name)
    end

    for z in [:fv, :mbr]
        for zdf in [:diff, :perc_diff]
            zvar = Symbol(ft, :_mp_fi_, z, :_, zdf)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z; zdf=zdf,
                                                        ft=Symbol(ft, :t))
            p2m.plot_2pm_iso_contour_curves(ext_fid[:X], ext_fid[:Y], ext_fid[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)

        end
    end
end


# Pooling Payoffs  {{{1
if plt_pool
    xvar = :q
    yvar = :mu_s
    eq_type = :pool
    for ft in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
            zvar = Symbol(ft, :_, eq_type, :_, z)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z;
                                                    ft=Symbol(ft, :t))
            p2m.plot_2pm_iso_contour_curves(ext_poold[:X], ext_poold[:Y], ext_poold[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end

    # Differentials
    for ft in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
           for zdf in [:diff, :perc_diff] 
                zvar = Symbol(ft, :_pool_fi_, z, :_, zdf)
                file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
                plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z; zdf=zdf,
                                                            ft=Symbol(ft, :t))
                p2m.plot_2pm_iso_contour_curves(ext_poold[:X], ext_poold[:Y], ext_poold[zvar][:zvals]';
                                                iso_cmap=p2m.iso_cmaps[eq_type],
                                                fig_title=plot_title,
                                                file_path_name=file_path_name)
            end
        end
    end
end

if plt_sep
    # Separating Payoffs {{{1
    mu_s_vals = unique(pooldf[:, :mu_s])
    xvar = :q
    yvar = :mu_s
    eq_type = :sep

    for ft in [:s, :r]
        for z in [x for x in p2m.zfuns if x != :bpr]
            zvar = Symbol(ft, :_, eq_type, :_, z)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z;
                                                    ft=Symbol(ft, :t))
            p2m.plot_2pm_iso_contour_curves(sepd[:X], sepd[:Y], sepd[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end


    # Differentials
    for z in [x for x in p2m.zfuns if x != :bpr]
       for zdf in [:diff, :perc_diff] 
            zvar = Symbol(:s_sep_fi_, z, :_, zdf)
            file_path_name = p2m.get_2pm_contour_plot_path_name(zvar)
            plot_title = p2m.get_2pm_contour_plot_title(sf, eq_type, z; zdf=zdf,
                                                        ft=:st)
            p2m.plot_2pm_iso_contour_curves(sepd[:X], sepd[:Y], sepd[zvar][:zvals]';
                                            iso_cmap=p2m.iso_cmaps[eq_type],
                                            fig_title=plot_title,
                                            file_path_name=file_path_name)
        end
    end
end

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

# JEQ Plots {{{1
cat_cmaps = ["gnuplot", "GnBu"]
zvars = [:s_fv, :r_mbr]

iso_cmap = "cool"
cat_cmap = "bone_r"

# EP Equilibria {{{2
if plt_ep
    for i in 1:2
        zz = zvars[i]
        
        ft = Symbol(Symbol(split(String(zz), "_")[1]), :t)
        zvar = Symbol(split(String(zz), "_")[2])
        plot_zvar = Symbol(Symbol(split(String(zz), "_")[1]), :_ep_, zvar)
        file_path_name = p2m.get_2pm_contour_plot_path_name(plot_zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, :ep, zvar;
                                                    ft=ft)
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

# EP v.s. OTC Equilibria {{{2
if plt_dual
    zvars = [:s_fv, :r_mbr]
    iso_cmap = "cool"
    cat_cmap = "bone_r"

    rbdisc_otc = sf.rt.rf + ilc_otc
    for i in 1:2
        zz = zvars[i]
        
        ft = Symbol(Symbol(split(String(zz), "_")[1]), :t)
        zvar = Symbol(split(String(zz), "_")[2])
        plot_zvar = Symbol(Symbol(split(String(zz), "_")[1]), :_dual_, zvar)
        file_path_name = p2m.get_2pm_contour_plot_path_name(plot_zvar)
        plot_title = p2m.get_2pm_contour_plot_title(sf, :dual, zvar;
                                                    ft=ft, rbdisc_otc=rbdisc_otc,
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
