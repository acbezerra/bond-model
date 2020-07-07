main_path = "/home/artur/BondPricing"
run_modularity = true
if run_modularity
    plot_script_path = string(main_path, "/Julia/Batch/plot_scripts/_contour")
    plot_script_name = "contour_plots_initializer.jl"
    include(string(plot_script_path, "/", plot_script_name))
else
    main_path = "/home/artur/BondPricing"
    module_path = string(main_path, "/", "Julia/modules/")
    modls = ["Batch", "ModelObj", 
             "AnalyticFunctions", "BondPrInterp",
             "EqFinDiff", "FullInfoEq",
             "ModelPlots", "JointEq"]
    for modl in modls
        include(string(joinpath(module_path, modl), ".jl"))
    end
end    


println("###############################################################################")
println(" FULL INFORMATION PLOTS #######################################################")
println("###############################################################################")
println(" ")
println(" FI - FIRM VALUE  #############################################################")
plot_title = ModelPlots.get_contour_plot_title(fidf, :firm_value)
file_path_name = ModelPlots.get_contour_plot_path_name(fidf, :fv)
fi_fv_fig = ModelPlots.plot_iso_contour_curves(fi_fd, fi_funs[:fv];
                                               fig_title=plot_title,
                                               file_path_name=file_path_name)
println(" ")
println(" FI - Market-to-Book Ratio  ###################################################")
plot_title = ModelPlots.get_contour_plot_title(fidf, :MBR)
file_path_name = ModelPlots.get_contour_plot_path_name(fidf, :mbr)
fi_mbr_fig = ModelPlots.plot_iso_contour_curves(fi_fd, fi_funs[:mbr];
                                                fig_title=plot_title,
                                                file_path_name=file_path_name)
println("###############################################################################")
println(" ")
println(" ")
println("###############################################################################")
println(" MISREPRESENTATION PLOTS ######################################################")
println("###############################################################################")
println(" ")
println(" MP - Risky Type's MBR   ######################################################")
plot_title = ModelPlots.get_contour_plot_title(misrepdf, :MBR, diff_fun=true)
file_path_name = ModelPlots.get_contour_plot_path_name(misrepdf, :mp_fi_mbr_diff;
                                                       firm_type=:risky)
mp_fi_mbr_diff_fig = ModelPlots.plot_iso_contour_curves(mp_fd, mp_funs[:mp_fi_mbr_diff];
                                                        rmp_diff_fun=mp_funs[:mp_fi_rm_diff],
                                                        fig_title=plot_title,
                                                        file_path_name=file_path_name,
                                                        iso_cmap=ModelPlots.iso_cmaps["misrep"])
println("###############################################################################")
println(" ")
println(" ")
println("###############################################################################")
println(" POOLING EQ. PLOTS ############################################################")
println("###############################################################################")
println(" ")
println(" Pool - Risky Type's MBR ######################################################")
plot_title = ModelPlots.get_contour_plot_title(pooldf, pool_funs, :MBR,
                                               ft=:risky)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :mbr; firm_type=:risky)
pool_r_mbr_fig = ModelPlots.plot_iso_contour_curves(pool_fd, pool_funs[:risky][:mbr];
                                                    rmp_diff_fun=pool_funs[:risky][:jeq_fi_rm_diff],
                                                    fig_title=plot_title,
                                                    file_path_name=file_path_name,
                                                    iso_cmap=ModelPlots.iso_cmaps["pooling"],
                                                    rmp_cmap="Greys", rmp_alpha=.3)
println(" ")
println(" Pool - Risky Type's MBR Differential #########################################")
plot_title = ModelPlots.get_contour_plot_title(pooldf, pool_funs, :MBR;
                                               ft=:risky, diff_fun=true)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :jeq_fi_mbr_diff;
                                                       firm_type=:risky)
pool_fi_r_mbr_diff_fig = ModelPlots.plot_iso_contour_curves(pool_fd,
                                                           pool_funs[:risky][:jeq_fi_mbr_diff];                                                          rmp_diff_fun=pool_funs[:risky][:jeq_fi_rm_diff],
                                                            fig_title=plot_title,
                                                            file_path_name=file_path_name,
                                                            iso_cmap=ModelPlots.iso_cmaps["pooling"],
                                                            rmp_cmap="Greys", rmp_alpha=.3)
println(" ")
println(" Pool - Safe Type's Firm Value  ################################################")
plot_title = ModelPlots.get_contour_plot_title(pooldf, pool_funs, :firm_value;
                                               ft=:safe, s_fi_fv=s_fi_fv)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :fv;
                                                       firm_type=:safe)
pool_s_fv_fig = ModelPlots.plot_iso_contour_curves(pool_fd, pool_funs[:safe][:fv];
                                                   fig_title=plot_title,
                                                   file_path_name=file_path_name,
                                                   iso_cmap=ModelPlots.iso_cmaps["pooling"])
println("###############################################################################")
println(" ")
println(" ")
println("###############################################################################")
println(" SEPARATING EQ. PLOTS #########################################################")
println("###############################################################################")
println(" ")
println(" Sep - Risky Type's MBR #######################################################")
plot_title = ModelPlots.get_contour_plot_title(sepdf, sep_funs, :MBR;
                                               ft=:risky)
file_path_name = ModelPlots.get_contour_plot_path_name(sepdf, :mbr; firm_type=:risky)
sep_r_mbr_fig = ModelPlots.plot_iso_contour_curves(sep_fd, sep_funs[:risky][:mbr];
                                                   fig_title=plot_title,
                                                   file_path_name=file_path_name,
                                                   iso_cmap=ModelPlots.iso_cmaps["separating"])
println(" ")
println(" Sep - Risky Type's MBR Differential ##########################################")
plot_title = ModelPlots.get_contour_plot_title(sepdf, sep_funs, :MBR;
                                               ft=:risky, diff_fun=true)
file_path_name = ModelPlots.get_contour_plot_path_name(sepdf, :jeq_fi_mbr_diff;
                                                       firm_type=:risky)
sep_fi_r_mbr_diff_fig = ModelPlots.plot_iso_contour_curves(sep_fd,
                                                           sep_funs[:risky][:jeq_fi_mbr_diff];
                                                           fig_title=plot_title,
                                                           file_path_name=file_path_name,
                                                           iso_cmap=ModelPlots.iso_cmaps["separating"])
println(" ")
println(" Sep - Safe Type's Firm Value  #################################################")
plot_title = ModelPlots.get_contour_plot_title(sepdf, sep_funs, :firm_value;
                                               ft=:safe, s_fi_fv=s_fi_fv)
file_path_name = ModelPlots.get_contour_plot_path_name(sepdf, :fv; firm_type=:safe)
sep_s_fv_fig = ModelPlots.plot_iso_contour_curves(sep_fd, sep_funs[:safe][:fv];
                                                  fig_title=plot_title,
                                                  file_path_name=file_path_name,
                                                  iso_cmap=ModelPlots.iso_cmaps["separating"])
println("###############################################################################")
println(" ")
println(" ")
println("###############################################################################")
println(" JOINT EQUILIBRIA PLOTS #######################################################")
println("###############################################################################")
println(" ")
println(" EP Market - Risky Type's MBR #################################################")
# fig, ax1, ax2 = ModelPlots.plot_iso_curves(pltd[:X], pltd[:Y], pltd[:r_MBR];
#                                             iso_cmap=ModelPlots.iso_cmaps["pooling"])
# ax1.contourf(pltd[:X], pltd[:Y], pltd[:bool_Z], cmap="GnBu", alpha=.4)
plot_title = ModelPlots.get_final_contour_plot_title(pooldf, :MBR, :risky)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :mbr; firm_type=:risky, 
                                                        fname_eq_type="ep_market")
ep_r_mbr_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[:r_MBR], 
                                                             pltd[:bool_Z];
                                                             fig_title=plot_title,
                                                             file_path_name=file_path_name,
                                                             iso_cmap=ModelPlots.iso_cmaps["pooling"], 
                                                             iso_levels=15, 
                                                             cat_cmap="GnBu")
println(" ")
println(" EP Market - Safe Type's Firm Value  ###########################################")
plot_title = ModelPlots.get_final_contour_plot_title(pooldf, :firm_value, :safe)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :fv; firm_type=:safe, 
                                                        fname_eq_type="ep_market")
ep_s_fv_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[:s_FV], 
                                                            pltd[:bool_Z];
                                                            fig_title=plot_title,
                                                            file_path_name=file_path_name,
                                                            iso_cmap=ModelPlots.iso_cmaps["pooling"], 
                                                            iso_levels=15, 
                                                            cat_cmap="gnuplot")
println(" ")
println(" Dual Market - Safe Type's Firm Value  #########################################")
plot_title = ModelPlots.get_final_contour_plot_title(pooldf, :firm_value, :safe; k_otc=32.5)
file_path_name = ModelPlots.get_contour_plot_path_name(pooldf, :fv; firm_type=:safe, 
                                                        fname_eq_type="dual_market")
ep_otc_s_fv_fig = ModelPlots.plot_equilibria_iso_contour_curves(pltd[:X], pltd[:Y], pltd[:s_FV], 
                                                                pltd[:bool_OTC_EP];
                                                                fig_title=plot_title,
                                                                file_path_name=file_path_name,
                                                                iso_cmap=ModelPlots.iso_cmaps["pooling"], 
                                                                iso_levels=15,
                                                                cat_cmap="viridis_r")
println("###############################################################################")

