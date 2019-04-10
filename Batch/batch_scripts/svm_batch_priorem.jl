
# Create Model Object
# bt = Batch.BatchObj()
# bt = Batch.set_par_dict(bt, comb_num)


# # ############################################
# comb_num = 2 
# bt = Batch.set_par_dict(bt, comb_num)

# # Adjust parameters
# bt.mi._svm_dict["lambda"] = .3
# bt.mi._svm_dict["sigmal"] = .2
# bt.mi._svm_dict["sigmah"] = .2
# bt.mi._svm_dict["kappa"] = 0.015

# Set Directories & Paths (Main, Batch, Maturity, Combination)
# bt = Batch.set_comb_res_paths(bt)

# Create Directories
# bt = Batch.mk_comb_res_dirs(bt)
# ############################################

# batch_obj_exists = false
# try
#     global batch_obj_exists = Batch.check_batch_results(bt)
# catch
#     println("Unable to load Batch object file. Recomputing... ")
# end

# println(string("Batch object exists: ", batch_obj_exists))

# # Check if results exist
# if batch_obj_exists
#     # Load SVM object
#     println("Loading SVM object...")
#     svm = Batch.load_batch_results(bt)["svm"]
    
#     # Check if Bond Pricing Surfaces exist
#     surfs_exist = Batch.check_bpr_surfaces(svm)
# else
#     # Create SVM object
#     svm = ModelObj.firm_constructor(bt.mi._svm_dict)
#     surfs_exist = false
# end

# # If Bond Pricing Surfaces do not exist, compute them
# if !surfs_exist | !skip_bpr_computation
#     # Compute Bond Price Inputs
#     svm = @time BondPrInterp.bpr_surfs(svm)
    
#     # Save Results
#     try 
#         Batch.save_batch_results(bt; svm=svm)
#         println("Batch Bond Surfaces saved.")
#     catch 
#         files = [x for x in readdir(bt.mi.comb_res_path)  
#                  if occursin(Batch.batch_file_name, x) | 
#                 occursin("bond_surfaces", x)]
#         if sum(files) > 0
#             [rm(string(bt.mi.comb_res_path, "/", x)) for x in files]
#         end
        
#         try
#             Batch.save_batch_results(bt; svm=svm)
#         catch
#             println("Unable to save batch file!")
#         end
#     end
# end

# # println("returning...")
# # return 

# # Interpolate Bond Pricing Surfaces
# println("Interpolating bond pricing surfaces...")
# svm = @time BondPrInterp.bpr_interp(svm)

