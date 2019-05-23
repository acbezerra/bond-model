using Distributions
using Interpolations
using Distributed
using DataFrames
using Printf
using Dierckx
using JLD
using CSV

start_tic = time_ns()

main_path = "/home/artur/BondPricing"
module_path = string(main_path, "/", "Julia/modules/")
include(string(module_path, "/", "TestFunctions.jl"))
modls = ["Batch", "ModelObj", "AnalyticFunctions", 
         "BondPrInterp", "EqFinDiff"]
for modl in modls
    include(string(joinpath(module_path, modl), "/", modl, ".jl"))
end


println(string("ARGUMENTS: ", ARGS))
# ################ SYS ARGUMENTS ################
# Capture Combination Number:
comb_num = parse(Int, ARGS[1])

# Skip Bond Surfaces Computation?
skip_bpr_computation = false
if size(ARGS, 1) > 1
    skip_bpr = parse(Int, ARGS[2])
    skip_bpr_computation = Bool(skip_bpr * (skip_bpr in [0,1]))
end
println(string("Skip bond pricing surfaces computation: ", skip_bpr_computation))

# ##########  COUPON VALUES ##########
# Get Coupon Value and set Coupon Counter:
# Position in Coupon Grid
c_pos = 999
if size(ARGS, 1) > 2
    c_pos = parse(Int, ARGS[3])
end

# Set Coupon Counter:
c_counter = 0

# Skip calculations if values already exist:
skip_c = false
if size(ARGS, 1) > 3
    skip_c_input = parse(Int, ARGS[4])
    skip_c = Bool(skip_c_input * (skip_c_input in [0, 1]))
end
println(string("Skip coupon calculations: ", skip_c))

# Skip Batch Part I
skip_sol = false
if size(ARGS, 1) > 4
    skip_sol_input = parse(Int, ARGS[5])
    skip_sol = Bool(skip_sol_input * (skip_sol_input in [0, 1]))
end
println(string("Skip Debt at Par calculations: ", skip_sol))

# Skip Finite Differences Method
skip_all_eqdf = false
if size(ARGS, 1) > 5
    skip_all_eqdf_input = parse(Int, ARGS[6])
    skip_all_eqdf = Bool(skip_all_eqdf_input * (skip_all_eqdf_input in [0, 1]))
end
println(string("Skip Equity Finite Differences Calculations: ", skip_all_eqdf))
# ###############################################

# ## Create Batch and SVM objects ########################
# if bond pricing surfaces do not exist, the function will
# compute and save them in the appropriate directory.
bt, svm = Batch.get_bt_svm(;comb_num=comb_num,
                           display_msgs=true,
                           compute_surfs=!skip_bpr_computation)


# ################ CREATE DIRECTORIES ################
println("Managing results directories...")
# Set Directories & Paths (Main, Batch, Maturity, Combination)
# bt = Batch.set_comb_res_paths(bt)

# Create Directories
bt = Batch.mk_comb_res_dirs(bt)

if c_pos > length(bt.coupon_grid)
    println("Computations will be done for all coupon values!")
    c_pos = 1
else
    println(string("Coupon position: ", c_pos,
                   ", value: ", bt.coupon_grid[c_pos]))
    c_counter = length(bt.coupon_grid) - 1
end

# If results already exist and skip_c is true,
# set counter to 999 to skip calculations
eqfiles = [x for x in readdir(bt.mi.comb_res_path) if
           occursin(bt.dfn.eq_fd_cp_fn_prefix, x)]
all_files = size(eqfiles, 1) == size(bt.coupon_grid, 1)
if all_files & skip_c
    println("All results have been computed. Skipping calculations...")
    c_counter = 999
end

# Check coupon condition:
while c_counter < size(bt.coupon_grid, 1)
    # Set Coupon
    svm.c = bt.coupon_grid[c_pos]

    println(" =============================================================")
    println(" =============================================================")
    println(string(" ===================== ",
                   "Coupon Value: ", svm.c,
                   " ====================="))
    println(" =============================================================")
    println(" =============================================================")

    # Start Timer
    tic = time_ns()

    # Coupon Files
    debt_at_par_cp_fname = string(bt.dfn.debt_at_par_cp_fn_prefix,
                                  @sprintf("%.2f", svm.c), ".csv")
    debt_at_par_cp_path_fname = joinpath(bt.mi.comb_res_path,  debt_at_par_cp_fname)
    eq_fd_cp_fname = string(bt.dfn.eq_fd_cp_fn_prefix, @sprintf("%.2f", svm.c), ".csv")
    all_eq_fd_cp_fname = string("all_",eq_fd_cp_fname)
    eq_fd_cp_path_fname = joinpath(bt.mi.comb_res_path,  eq_fd_cp_fname)
    all_eq_fd_cp_path_fname = joinpath(bt.mi.comb_res_path,  all_eq_fd_cp_fname)
    
    # If coupon folder does not contain final results or
    # if the coupon skip condition is false,
    # run the calculations below:
    if !isfile(eq_fd_cp_path_fname) | (skip_c == false)
        # If first part solution is missing or if skip
        # condition for the 1st part of the calculations
        # is false, run the 1st part:
        if !isfile(debt_at_par_cp_path_fname) | (skip_sol == false)
            global skip_all_eqdf=false
            # ############ Form VB and p grids: ############
            # 1. Form P Grid 
            # 2. Compute CVM VB value for sigmal and sigmah
            # 3. Interpolate VB functions
            # 4. Compute VB s.t. Debt is issued at par
            # 5. Form Adjusted p Grid: [.85 x min(Opt P), 1.15 x max(Opt P)]
            count = Batch.p_candidates(svm; c=svm.c)

            # Bound VB Values: 
            # For each p in Adjusted p-grid, find min(VB) and max(VB) s.t. .6 <= vbhl <= 1.4 
            # that is, ratio of pre and post default barriers must be whitin bounds.
            # Now, for each p in adjusted p-grid, I have an interval of VB candidate values
            # for which debt is issued at par.
            BVB = [vcat(p, [i for i in Batch.bounded_vbs(svm, count, p)])
                   for p in count["adj_pgrid"]]

            # Store results in a DataFrame
            df = DataFrame(transpose(hcat(BVB...)), [:p, :vblmin, :vblmax])

            println("Printing p VB dataframe...")
            println(df)

            # ################ Debt at Par #################
            # For each (p, VB) candidate, compute
            # i. Debt
            # ii. Debt - P
            # iii. abs(Debt - P)
            # iv. abs(Debt -P)/P
            # v. (Debt - P)/P
            LL = @time fetch(@spawn [Batch.svm_debt_at_par(svm,
                                                           df[i, :vblmin],
                                                           df[i, :vblmax],
                                                           svm.c; p=df[i, :p]) 
                                     for i in range(1, stop=size(df, 1))])

            # Store results in a DataFrame
            DFI = vcat(LL...)

            # Save Results
            CSV.write(debt_at_par_cp_path_fname, DFI)
            
            print(string("Time taken to compute first part: ", time_ns() - tic))

            println(" ===================================================")
            println(" =============== First part is done! ===============")
            println(" ===================================================")
            println(" ")
        end
        
        # ###########################################################################
        # ########################## Compute Equity Values ##########################
        # ###########################################################################
        println(" ===================================================")
        println(" ============= Starting second part... =============")
        println(" ===================================================")
        println(" ")

        # Start Timer
        tic2 = time_ns()
        
        println(debt_at_par_cp_path_fname)
        
        if !skip_all_eqdf | !(all_eq_fd_cp_fname in readdir(bt.mi.comb_res_path))
            # Load DFI 
            df = Batch.conditionally_load_df(:DFI, debt_at_par_cp_path_fname)

            # Filter DataFrame
            fdf = Batch.filter_batch_I_df(bt, svm, df; tol=.05)
            
            # DataFrame Size:
            println(string("Part I DF length: ", size(df, 1)))
            println(string("Part I Filtered DF length: ", size(fdf, 1)))

            # Apply Finite Differences Method
            all_eqdf = Batch.eq_fd_method(bt, svm, fdf)

            # Save All Equity Finite Differences Results
            CSV.write(all_eq_fd_cp_path_fname, all_eqdf)
        else
            all_eqdf = CSV.read(all_eq_fd_cp_path_fname; 
                                types=vcat(fill(Float64, 16), [Bool], fill(Float64, 14)))
        end
        
        # Compute candidate vb for each p value
        eqdf_final = Batch.eq_fd_processing(bt, svm, all_eqdf)

        # Compute (p, vb) solution for c value
        res = Batch.eq_fd_sol(bt, svm, eqdf_final)

        # Save results
        Batch.save_eq_results(bt, res, eqdf_final, eq_fd_cp_path_fname)
    end

    global c_pos += 1
    global c_counter += 1
end

println(string("Total Script Run Time: ", (time_ns() - start_tic)/1e9/60., " minute(s)."))
    

    
    
