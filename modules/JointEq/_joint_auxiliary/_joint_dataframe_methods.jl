

# Match DataFrames' Rows ###########################################################
function identify_match(df::DataFrame, dfrow::DataFrame)
    compf = var -> abs.(df[:, var] .- dfrow[1, var]) .< 1e-6
    fspf = var -> isnan.(dfrow[1, var]) ? isnan.(df[:, var]) : compf(var)
    
    fcp_cond  = [compf(var) for var in commoncols]
    fsp_cond = [fspf(var) for var in fspcols]
    
    return .&([.&(x...) for x in [fcp_cond, fsp_cond]]...)
end

    
function identify_matching_row(df::DataFrame, svm)
    compf = var -> abs.(df[:, var] .- get_param(svm, var)) .< 1e-6
    fspf = var -> isnan.(get_param(svm, var)) ? isnan.(df[:, var]) : compf(var)
    
    fcp_cond  = [compf(var) for var in commoncols]
    fsp_cond = [fspf(var) for var in fspcols]
    
    return .&([.&(x...) for x in [fcp_cond, fsp_cond]]...)
end


function identify_match2(df::DataFrame, dfrow::DataFrame)
    compf = var -> abs.(df[:, var] .- dfrow[1, var]) .< 1e-6
    scompf = (var, prefix) -> (abs.(df[:, Symbol(prefix, var)] .- 
                               dfrow[1, Symbol(prefix, var)]) .< 1e-6)
    fspf = (var, prefix) -> (isnan.(dfrow[1, Symbol(prefix, var)]) ? 
                             isnan.(df[:, Symbol(prefix, var)]) : scompf(var, prefix))

    # Check Parameters
    fcp_cond  = [compf(var) for var in vcat(:kappa, commoncols)]
    sfsp_cond = [fspf(var, :s_) for var in fspcols]
    rfsp_cond = [fspf(var, :r_) for var in fspcols]
    epm_cond = [compf(var) for var in epmcols if var != :kappa]
    
    return .&([.&(x...) for x in [fcp_cond, sfsp_cond, rfsp_cond, epm_cond]]...)
end


function identify_matching_row2(df::DataFrame, jks, jf)
    compf = (svm, var) -> abs.(df[:, var] .- get_param(svm, var)) .< 1e-6
    scompf = (svm, var, prefix) -> abs.(df[:, Symbol(prefix, var)] .- get_param(svm, var)) .< 1e-6
    fspf = (svm, var, prefix) -> isnan.(get_param(svm, var)) ? isnan.(df[:, Symbol(prefix, var)]) : scompf(svm, var, prefix)

    # Parameter Comparison Function
    pcf = var -> abs.(df[:, var] .- getfield(jks, var)) .< 1e-6

    # Check Parameters
    fcp_cond  = [compf(jf.sf, var) for var in vcat(:kappa, commoncols)]
    sfsp_cond = [fspf(jf.sf, var, :s_) for var in fspcols]
    rfsp_cond = [fspf(jf.rf, var, :r_) for var in fspcols]
    epm_cond = [pcf(var) for var in epmcols if var != :kappa]
    
    return .&([.&(x...) for x in [fcp_cond, sfsp_cond, rfsp_cond, epm_cond]]...)
end
# ##################################################################################


# Reshape Misrep, Pool and Separating Eq DFs #######################################
function reshape_sf_rf_df(df::DataFrame)
    specific_cols = [col for col in names(df) if 
                     !(col in vcat(:sf_defaults_first,
                                   epmcols,
                                   commoncols, vbcols))]

    # Capital Structure
    ksdf = DataFrame(df[1, [:sf_defaults_first, epmcols..., :vb]])

    # Common Parameters
    commondf =  DataFrame(df[1, commoncols])

    # Safe Firm Results
    sfdf = DataFrame(df[1, vcat([:sf_vb], specific_cols)])
    names!(sfdf, vcat([:sf_vb], [Symbol(:s_, col) for col in specific_cols
                                 if (col != :sf_defaults_first)]))

    # Risky Firm Results
    rfdf = DataFrame(df[2, vcat([:rf_vb], specific_cols)])
    names!(rfdf, vcat([:rf_vb], [Symbol(:r_, col) for col in specific_cols
                                 if (col != :sf_defaults_first)]))

    return hcat(ksdf, sfdf, rfdf, commondf)
end
# ###############################################################################


# Slice DataFrames ##############################################################
function slice_df_cond(df::DataFrame, svm, rerun::Bool)
    row_cond = false
    if !isempty(df)
        row_cond = identify_matching_row(df, svm)
    end

    if (sum(row_cond) == 0)
        println("No matches found!")
        rerun = true
    elseif .&(sum(row_cond) == 1, !rerun)
        println("Match found!")
    end
        
    return row_cond, rerun
end


function slice_full_info_dataframe(df::DataFrame, jks, svm; rerun::Bool=true)
    row_cond, rerun = slice_df_cond(df, svm, rerun)

    if rerun
        println("Generating results...")
        dfrow = find_optimal_bond_measure(svm; jks=jks) 
    else
        println("Extracting row...")
        dfrow = DataFrame(df[row_cond, :])
    end
        
    return dfrow
end


function slice_mps_dataframe(df::DataFrame, jks, jf; rerun::Bool=true)
    row_cond = false
    if !isempty(df)
        row_cond = identify_matching_row2(df, jks, jf)
    end

    if (sum(row_cond) == 0)
        println("No matches found!")
        rerun = true
    elseif .&(sum(row_cond) == 1, !rerun)
        println("Match found!")
    end

    if rerun
        println("Generating results...")
        dfrow = find_joint_optimal_vb(jf, jks;
                                      rerun_fi_vb=true)
    else
        println("Extracting row...")
        dfrow = DataFrame(df[row_cond, :])
    end
       
    return dfrow
end
# ###############################################################################


# Save DataFrames ###############################################################
function update_file(fpath::String, fname::String, df::DataFrame)
    if !JointEq.exists_ep_df(fpath, fname)
        CSV.write(string(fpath, "/", fname, ".csv"), df)
    else
        if fname == fidf_name
            resdf = CSV.read(string(fpath, "/", fname, ".csv"); types=fidf_col_types)
        else
            resdf = CSV.read(string(fpath, "/", fname, ".csv"); types=mps_col_types)
        end
        
        for i in 1:nrow(df)
            # check if there is a match
            if fname == fidf_name
                row_cond = identify_match(resdf, DataFrame(df[i, :]))
            else
                row_cond = identify_match2(resdf, DataFrame(df[i, :]))
            end
            
            if sum(row_cond) == 1
                println("Match found! Replacing row...")
                resdf[row_cond, :] = DataFrame(df[i, :])
            elseif sum(row_cond) == 0
                println("No match found. Appending row...")
                resdf = vcat(resdf, DataFrame(df[i, :]))
            else
                println(string("Multiple matches for row ", i, 
                               "! Please refine ID columns. Exiting..."))
                return
            end
        end
        
        # Save File
        println(string("Saving File ", fname," ..."))
        CSV.write(string(fpath, "/", fname, ".csv"), resdf)
    end
end
# ###############################################################################