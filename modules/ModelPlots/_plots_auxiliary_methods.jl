

function set_cvm_data(pt::PlotStruct, df::DataFrame)
    if isempty(pt.cvm_data)
        println("Setting the Data")
    else
        println("Updating the Data")
    end

    pt.cvm_data=df
end


function plots_form_combinations(bt, fig_name_vars::Array{Symbol,1})
    value_lists = [bt.bp._param_values_dict[x] for x in fig_name_vars]

    # Remaining Params
    # rparams = [x for x in bt._params_order if !(x in fig_name_vars)]

    combs  = [[x1, x2, x3, x4, x5, x6] for 
              x1 in bt.bp._param_values_dict[fig_name_vars[1]],
              x2 in bt.bp._param_values_dict[fig_name_vars[2]],
              x3 in bt.bp._param_values_dict[fig_name_vars[3]],
              x4 in bt.bp._param_values_dict[fig_name_vars[4]],
              x5 in bt.bp._param_values_dict[fig_name_vars[5]],
              x6 in bt.bp._param_values_dict[fig_name_vars[6]]]

    # combinations  = hcat(combs...)'

    return combs
end


function plot_svmdf_slicer(pt, fixed_params::Dict{Symbol, Float64})
    locs = [abs.(pt.svm_data[:, x] .- fixed_params[x]) .< 1e-6 for x in keys(fixed_params)]
    srows = sum(hcat(locs...), dims=2) .== size(hcat(locs...), 2)
    return pt.svm_data[vcat(srows...), :]
end


function get_cvm_svm_dfs(cvmdict::Dict{Symbol,Array{Float64,1}},
                         svmdict::Dict{Symbol,Array{Float64,1}};
                         firm_obj_fun::Symbol=:firm_value)
    m = (:m in keys(cvmdict)) ? cvmdict[:m][1] : NaN
    
    bt = BatchObj()
    pt = PlotsObj(bt; firm_obj_fun=firm_obj_fun, cvm_m=m, svm_m=m)
    sbt = get_bt(; model="svm", m=m, m_comb_num=1)
    cbt = get_bt(; model="cvm", m=m, m_comb_num=1)
    
    # Set Batch Objects
    cbt = get_bt(;model="cvm", m=m, m_comb_num=1)
    sbt = get_bt(;model="svm", m=m, m_comb_num=1)

    # Get Combination Numbers
    cvm_combs = get_batch_comb_numbers(cbt, cvmdict)[:comb_num]
    svm_combs = get_batch_comb_numbers(sbt, svmdict)[:comb_num]
    # #########################################################


    # Safe and Risky Firms' DFs ###############################
    cvmdf = pt.cvm_data[[(x in cvm_combs) for x in pt.cvm_data[:comb_num]], :]
    svmdf = pt.svm_data[[(x in svm_combs) for x in pt.svm_data[:comb_num]], :]
    # #########################################################

    # x-axis variable:
    xvar =  [x for x in cvs_xvars if size(unique(svmdf[x]), 1)  == size(unique(pt.svm_data[x]), 1)]
    if !isempty(xvar)
        if size(xvar, 1) > 1
            println("Multiple xvars!")
        end
        xvar = xvar[1]
    end
    
    return cvmdf, svmdf, xvar
end

