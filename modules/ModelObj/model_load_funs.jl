
function load_firm_params_bpr_inputs(svm_input_path::String; file_name::String="_svm_inputs")
    df = CSV.read(string(svm_input_path, "/", file_name, ".csv"))

    for var in obj_params
        if !(var in names(df))
            df[var] = NaN
        end
    end
    
    # ####### Generate Firm Params Dictionary #######
    params_list = vcat(obj_params, firm_params)    
    params = Dict{Symbol, Any}(cn => df[1, cn] for cn in params_list if
                                !occursin(:Column, cn))
    
    # ####### Generate BPr Inputs Dictionary #######
    cols = [:vtN, :ttmN, :vbhlN, :vN, :uN]
    bpr_inputd = Dict{Symbol, Any}(cn => df[1, cn] for cn in names(df) if
                                   !occursin(:Column, cn) &
                                   !(cn in cols) & !(cn in params_list))

    # Grid Length Numbers Should be Interger
    for cn in cols
        bpr_inputd[cn] = trunc(Integer, df[1, cn])
    end

    return  params, bpr_inputd
end


function load_bpr_surfs(svm_input_path::String)
    fsurfs = [string("f", f, "_surf") for f in ["11", "12", "13", "21", "22"]]
    bprf = Dict()
    for f in fsurfs
        bprf[Symbol(f)] = load(string(svm_input_path, "/", f, ".jld"))[f]
    end

    return bprf
end

