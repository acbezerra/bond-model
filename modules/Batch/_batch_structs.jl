

@with_kw mutable struct DirFileNames
    # Directories
    main_dir::String
    res_dir::String
    mat_dir_prefix::String
    coupon_dir_prefix::String

    # File Names
    debt_at_par_cp_fn_prefix::String
    eq_fd_cp_fn_prefix::String
end


@with_kw mutable struct ModelInputs
    _svm_dict::Dict{Symbol, Float64}

    # Directories
    batch_res_dir::String
    maturity_dir::String
    comb_res_dir::String

    # Paths
    main_dir_path::String
    batch_res_path::String
    maturity_path::String
    comb_res_path::String
end


@with_kw mutable struct BatchDFColumns
    main_params::Array{Symbol,1}
    k_struct_params::Array{Symbol,1}
    fixed_params::Array{Symbol,1}
    debt_vars::Array{Symbol,1}
    equity_vars::Array{Symbol,1}
    share_values::Array{Symbol,1}
    dfcols::Array{Symbol,1}
end


@with_kw mutable struct BatchParams
    _param_values_dict::Dict{Symbol,Array{T,1} where T}
    _common_params::Dict{Symbol,Float64}
    _params_order::Array{Symbol,1}
    _params_combs::Array{Array{Float64,1},9}
    df::DataFrame
end



@with_kw mutable struct BatchStruct
    # Model Type
    model::String
    
    # Paths & File Names
    dfn::DirFileNames

    # DataFrame Columns
    dfc::BatchDFColumns

    # Model Object Inputs
    mi::ModelInputs

    # Parameters
    bp::BatchParams

    # Coupon
    # coupon_grid::StepRangeLen{Float64,Base.TwicePrecision{Float64},
    #                           Base.TwicePrecision{Float64}}
    coupon_grid::Array{Float64,1}
    pvb_tol_vec::Array{Float64,1}
end

