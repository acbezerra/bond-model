@with_kw mutable struct KStruct
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64
end


@with_kw mutable struct FirmParams
    V0::Float64
    alpha::Float64
    pi::Float64
    
    r::Float64
    gross_delta::Float64
    iota::Float64
    
    xi::Float64
    kappa::Float64
    
    lambda::Float64
    sigmal::Float64
    sigmah::Float64
end


# Bond Pricing Inputs
@with_kw mutable struct BPrInputs
    # Grid Bounds & Lengths
    vtmax::Float64
    vtN::Int64
    ttm_max::Float64
    ttmN::Int64
    vbhlmin::Float64
    vbhlmax::Float64
    vbhlN::Int64

    vmax::Float64
    vN::Int64
    uN::Int64

    vtN_ref::Int64
    ttmN_ref::Int64
    vbhlN_ref::Int64
end


# Bond Pricing Inputs
@with_kw mutable struct BPrFixedTTMInputs
    # Grid Bounds & Lengths
    vtmax::Float64
    vtN::Int64
    ttm::Float64
    vbhlmin::Float64
    vbhlmax::Float64
    vbhlN::Int64
end


# Bond Pricing Grids & Surfaces
@with_kw mutable struct BPrSurfs
    # Grids
    vtgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    ttmgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    vbhlgrid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    
    # Surfs
    f11_surf::Array{Float64,2}
    f12_surf::Array{Float64,3}
    f13_surf::Array{Float64,3}
    f21_surf::Array{Float64,2}
    f22_surf::Array{Float64,2}
end


# Bond Pricing Interpolated Functions
@with_kw mutable struct BPrInterpFuns
    f11
    f12
    f13
    f21
    f22
end


@with_kw mutable struct Firm
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64
    vbh::Float64
    pm::FirmParams
    model::String
    bi::BPrInputs
    bs::BPrSurfs
    bf::BPrInterpFuns
    bit::BPrFixedTTMInputs
    bft::BPrInterpFuns
    optKS::KStruct
end
