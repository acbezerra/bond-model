
mutable struct JointKStruct
    mu_s::Float64
    mu_b::Float64
    m::Float64
    c::Float64
    p::Float64
    vbl::Float64

    # Safe Firm
    fi_sf_vb::Float64
    sf_vb::Float64

    # Risk Firm
    fi_rf_vb::Float64
    rf_vb::Float64
end


@with_kw mutable struct FirmObj
    bt
    svm
end

mutable struct JointFirms
    jks::JointKStruct
    sf::Firm
    rf::Firm
    cvm_bt::BatchStruct
    svm_bt::BatchStruct
    cvmdf::DataFrame
    svmdf::DataFrame
end


mutable struct JointFDParams
    # Safe Firm
    sf_eq_vbl::Float64
    sf_eq_max::Float64

    # Risk Firm
    rf_eq_vbl::Float64
    rf_eq_max::Float64
end






# #############################################################
mutable struct FirmSpecificParams
    iota::Float64
    lambda::Float64
    sigmah::Float64
end


mutable struct FirmCommonParams
    V0::Float64
    alpha::Float64
    pi::Float64
    r::Float64
    gross_delta::Float64
    xi::Float64
    sigmal::Float64
end


mutable struct JointEqParams
    # Measure of Bonds
    mu_s::Float64
    
    # Transaction Costs
    kep::Float64
    kotc::Float64

    # Safe Firm Params
    sfp::FirmSpecificParams

    # Risky Firm Params
    rfp::FirmSpecificParams

    # Firm Common Params
    fcp::FirmCommonParams
end


mutable struct EPStruct
    ep_ks::JointKStruct
    sf::Firm
    rf::Firm
    sfdf::DataFrame
    rfdf::DataFrame
    jedf::DataFrame
end


mutable struct OTCStruct
    sf::Firm
    rf::Firm
    cvm_bt::BatchStruct
    svm_bt::BatchStruct
    cvmdf::DataFrame
    svmdf::DataFrame
end


mutable struct JointEquilibrium
    # Joint Equilibrium Params
    jep::JointEqParams

    # Electronic Platform Data
    ep::EPStruct

    # OTC Markets Data
    otc::OTCStruct
end
