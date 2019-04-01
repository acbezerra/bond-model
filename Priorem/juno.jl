push!(LOAD_PATH, "/home/artur/BondPricing/Julia/modules/")
include(string(LOAD_PATH[end], "AnalyticFunctions.jl"))

include(string(LOAD_PATH[end], "AnalyticFunctions.jl"))
# reload("AnalyticFunctions")

# nprocs()
# nworkers()

pardict = Dict("V0" => 100,
               "vbl" => 60,
               "vbh" => 62,
                "m" => 1.,
                "c" => 6.39,
                "p" => 61.68,
                "sigmal" => .25,
                "sigmah" => .25,
                "r" => .08,
                "gross_delta" => .02,
                "iota" => .0,
                "xi" => 1.,
                "kappa" => .015,
                "alpha" => .6,
                "pi" => .27,
                "lambda" => .3)

using Parameters

@with_kw struct FirmParams
    m::Float16
    alpha::Float16
    pi::Float16

    r::Float16
    gross_delta::Float16
    iota::Float16

    xi::UInt8
    kappa::Float16

    lambda::Float16
    sigmal::Float16
    sigmah::Float16
end

@with_kw mutable struct Firm
    V0::Float16
    c::Float16
    p::Float16
    vbl::Float32
    vbh::Float32

    pm::FirmParams
end


function firm_constructor(dict)

    return Firm(dict["V0"], dict["c"], dict["p"],
                   dict["vbl"], dict["vbh"],
                FirmParams(dict["m"], dict["alpha"], dict["pi"],
                dict["r"], dict["gross_delta"], dict["iota"],
                dict["xi"], dict["kappa"],
                dict["lambda"], dict["sigmal"], dict["sigmah"]))
end

arm = firm_constructor(pardict)

AnalyticFunctions


tic()

u = .5
ttm = arm.pm.m
vmax=2
vbl = 55
vt = log(arm.V0/vbl)
v = .8
#vbh = AnalyticFunctions.zhi_vb(arm.rm.m, c, p, sigma, r, gross_delta, iota, xi, k, alpha, pi)
vbh = 60

tauN = 10^3
vN = 10^3
dt = (ttm - 1e-4)/tauN
dv = vmax/vN

ugrid = linspace(1e-4, ttm, tauN)
vgrid = linspace(0, vmax, vN)

Psi = AnalyticFunctions.bond_int_Psi(vt, v,
                                u, ttm, # vmax,
                                arm.pm.sigmal,
                                arm.pm.r, arm.pm.gross_delta,
                                arm.pm.xi, arm.pm.kappa, arm.pm.lambda,
                                ugrid, vgrid) * dt * dv

FPsi = AnalyticFunctions.bond_int_FPsi(vt, v, vbl, vbh,
                                        u, ttm, # vmax,
                                        arm.pm.sigmal, arm.pm.sigmah,
                                        arm.pm.r, arm.pm.gross_delta,
                                        arm.pm.xi, arm.pm.kappa, arm.pm.lambda,
                                        ugrid, vgrid) * dt * dv

GPsi= AnalyticFunctions.bond_int_GPsi(vt, v, vbl, vbh,
                                        u, ttm, # vmax,
                                        arm.pm.sigmal, arm.pm.sigmah,
                                        arm.pm.r, arm.pm.gross_delta,
                                        arm.pm.xi, arm.pm.kappa, arm.pm.lambda,
                                        ugrid, vgrid) * dt * dv


toc()

println(string("Psi: ", Psi))
println(string("FPsi: ", FPsi))
println(string("GPsi: ", GPsi))
