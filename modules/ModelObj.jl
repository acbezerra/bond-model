module ModelObj
	using Parameters
	
	@with_kw struct FirmParams
	    m::Float32
	    alpha::Float32
	    pi::Float32
	
	    r::Float32
	    gross_delta::Float32
	    iota::Float32
	
	    xi::Float32
	    kappa::Float32
	
	    lambda::Float32
	    sigmal::Float32
	    sigmah::Float32
	end
	
	@with_kw mutable struct Firm
	    V0::Float32
	    c::Float32
	    p::Float32
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
end
