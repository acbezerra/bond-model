

# Function below was adapted to run in Julia 1.1.0. It requires the packages
# 1. LinearAlgebra
# 2. DSP

#Polynomial smoothing with the Savitsky Golay filters
#
# Sources
# ---------
# Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf
# Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay
function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)
    #Some error checking
    @assert isodd(windowSize) "Window size must be an odd integer."
    @assert polyOrder < windowSize "Polynomial order must me less than window size."

    halfWindow = Int((windowSize-1)/2)

    #Setup the S matrix of basis vectors. 
    S = zeros(windowSize, polyOrder+1)
    for ct = 0:polyOrder
	#S[:,ct+1] = [-halfWindow:halfWindow].^(ct)
        S[:,ct+1] = range(-halfWindow, stop=halfWindow).^(ct)
    end

    #Compute the filter coefficients for all orders
    #From the scipy code it seems pinv(S) and taking rows should be enough
    G = S*LinearAlgebra.pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv)

    #Pad the signal with the endpoints and convolve with filter 
    # paddedX = [x[1]*ones(halfWindow), x, x[end]*ones(halfWindow)]
    paddedX = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
    # y = conv(filterCoeffs[end:-1:1], paddedX)
    y = DSP.conv(reverse(filterCoeffs), paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]

end


function findlocalmaxima(signal::Vector)
    inds = Int[]
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end
        
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(inds,i)
            end
        end
        
        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end
    
    return inds
end


function filter_k_struct(df; interp_polyOrder::Integer=3,
                             filter_windowSize::Integer=5 * 10^2 + 1, 
                         filter_polyOrder::Integer=3,
                         interp_vars=[:p, :opt_vb, :cvml_vb, :cvmh_vb, :debt, :equity])

    # Interpolate and Filter p, VB, Debt, Equity and Firm Value
    sgF = Dict()
    sgF[:cgrid] = range(minimum(df[:c]), stop=maximum(df[:c]), length=10^4)
    for x in interp_vars 
        tmp = Dierckx.Spline1D(df[:c], df[x], k=interp_polyOrder, bc="extrapolate")
        sgF[x] = savitsky_golay(tmp(sgF[:cgrid]),  filter_windowSize, filter_polyOrder)
    end
    # sgF[:firm_value] = sgF[:debt] .+ sgF[:equity]
    
    return sgF
end
