module Gauss 

using Random

function experimental(m; mu=0., sigma=1.) # default parameters from jctc 2019
    y = randn(m).*sigma .+ mu
    return y
end

function theoretical(n, Y; sigma=2., offset=1.) # default parameters from jctc 2019
    y = randn((n, size(Y,1)))
    for alpha in axes(y,1)
        y[alpha,:] .=  y[alpha,:].*sigma .+ Y .+ offset
    end
    return y
end

end # module 
