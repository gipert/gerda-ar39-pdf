"""
    reshape(h, range)

Change histogram format (min, max, bin size).
"""
function reshape(h::Histogram{<:Real,1}, format::AbstractRange; isdensity::Bool=false)
    hr = fit(Histogram,
             StatsBase.midpoints(h.edges[1]),
             Weights(h.weights),
             format)
    isdensity && (hr.weights .*= step(h.edges[1])/step(format))
    hr
end

"""
    kde(h)

Compute kernel density estimate of histogram.
"""
function kde(h::Histogram{<:Real,1})
    nbins = length(h.weights)
    xmin = first(h.edges[1])
    xmax = last(h.edges[1])

    # use bin contents to reproduce counts (uniform distribution in bin)
    data = Vector{Float64}()
    for i in 1:nbins
        for _ in 1:true_ar39_counts(h, i)
            energy = xmin+(i-1+rand())*(xmax-xmin)/nbins
            push!(data, energy)
        end
    end

    return kde(data)
end

"""
    interpolate(h)

Perform cubic spline interpolation of an histogram (can be rebinned first).
"""
function interpolate(h::Histogram{<:Real,1}; binsize=-1)
    # rebin first
    hr = h
    if binsize > 0
        hr = reshape(h, first(h.edges[1]):binsize:last(h.edges[1]), isdensity=true)
    end

    # interpolate
    itp = interpolate(hr.weights, BSpline(Cubic(Line(OnCell()))))
    return extrapolate(Interpolations.scale(itp, StatsBase.midpoints(hr.edges[1])), 0)

    # or with Dierckx.jl
    # return Spline1D(StatsBase.midpoints(hr.edges[1]), hr.weights, k=3, bc="zero")
end

"""
    optimal_ar39_bw(h)
Compute optimal number of bins for an Ar39 histogram. Rice rule tuned by hand.
"""
function optimal_ar39_bw(h::Histogram{<:Real,1})
    nbins = ceil(Integer, 2*cbrt(true_ar39_counts(h)))
    binw = 1.5*(last(h.edges[1]) - first(h.edges[1]))/nbins
    return binw > 13 ? binw : 13
end
