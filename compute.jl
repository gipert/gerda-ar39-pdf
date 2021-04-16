using StatsBase
using Random
using LinearAlgebra
using Interpolations
using KernelDensity
using UpROOT
using DelimitedFiles
using Printf
using ProgressMeter

import KernelDensity: kde
import Interpolations: interpolate

"""
    reshape(h, range)

Change histogram format (min, max, bin size).
"""
function reshape(h::Histogram{<:Any,1}, format::AbstractRange)
    fit(Histogram,
        StatsBase.midpoints(h.edges[1]),
        Weights(h.weights),
        format)
end

"""
    kde(h)

Compute kernel density estimate of histogram.
"""
function kde(h::Histogram{<:Any,1})
    nbins = length(h.weights)
    xmin = first(h.edges[1])
    xmax = last(h.edges[1])

    # use bin contents to reproduce counts (uniform distribution in bin)
    data = Vector{Float64}()
    for i in 1:nbins
        for _ in 1:Int64(h.weights[i])
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
function interpolate(h::Histogram{<:Any,1}; binsize=1)
    hr = normalize(reshape(h, first(h.edges[1]):binsize:last(h.edges[1])))
    itp = interpolate(hr.weights, BSpline(Cubic(Line(OnCell()))))
    return extrapolate(scale(itp, StatsBase.midpoints(hr.edges[1])), 0)
end

"Default Ar39 pdf estimation method"
estimate_ar39_pdf(h::Histogram{<:Any,1}) = interpolate(h, binsize=15)

"""
    read_root_histogram(file, hist)

Open ROOT file `file`, get histogram object `hist` and convert it to StatsBase.Histogram.
"""
function read_root_histogram(filename::String, histname::String)

    @debug "Reading '$histname' from '$filename'..."

    # the following I unfortunately need to hardcode, UpROOT.jl does not
    # support Histograms yet
    hist_format = 5:0.1:650

    if (!isfile(filename))
        @warn "'$filename' does not exist"
        hist = StatsBase.Histogram(0:step(hist_format):last(hist_format), zeros(6500))
        return reshape(hist, hist_format)
    end

    # get histogram weights from ROOT file
    file = UpROOT.TFile(filename)
    # discard under/over-flow bins, make sure bin contents are integer
    weights = convert(Vector{Int64}, file[histname][2:end-1])

    # this is not native in UpROOT.jl yet
    file.pyobj.close()

    # create Julia histogram out of it, use floats so one can e.g. normalize it afterwards
    hist = StatsBase.Histogram(0:step(hist_format):last(hist_format), convert(Vector{Float64}, weights))
    return reshape(hist, hist_format)
end

"""
    get_ar39_histogram(channel, fccd_mm, dlf)

Get Ar39 histogram from gerda-pdfs store.
"""
function get_ar39_histogram(channel::Int, fccd_mm::Float64, dlf::Float64)
    # convert numeric input to string
    fccd_str = Int(fccd_mm*1e3)
    dlf_str = lpad(string(Int(dlf*1e2)), 3, '0')
    return read_root_histogram("histos/nplus-fccd$(fccd_str)um-dlf$(dlf_str)/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root", "raw/M1_ch$channel")
end

"""
    estimate_ar39_pdf(channel, fccd_mm, dlf)

Returns the estimation of the Ar39 pdf for a germanium channel with given FCCD and DLF.
"""
estimate_ar39_pdf(channel::Int, fccd_mm::Float64, dlf::Float64) = estimate_ar39_pdf(get_ar39_histogram(channel, fccd_mm, dlf))

energy_kev_grid = 0:0.1:565
fccd_mm_grid = 0.65:0.05:2.4
dlf_grid = 0:0.1:1

"""
    estimate_ar39_pdf(channel)

Compute and get Ar39 pdfs for all (fccd × dlf) values present in store via
broadcasting [`estimate_ar39_pdf`](@ref).
"""
estimate_ar39_pdf(channel::Int) = estimate_ar39_pdf.(channel, fccd_mm_grid', dlf_grid)

"""
    serialize(channel, pdf_matrix)

Write to disk (one single file) a matrix (fccd × dlf) of Ar39 pdf estimates.

See also: [`estimate_ar39_pdf`](@ref)
"""
function serialize(channel::Int, container#=::Matrix{Function}=#)
    open("lookup/ar39-pdf-ch$channel.dat", "w") do file
        for energy in energy_kev_grid
            for pdf in container
                p = pdf(energy)
                @printf file "%.5g\n" (!isfinite(p) || p < 0) ? 0 : p
            end
        end
    end
end

function run_pipeline()
    @showprogress for ch in [0:5; 7; 9:35; 37:40]
        serialize(ch, estimate_ar39_pdf(ch))
    end
end
