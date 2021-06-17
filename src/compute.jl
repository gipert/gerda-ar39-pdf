using StatsBase
using Random
using LinearAlgebra
using Interpolations
using KernelDensity
using BayesianBlocks
using UpROOT
using DelimitedFiles
using ProgressMeter
using Printf
using Plots

import KernelDensity: kde
import Interpolations: interpolate

n_ar39_primaries = 299999999998

energy_lower_thr_keV = 45

energy_kev_grid = 45:20:565
fccd_mm_grid = 0.65:0.05:2.4
dlf_grid = 0:0.1:1

include("disk.jl")
include("nonparametric.jl")
# include("parametric.jl")

"Default Ar39 pdf estimation method"
estimate_ar39_pdf(h::Histogram{<:Real,1}) = interpolate(h, binsize=20)

"""
    estimate_ar39_pdf(channel, fccd_mm, dlf)

Returns the estimation of the Ar39 pdf for a germanium channel with given FCCD and DLF.
"""
estimate_ar39_pdf(channel::Int, fccd_mm::Float64, dlf::Float64) = estimate_ar39_pdf(get_ar39_histogram(channel, fccd_mm, dlf))

"""
    estimate_ar39_pdf(channel)

Compute and get Ar39 pdfs for all (fccd × dlf) values present in store via
broadcasting [`estimate_ar39_pdf`](@ref).
"""
estimate_ar39_pdf(channel::Int) = estimate_ar39_pdf.(channel, fccd_mm_grid', dlf_grid)

function run_pipeline()
    @showprogress for ch in [0:5; 7; 9:35; 37:40]
        serialize(ch, estimate_ar39_pdf(ch))
    end
end

function display_ar39_pdf(channel::Int, fccd_mm::Float64, dlf::Float64)

    hist = get_ar39_histogram(channel, fccd_mm, dlf)
    pdf = estimate_ar39_pdf(channel, fccd_mm, dlf)

    plot(size=(700,400),
         ylim=(0,Inf), xlim=(0,250),
         title=@sprintf("Channel %i / FCCD = %.2f mm / DLF = %.2f", channel, fccd_mm, dlf),
         xlabel="Energy (keV)",
         ylabel="Counts / decay / keV")

    plot!(hist, st=:step, label="0.1 keV MC hist")
    plot!(pdf.itp.ranges[1], [pdf(e) for e in pdf.itp.ranges[1]], st=:scatter, color=:red, label="Interpolation nodes")
    plot!(0:0.1:565, [pdf(e) for e in 0:0.1:565], linewidth=2, color=:red, label="Cubic spline")
end
