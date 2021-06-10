using StatsBase
using Random
using LinearAlgebra
using Interpolations
using KernelDensity
using BayesianBlocks
using UpROOT
using DelimitedFiles
using ProgressMeter
using EmpiricalDistributions
using BAT
using LsqFit
using Distributions
using IntervalSets
using ValueShapes
using Printf

import KernelDensity: kde
import Interpolations: interpolate

n_ar39_primaries = 299999999998

energy_lower_thr_keV = 30

energy_kev_grid = 0:0.1:565
fccd_mm_grid = 0.65:0.05:2.4
dlf_grid = 0:0.1:1

include("disk.jl")
include("nonparametric.jl")
# include("parametric.jl")

"Default Ar39 pdf estimation method"
estimate_ar39_pdf(h::Histogram{<:Real,1}) = interpolate(h, binsize=optimal_ar39_bw(h))

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
