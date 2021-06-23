true_ar39_counts(h::Histogram{<:Any,1}) = ceil(Integer, norm(h)*n_ar39_primaries*step(h.edges[1]))
true_ar39_counts(h::Histogram{<:Any,1}, binindex::Integer) = ceil(Integer, h.weights[binindex]*n_ar39_primaries*step(h.edges[1]))

"""
    read_root_histogram(file, hist)

Open ROOT file `file`, get histogram object `hist` and convert it to StatsBase.Histogram.
"""
function read_root_histogram(filename::String, histname::String)

    @debug "Reading '$histname' from '$filename'..."

    if (!isfile(filename))
        @warn "'$filename' does not exist"
        error()
    end

    file = UpROOT.TFile(filename)
    hist = file[histname]

    # this is not native in UpROOT.jl yet
    file.pyobj.close()

    return hist
end

"""
    get_ar39_histogram(channel, fccd_mm, dlf)

Get Ar39 histogram from gerda-pdfs store. Normalize it as a probability density
scaled by number of Monte Carlo primaries.
"""
function get_ar39_histogram(channel::Int, fccd_mm::Float64, dlf::Float64)
    # convert numeric input to string
    fccd_str = Int(fccd_mm*1e3)
    dlf_str = lpad(string(Int(dlf*1e2)), 3, '0')
    h = read_root_histogram("histos/nplus-fccd$(fccd_str)um-dlf$(dlf_str)/lar/sur_array_4/Ar39/pdf-lar-sur_array_4-Ar39.root", "raw/M1_ch$channel")
    h = reshape(h, energy_lower_thr_keV:step(h.edges[1]):last(h.edges[1]))
    h.weights ./= n_ar39_primaries
    normalize!(h, mode=:density)
    return h
end

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

