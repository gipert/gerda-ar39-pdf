#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <cmath>

#ifndef _GERDA_AR39_PDF_HPP_
#define _GERDA_AR39_PDF_HPP_

namespace gerda {

    struct range_fmt {
        double min;
        double step;
        double max;
    };

    // Declare here the parameters of the 3D grid (must match what is written by
    // the compute.jl program)
    constexpr range_fmt energy_kev_fmt = {0., 0.1, 565.};
    constexpr range_fmt fccd_mm_fmt = {0.65, 0.05, 2.4};
    constexpr range_fmt dlf_fmt = {0., 0.1, 1.};
    // plus the allowed channel values
    constexpr range_fmt channel_fmt = {0, 1, 41};

    constexpr size_t n_energy = (energy_kev_fmt.max-energy_kev_fmt.min)/energy_kev_fmt.step+1;
    constexpr size_t n_fccd = (fccd_mm_fmt.max-fccd_mm_fmt.min)/fccd_mm_fmt.step+1;
    constexpr size_t n_dlf = (dlf_fmt.max-dlf_fmt.min)/dlf_fmt.step+1;

    /* Ar-39 pdf lookup table
     *
     * [*][…][…] -> energy (keV)
     * […][*][…] -> full charge collection depth (mm)
     * […][…][*] -> dead layer fraction
     *
     * The chosen data structure is a heap-array of stack-matrices (one for each
     * energy point)
     */

    using model_matrix = double[n_fccd][n_dlf];
    using lookup_table = model_matrix*;

    // germanium channel -> lookup_table map
    std::map<int, lookup_table> lookup_table_map;

    /* load_lookup_table(map, channel)
     *
     * Read Ar39 lookup table for a germanium channel from disk and store it inside `map`
     */
    void load_lookup_table(std::map<int, lookup_table>& table_map, int channel) {

        if (channel < channel_fmt.min or channel > channel_fmt.max) {
            throw std::domain_error("channel value " + std::to_string(channel) + " is out of the allowed range");
        }

        table_map.emplace(channel, new model_matrix[n_energy]);

        auto location = std::string(__FILE__);
        location = location.substr(0, location.find_last_of('/'));
        std::string filename = location + "/../lookup/ar39-pdf-ch" + std::to_string(channel) + ".dat";
        std::ifstream fstr(filename);
        if (!fstr.is_open()) throw std::runtime_error(filename + " not found");

        for (size_t ene_i = 0; ene_i < n_energy; ++ene_i) {
            for (size_t fccd_i = 0; fccd_i < n_fccd; ++fccd_i) {
                for (size_t dlf_i = 0; dlf_i < n_dlf; ++dlf_i) {
                    fstr >> table_map[channel][ene_i][fccd_i][dlf_i];
                }
            }
        }
    }

    /* lower_edge_idx(x, frame)
     *
     * determine `x` closest lower point index in the unidimensional grid defined by `frame`
     */
    size_t lower_edge_idx(const double point, const range_fmt& frame) {
        return size_t((point-frame.min)/frame.step);
    }

    /* lower_edge(x, frame)
     *
     * determine `x` closest lower point in the unidimensional grid defined by `frame`
     */
    double lower_edge(const double point, const range_fmt& frame) {
        return lower_edge_idx(point, frame) * frame.step + frame.min;
    }

    /* ar39_pdf(channel, energy, fccd, dlf)
     *
     * pdf for Ar39 decays as detected in GERDA PhaseII+ for each channel, energy
     * (keV), full charge-collection depth (fccd) and dead-layer fraction (dlf).
     * The initial pdf grid, estracted from Monte Carlo simulations (MaGe framework)
     * is read from disk. A tri-linear interpolation is used to determine pdf values
     * outside the grid.
     */
    double ar39_pdf(int channel, double energy_kev, double fccd_mm, double dlf, bool debug=false) {

        if (debug) std::cout << "gerda::ar39_pdf(" << channel << ", " << energy_kev << ", "
                             << fccd_mm << ", " <<  dlf << ") = " << std::flush;

        if (energy_kev < energy_kev_fmt.min or energy_kev > energy_kev_fmt.max) {
            throw std::domain_error("energy value out of allowed range");
        }

        if (fccd_mm < fccd_mm_fmt.min or fccd_mm > fccd_mm_fmt.max) {
            throw std::domain_error("FCCD value out of allowed range");
        }

        if (dlf < dlf_fmt.min or dlf > dlf_fmt.max) {
            throw std::domain_error("DLF value out of allowed range");
        }

        // load lookup table, if not done yet
        if (lookup_table_map.find(channel) == lookup_table_map.end()) {
            load_lookup_table(lookup_table_map, channel);
        }

        const auto& table = lookup_table_map[channel];

        // perform tri-linear interpolation
        // ref: https://en.wikipedia.org/wiki/Trilinear_interpolation

        auto e_i = lower_edge_idx(energy_kev, energy_kev_fmt);
        auto e_0 = lower_edge(energy_kev, energy_kev_fmt);
        auto e_d = (energy_kev - e_0) / energy_kev_fmt.step;

        auto f_i = lower_edge_idx(fccd_mm, fccd_mm_fmt);
        auto f_0 = lower_edge(fccd_mm, fccd_mm_fmt);
        auto f_d = (fccd_mm - f_0) / fccd_mm_fmt.step;

        auto d_i = lower_edge_idx(dlf, dlf_fmt);
        auto d_0 = lower_edge(dlf, dlf_fmt);
        auto d_d = (dlf - d_0) / dlf_fmt.step;

        auto c_000 = table[e_i  ][f_i  ][d_i  ];
        auto c_100 = table[e_i+1][f_i  ][d_i  ];
        auto c_001 = table[e_i  ][f_i  ][d_i+1];
        auto c_101 = table[e_i+1][f_i  ][d_i+1];
        auto c_010 = table[e_i  ][f_i+1][d_i  ];
        auto c_110 = table[e_i+1][f_i+1][d_i  ];
        auto c_011 = table[e_i  ][f_i+1][d_i+1];
        auto c_111 = table[e_i+1][f_i+1][d_i+1];

        auto c_00 = c_000*(1-e_d)+c_100*e_d;
        auto c_01 = c_001*(1-e_d)+c_101*e_d;
        auto c_10 = c_010*(1-e_d)+c_110*e_d;
        auto c_11 = c_011*(1-e_d)+c_111*e_d;

        auto c_0 = c_00*(1-f_d)+c_10*f_d;
        auto c_1 = c_01*(1-f_d)+c_11*f_d;

        if (debug) std::cout << c_0*(1-d_d)+c_1*d_d << std::endl;

        return c_0*(1-d_d)+c_1*d_d;
    }
}

#endif
