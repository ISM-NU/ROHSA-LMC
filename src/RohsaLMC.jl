module RohsaLMC

# ROHSA Gaussian decomposition products, velocity-layer cubes and disk
# de-rotation for the LMC H I cubes. Generic FITS/VGT helpers come from the VGT
# package (`Astronomy`, deps/VGT).

using Astronomy
using Astronomy: write_fits, read_fits_image, normalize_fits_extension, robust_negative_noise,
                 velocity_axis_kms
using FFTW
using FITSIO
using Images
using LinearAlgebra: mul!
using Random
using Statistics

include("rohsa.jl")
include("rohsa_layers.jl")
include("derotate.jl")

export HI_COLUMN_DENSITY_FACTOR, DOPPLER_TEMPERATURE_FACTOR
export read_rohsa_parameters, reconstruct_rohsa_components
export doppler_temperature_k, fwhm_doppler_temperature, rohsa_component_rows, rohsa_fit_metrics
export write_rohsa_decomposition_products
export ROHSA_LAYER_NAMES, ROHSA_PHASE_NAMES, velocity_layer_weights, reference_centroid_channels
export build_layer_cubes, gaussian_phase_weights, build_phase_cubes
export derotation_offset_grid, derotate_cube, derotated_header, noise_sigma_map, fill_padding!
export sanitize_reference, clamp_reference

end
