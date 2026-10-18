using Astronomy
using Statistics
using Test

@testset "De-rotation" begin
    nv = 60
    gauss(mu, s) = [exp(-((k - 1 - mu)^2) / (2s^2)) for k in 1:nv]
    cube = zeros(3, 1, nv)
    cube[1, 1, :] = gauss(20.0, 3.0)     # line at the reference velocity
    cube[2, 1, :] = gauss(30.5, 3.0)
    cube[3, 1, :] .= NaN                  # no data
    reference = reshape([20.0, 30.5, 25.0], 3, 1)
    out, omin = Astronomy.derotate_cube(cube, reference)
    offsets = omin .+ (0:size(out, 3) - 1)

    # Integer shift: an exact copy of the input, placed at offset = channel - reference.
    j0 = findfirst(==(-20), offsets)
    @test out[1, 1, j0:j0 + nv - 1] == Float32.(cube[1, 1, :])
    # Outside the observed band: NaN, so VGT sees no artificial edge there.
    @test all(isnan, out[1, 1, 1:j0 - 1]) && all(isnan, out[1, 1, j0 + nv:end])
    # Fractional shift: flux conserved and the line lands at offset 0.
    @test sum(v -> isfinite(v) ? Float64(v) : 0.0, out[2, 1, :]) ≈ sum(cube[2, 1, :]) rtol = 1e-5
    @test offsets[argmax(replace(out[2, 1, :], NaN32 => -Inf32))] == 0
    @test sum(v -> isfinite(v) ? abs2(v) : 0.0, out[2, 1, :]) ≈ sum(abs2, cube[2, 1, :]) rtol = 1e-4  # all-pass
    @test all(isnan, out[3, 1, :])

    # A wild reference far outside the band (no-emission pixels) does not blow up the grid.
    wild = reshape([20.0, 30.5, 9000.0], 3, 1)
    @test size(Astronomy.derotate_cube(cube, wild)[1], 3) <= 2nv - 1

    # Integer mode: whole-channel copy even for a fractional reference.
    outi, ominI = Astronomy.derotate_cube(cube, reference; integer_shift=true)
    ji = findfirst(==(-30), ominI .+ (0:size(outi, 3) - 1))   # round(30.5) == 30 (round-to-even)
    @test outi[2, 1, ji:ji + nv - 1] == Float32.(cube[2, 1, :])

    # Padding fill: zero for model cubes, noise of the local level for data; no-data pixels stay NaN.
    observed = reshape([true, true, false], 3, 1)
    z = Astronomy.fill_padding!(copy(out), observed)
    @test all(z[1, 1, 1:j0 - 1] .== 0) && all(isnan, z[3, 1, :])
    big = fill(NaN32, 64, 64, 4)
    noisy = Astronomy.fill_padding!(big, trues(64, 64); sigma=fill(0.25, 64, 64))
    @test all(isfinite, noisy)
    @test isapprox(std(Float64.(vec(noisy))), 0.25; rtol=0.1)
    @test minimum(Astronomy.fill_padding!(fill(NaN32, 64, 64, 2), trues(64, 64); sigma=fill(0.25, 64, 64), clip=true)) == 0

    # Reference sanitizing: trusted where the model has emission, constant where it has none.
    refmap = fill(50.0, 40, 40); refmap[:, 31:40] .= 5000.0          # junk where there is no emission
    colmap = fill(100.0, 40, 40); colmap[:, 31:40] .= 0.0
    clean, wt = Astronomy.sanitize_reference(refmap, colmap, 159; smooth_px=1.0)
    @test clean[20, 10] ≈ 50 && clean[20, 40] ≈ 50 && maximum(clean) <= 158
    @test wt[20, 10] ≈ 1 && wt[20, 40] ≈ 0

    # Header: spatial WCS copied unchanged, and the whole header parses in WCSLIB.
    h = Dict{String, Any}("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 975.0,
                          "CRPIX2" => 1087.0, "CRVAL1" => 81.19, "CRVAL2" => -68.53,
                          "CDELT1" => -0.00555555, "CDELT2" => 0.00555555, "EQUINOX" => 2000.0,
                          "CTYPE3" => "VELO-LSR", "CUNIT3" => "m/s", "CRPIX3" => 1.0,
                          "CRVAL3" => 181689.5, "CDELT3" => 1000.0)
    hd = Astronomy.derotated_header(h, 317, -158)
    for key in ("CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "EQUINOX")
        @test hd[key] == h[key]
    end
    @test hd["CRVAL3"] == -158 * 1000.0
    @test length(Astronomy.WCS.from_header(Astronomy._wcs_header_string(hd))) == 1
end
