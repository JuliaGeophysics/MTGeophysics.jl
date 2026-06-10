# Core/padding detection regression tests (headless, no solver).
# Author: @pankajkmishra
# Guards against the narrow-core misdetection bug: when the uniform fine-cell
# survey area is a minority of cells along an axis, core_indices must still
# return the full uniform plateau, not a sliver at the padding boundary.

@testset "Core detection (CoreUtils3D)" begin

    # Build cell centers from a width vector, mirroring load_ws3d_model.
    centers(dw) = begin
        edges = vcat(0.0, cumsum(dw))
        (edges[1:end-1] .+ edges[2:end]) ./ 2
    end

    # Geometrically-expanding padding either side of a uniform fine core.
    function widths(; npad_lo, ncore, npad_hi, core_w, factor)
        lo = reverse([core_w * factor^k for k in 1:npad_lo])
        hi =          [core_w * factor^k for k in 1:npad_hi]
        vcat(lo, fill(core_w, ncore), hi)
    end

    @testset "narrow core, symmetric padding" begin
        # Core is only 14 of 50 cells in x (the failing geometry).
        dx = widths(npad_lo = 18, ncore = 14, npad_hi = 18, core_w = 2062.0, factor = 1.3)
        ix = core_indices(centers(dx))
        @test ix == 19:32
        @test length(ix) == 14
    end

    @testset "wide core" begin
        dy = widths(npad_lo = 18, ncore = 53, npad_hi = 18, core_w = 2062.0, factor = 1.3)
        iy = core_indices(centers(dy))
        @test iy == 19:71
        @test length(iy) == 53
    end

    @testset "offset (asymmetric) padding" begin
        # Different padding counts each side; core must still be found exactly.
        dz = widths(npad_lo = 6, ncore = 20, npad_hi = 22, core_w = 500.0, factor = 1.4)
        iz = core_indices(centers(dz))
        @test iz == 7:26
        @test length(iz) == 20
    end
end
