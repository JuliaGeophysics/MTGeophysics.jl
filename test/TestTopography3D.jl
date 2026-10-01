# GeoTIFF elevations must retain their geographic positions when used by the 3D mesh builder.
using ArchGDAL

@testset "3D GeoTIFF topography orientation" begin
    # A sloping surface with unequal x/y gradients exposes transposition even on square rasters.
    elevation(x, y) = 100.0 + 2.0 * x + 3.0 * y
    for (nx, ny) in ((3, 2), (3, 3)), (dx, dy) in ((5.0, -10.0), (-5.0, 10.0))
        @testset "$(nx)x$(ny), pixel spacing ($dx, $dy)" begin
            mktempdir() do dir
                path = joinpath(dir, "terrain.tif")
                xs = 10.0 .+ ((0:nx-1) .+ 0.5) .* dx
                ys = 100.0 .+ ((0:ny-1) .+ 0.5) .* dy
                # GDAL's first array dimension is the raster width (x).
                pixels = [elevation(x, y) for x in xs, y in ys]
                ArchGDAL.create(path; driver = ArchGDAL.getdriver("GTiff"),
                                width = nx, height = ny, nbands = 1, dtype = Float64) do ds
                    ArchGDAL.setgeotransform!(ds, [10.0, dx, 0.0, 100.0, 0.0, dy])
                    ArchGDAL.write!(ArchGDAL.getband(ds, 1), pixels)
                end

                topo = MTGeophysics._load_topography_geotiff(path)
                @test topo.x == sort(xs)
                @test topo.y == sort(ys)
                @test size(topo.z) == (ny, nx)
                @test topo.bbox == (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
                for x in xs, y in ys
                    @test MTGeophysics._sample_topography(topo, x, y) ≈ elevation(x, y)
                end
                # Check interpolation away from pixel centres as well as exact pixel samples.
                xq = 0.25 * minimum(xs) + 0.75 * maximum(xs)
                yq = 0.75 * minimum(ys) + 0.25 * maximum(ys)
                @test MTGeophysics._sample_topography(topo, xq, yq) ≈ elevation(xq, yq)
            end
        end
    end
end
