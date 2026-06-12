@testset "Shapefile overlay utilities" begin

    @testset "prj sidecar and CRS detection" begin
        mktempdir() do dir
            shp = joinpath(dir, "test.shp")
            @test MTGeophysics.prj_path_from_shp(shp) == joinpath(dir, "test.prj")
            @test detect_shapefile_crs(shp) == (:unknown, "")

            write(joinpath(dir, "test.prj"), "PROJCS[\"ETRS89 / TM35FIN\"]")
            crs_type, wkt = detect_shapefile_crs(shp)
            @test crs_type == :projected
            @test !isempty(wkt)

            write(joinpath(dir, "test.prj"), "GEOGCS[\"WGS 84\"]")
            crs_type, _ = detect_shapefile_crs(shp)
            @test crs_type == :geographic
        end
    end

    @testset "coordinate transform fallbacks" begin
        mktempdir() do dir
            shp = joinpath(dir, "test.shp")

            f, transformed, info = shapefile_coord_transform(shp, "model")
            @test !transformed
            @test f(3.0, 4.0) == (3.0, 4.0)

            f, transformed, info = @test_logs (:warn,) match_mode = :any begin
                shapefile_coord_transform(shp, "EPSG:3067")
            end
            @test !transformed
            @test info == "missing .prj"
            @test f(3.0, 4.0) == (3.0, 4.0)
        end
    end

    @testset "geographic to projected transform" begin
        mktempdir() do dir
            shp = joinpath(dir, "test.shp")
            wgs84 = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\"," *
                    "SPHEROID[\"WGS 84\",6378137,298.257223563]]," *
                    "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
            write(joinpath(dir, "test.prj"), wgs84)

            f, transformed, _ = shapefile_coord_transform(shp, "EPSG:3067")
            @test transformed
            e, n = f(27.0, 65.0)
            @test isapprox(e, 500000.0; atol = 10.0)
            @test n > 7.2e6
        end
    end

    @testset "prepare_shapefiles skips bad entries" begin
        defs = [
            (path = "missing_a.shp", enabled = false, color = :red, alpha = 1.0,
             point_size = 5, line_width = 1.0),
            (path = "missing_b.shp", enabled = true, color = :red, alpha = 1.0,
             point_size = 5, line_width = 1.0),
        ]
        loaded = @test_logs (:warn,) match_mode = :any begin
            prepare_shapefiles(defs, "EPSG:3067")
        end
        @test isempty(loaded)
    end
end
