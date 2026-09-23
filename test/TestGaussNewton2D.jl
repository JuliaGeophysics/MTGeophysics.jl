using Test, LinearAlgebra, Random

@testset "2D skin-depth mesh" begin
    f = [0.1, 1.0, 10.0]
    δ = mt2d_skin_depth(100.0, 0.1)
    @test δ ≈ 503.29*sqrt(100/0.1) rtol=1e-4
    layers = mt2d_skin_depth_layers(f; background_resistivity=100.0, z_core_cell=250.0)
    core = findall(==(250.0), layers)
    @test core == 1:length(core)
    @test sum(layers[core]) >= δ && sum(layers[core]) < δ + 250
    @test sum(layers) >= 4δ
    @test all(diff(layers[length(core):end]) .> 0)
    mesh = BuildMesh2D(frequencies=f)
    @test mesh.z_cell_sizes[mesh.n_air_cells+1:end] ≈ mt2d_skin_depth_layers(f)
end

@testset "2D Fréchet derivatives and Gauss-Newton" begin
    mesh = BuildMesh2D(frequencies=[0.3, 3.0], y_core_range=(-500.,500.),
        y_core_cell=250., y_padding=500., air_cells=2, air_top=-1000.,
        ground_layers=[100.,300.,900.], receiver_positions=[-375.,125.,400.])
    rho = build_mt2d_halfspace_model(mesh)
    rng = MersenneTwister(127)
    rho[3:end,:] .*= exp.(0.3randn(rng, size(rho,1)-2, size(rho,2)))
    direction = randn(rng, size(rho))
    fields = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)
    h = 1e-4

    @testset "Directional derivatives, all response fields and modes" begin
        for mode in (:TE, :TM, :TETM)
            tangent = ApplyFrechet2D(mesh, rho, direction; mode, parameterization=:log_resistivity)
            rp = run_mt2d_forward(mesh, rho .* exp.(h*direction); mode)
            rm = run_mt2d_forward(mesh, rho .* exp.(-h*direction); mode)
            for key in fields
                fd = (getproperty(rp,key) - getproperty(rm,key))/(2h)
                @test getproperty(tangent,key) ≈ fd rtol=2e-5 atol=1e-10
            end
            weights = (rho_xy=randn(rng,2,3), phase_xy=randn(rng,2,3),
                z_xy=randn(rng,ComplexF64,2,3), rho_yx=randn(rng,2,3),
                phase_yx=randn(rng,2,3), z_yx=randn(rng,ComplexF64,2,3))
            gradient = ApplyFrechetTranspose2D(mesh, rho, weights; mode, parameterization=:log_resistivity)
            @test dot(gradient,direction) ≈ sum(real(dot(getproperty(weights,k),getproperty(tangent,k))) for k in fields) rtol=1e-8
            @test all(iszero,gradient[1:mesh.n_air_cells,:])
        end
    end

    @testset "Fréchet ordering, boundaries, parameterizations and fixed air" begin
        cells = [CartesianIndex(3,1), CartesianIndex(5,size(rho,2)),
                 CartesianIndex(3,4), CartesianIndex(4,5), CartesianIndex(1,2)]
        for parameterization in (:resistivity,:log_resistivity,:log10_resistivity)
            sens = FrechetDerivative2D(mesh,rho;active_cells=cells,parameterization)
            @test sens.cells == cells
            @test size(sens.z_xy) == (6,5)
            for (j,cell) in enumerate(cells)
                plus, minus = copy(rho), copy(rho)
                if parameterization == :resistivity
                    plus[cell] += h; minus[cell] -= h
                elseif parameterization == :log_resistivity
                    plus[cell] *= exp(h); minus[cell] *= exp(-h)
                else
                    plus[cell] *= 10.0^h; minus[cell] *= 10.0^-h
                end
                rp,rm = run_mt2d_forward(mesh,plus),run_mt2d_forward(mesh,minus)
                for key in (:z_xy,:z_yx)
                    fd = vec(getproperty(rp,key)-getproperty(rm,key))/(2h)
                    @test getproperty(sens,key)[:,j] ≈ fd rtol=2e-4 atol=1e-10
                end
            end
            @test all(iszero,sens.z_xy[:,end])
            @test all(iszero,sens.z_yx[:,end])
        end
        air = zeros(size(rho)); air[1:2,:] .= 1
        @test all(k -> all(iszero,getproperty(ApplyFrechet2D(mesh,rho,air),k)),fields)
        @test all(iszero,ApplyFrechetTranspose2D(mesh,rho,(;)))
        @test_throws ArgumentError run_mt2d_forward(mesh,rho;mode=:bad)
        @test_throws ArgumentError ApplyFrechet2D(mesh,rho,direction;parameterization=:bad)
        @test_throws DimensionMismatch ApplyFrechet2D(mesh,rho,zeros(2,2))
        @test_throws DimensionMismatch FrechetDerivative2D(mesh,rho;active_cells=falses(2,2))
        bad = copy(rho); bad[3,2] = -1
        @test_throws ArgumentError run_mt2d_forward(mesh,bad)
    end

    @testset "Synthetic recovery, objective descent, bounds and fixed cells" begin
        initial = build_mt2d_halfspace_model(mesh)
        truth = copy(initial); truth[3:4,4:5] .= 40
        observed = data_from_response2d(run_mt2d_forward(mesh,truth);impedance_error_fraction=0.03)
        mask = falses(size(rho)); mask[3:4,4:5] .= true
        opts = Inv2DOptions(beta=0.001,max_iter=10,target_rms=0.01,verbose=false)
        result = GaussNewton2D(mesh,initial,observed;options=opts,active_cells=mask)
        @test result.converged
        @test result.reason == :target_rms
        @test result.fit.rms < 0.01
        @test result.fit.rms < result.history[1].rms/100
        @test maximum(abs.(result.resistivity[mask] ./ truth[mask] .- 1)) < 0.03
        @test result.resistivity[.!mask] == initial[.!mask]
        @test all(diff([h.objective for h in result.history]) .< 0)
        @test all(1 .<= result.resistivity[mask] .<= 1e5)
        @test initial == build_mt2d_halfspace_model(mesh)
        @test observed.z_xy == run_mt2d_forward(mesh,truth).z_xy
        predicted = data_from_response2d(result.response)
        fit = chi2_rms2d(observed,predicted)
        @test fit.chi2 ≈ result.fit.chi2
        @test fit.count == result.fit.count == 24

        for mode in (:TE,:TM)
            partial = GaussNewton2D(mesh,initial,observed;active_cells=mask,
                options=Inv2DOptions(mode=mode,beta=0.001,max_iter=5,target_rms=0.,verbose=false))
            @test partial.fit.count == 12
            @test partial.fit.rms < partial.history[1].rms/10
        end
        regularized = GaussNewton2D(mesh,initial,observed;active_cells=mask,
            options=Inv2DOptions(beta=1e5,max_iter=5,target_rms=0.,verbose=false))
        @test norm(log10.(regularized.resistivity[mask] ./ initial[mask])) <
              norm(log10.(result.resistivity[mask] ./ initial[mask]))/10
        @test all(diff([h.objective for h in regularized.history]) .< 0)
        @test regularized.resistivity[.!mask] == initial[.!mask]
        bounded = GaussNewton2D(mesh,initial,observed;active_cells=mask,
            options=Inv2DOptions(beta=0.,max_iter=5,target_rms=0.,log_bounds=(log10(70.),3.),verbose=false))
        @test minimum(bounded.resistivity[mask]) >= 70-1e-10
        @test bounded.fit.rms < bounded.history[1].rms

        masked = deepcopy(observed)
        masked.z_xy[1,1] = complex(NaN,NaN)
        masked.z_yx_error[2,2] = 0
        zeroiter = Inv2DOptions(max_iter=0,verbose=false)
        masked_result = GaussNewton2D(mesh,initial,masked;active_cells=mask,options=zeroiter)
        @test masked_result.fit.count == 20
        @test length(masked_result.history) == 1
        @test masked_result.reason == :max_iter
        @test !masked_result.converged
        @test_throws ArgumentError GaussNewton2D(mesh,initial,observed;active_cells=falses(size(mask)))
        @test_throws ArgumentError GaussNewton2D(mesh,initial,observed;active_cells=[CartesianIndex(1,1)])
        @test_throws ArgumentError GaussNewton2D(mesh,initial,observed;active_cells=[CartesianIndex(3,1),CartesianIndex(3,1)])
        @test_throws ArgumentError GaussNewton2D(mesh,initial,observed;options=Inv2DOptions(beta=-1))
        @test_throws ArgumentError GaussNewton2D(mesh,initial,observed;config=GaussNewton2DConfig(damping=0))
        invalid = deepcopy(observed); invalid.z_xy_error .= NaN; invalid.z_yx_error .= NaN
        @test_throws ArgumentError GaussNewton2D(mesh,initial,invalid)

        @testset "Generic driver and shared derivatives" begin
            generic = Invert2D(mesh,initial,observed;algorithm=GaussNewton2DConfig(),options=opts,active_cells=mask)
            @test generic.resistivity == result.resistivity
            @test result isa GaussNewton2DResult
            # adjoint gradient (one solve per frequency/mode) matches the explicit J'r
            problem_rho = copy(initial); problem_rho[mask] .*= 1.3
            probe = Invert2D(mesh,problem_rho,observed;options=Inv2DOptions(max_iter=0,beta=0.1,verbose=false),active_cells=mask)
            rho0 = copy(problem_rho); rho0[1:mesh.n_air_cells,:] .= 1e9
            opt = Inv2DOptions(beta=0.1)
            R = MTGeophysics._inv2d_regularizer(mesh,rho0,opt)
            cells = findall(mask)
            problem = MTGeophysics.Inv2DProblem(mesh,MTGeophysics._inv2d_data(mesh,observed,:TETM),cells,
                R,R[:,LinearIndices(rho0)[cells]],log10.(rho0),opt)
            state = MTGeophysics._inv2d_evaluate(problem,rho0,log10.(rho0[cells]))
            @test state.objective ≈ probe.history[1].objective
            J = inv2d_frechet(problem,state;method=:forward)
            @test inv2d_frechet(problem,state;method=:adjoint) ≈ J rtol=1e-8 atol=1e-12
            @test inv2d_gradient(problem,state) ≈ inv2d_gradient(problem,state,J) rtol=1e-8
        end

        @testset "NLCG" begin
            nlcg = NLCG2D(mesh,initial,observed;active_cells=mask,
                options=Inv2DOptions(beta=0.001,max_iter=40,target_rms=0.05,verbose=false))
            @test nlcg isa NLCG2DResult
            @test all(diff([h.objective for h in nlcg.history]) .< 0)
            @test nlcg.fit.rms < nlcg.history[1].rms/20
            @test nlcg.resistivity[.!mask] == initial[.!mask]
            @test haskey(nlcg.history[1], :restart)
            @test_throws ArgumentError NLCG2D(mesh,initial,observed;config=NLCG2DConfig(restart=0))
        end
    end
end
