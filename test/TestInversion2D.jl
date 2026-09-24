# 2D deterministic inversion
# Author: @pankajkmishra
# Ensures Gauss-Newton and NLCG recover a buried conductor with a falling objective, honour fixed cells, modes,
# regularization and masked data, refuse bad settings, share one driver whose adjoint Jacobian and gradient match
# the explicit ones, and keep topographic air and water out of the model and its regularization

using Test, LinearAlgebra, Random, SparseArrays

@testset "2D deterministic inversion" begin
    mesh = BuildMesh2D(frequencies=[0.3, 3.0], y_core_range=(-500.,500.),
        y_core_cell=250., y_padding=500., air_cells=2, air_top=-1000.,
        ground_layers=[100.,300.,900.], receiver_positions=[-375.,125.,400.])

    @testset "Synthetic recovery, objective descent, bounds and fixed cells" begin
        initial = build_mt2d_halfspace_model(mesh)
        truth = copy(initial); truth[3:4,4:5] .= 40
        observed = data_from_response2d(run_mt2d_forward(mesh,truth);impedance_error_fraction=0.03)
        mask = falses(size(initial)); mask[3:4,4:5] .= true
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

@testset "2D inversion with topographic air and water" begin
    mesh = _hill_mesh()
    air = mt2d_air_mask(mesh)
    rng = MersenneTwister(4)
    ρ = build_mt2d_halfspace_model(mesh)
    ρ[3:end, :] .*= exp.(0.3randn(rng, size(ρ, 1) - 2, size(ρ, 2)))
    water = falses(size(ρ)); water[5:6, 1] .= true
    excluded = air .| water
    R = MTGeophysics._inv2d_regularizer(mesh, ρ, Inv2DOptions(); excluded)
    @test nnz(R[:, findall(vec(excluded))]) == 0
    @test nnz(R[:, findall(vec(.!excluded))]) > 0

    observed = data_from_response2d(run_mt2d_forward(mesh, ρ))
    start = build_mt2d_halfspace_model(mesh)
    result = Invert2D(mesh, start, observed; algorithm = GaussNewton2DConfig(),
                      options = Inv2DOptions(max_iter = 2, target_rms = 0.0, verbose = false), water_cells = water)
    @test all(c -> !air[c] && !water[c], result.active_cells)
    @test result.resistivity[water] == start[water]
    @test result.history[end].rms < result.history[1].rms
end
