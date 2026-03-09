using CairoMakie
using Printf
using Random

"""
    MT1DMesh

Inputs:
- Layer interfaces, layer thicknesses, finite-difference nodes, and finite-difference resistivities.

Output:
- `MT1DMesh`: Container for the layered 1D model and its finite-difference discretization.

Description:
- Stores the 1D earth model in both layered and finite-difference form.
"""
Base.@kwdef struct MT1DMesh
    interfaces::Vector{Float64}
    thicknesses::Vector{Float64}
    fd_nodes::Vector{Float64}
    fd_resistivities::Vector{Float64}
end

"""
    MT1DDataSpec

Inputs:
- Frequencies, apparent-resistivity error fractions, phase errors, and an optional source path.

Output:
- `MT1DDataSpec`: Container for a 1D MT survey specification.

Description:
- Stores the frequency axis and data uncertainties used by the 1D workflows.
"""
Base.@kwdef struct MT1DDataSpec
    frequencies::Vector{Float64}
    rho_error_fraction::Vector{Float64}
    phase_error_deg::Vector{Float64}
    path::String = ""
end

"""
    MT1DResponse

Inputs:
- Solver name, frequencies, impedance values, apparent resistivities, and phases.

Output:
- `MT1DResponse`: Container for a 1D MT forward response.

Description:
- Stores the complex impedance and the derived apparent-resistivity and phase curves.
"""
Base.@kwdef struct MT1DResponse
    method::Symbol
    frequencies::Vector{Float64}
    impedance::Vector{ComplexF64}
    apparent_resistivity::Vector{Float64}
    phase::Vector{Float64}
end

const μ₀_1D = 4π * 1e-7
const DEFAULT_MT1D_LAYER_THICKNESSES = [120.0, 280.0, 650.0, 1400.0]
const DEFAULT_MT1D_LAYER_RESISTIVITIES = [100.0, 20.0, 350.0, 40.0, 800.0]
const DEFAULT_MT1D_FIRST_CELL_M = 20.0
const DEFAULT_MT1D_GROWTH_FACTOR = 1.16
const DEFAULT_MT1D_PADDING_CELLS = 30

"""
    build_mt1d_mesh(t, ρ; first_cell=25.0, growth_factor=1.12, padding_cells=25)

Inputs:
- `t`: Layer thicknesses in metres.
- `ρ`: Layer resistivities in ohm metres, including the basement half-space.
- `first_cell`, `growth_factor`, `padding_cells`: Finite-difference mesh controls.

Output:
- `MT1DMesh`: Layered model and refined finite-difference mesh.

Description:
- Builds the 1D finite-difference mesh used by the analytical and numerical forward solvers.
"""
function build_mt1d_mesh(
    t::AbstractVector{<:Real},
    ρ::AbstractVector{<:Real};
    first_cell::Real = 25.0,
    growth_factor::Real = 1.12,
    padding_cells::Integer = 25,
)
    length(ρ) == length(t) + 1 || error("resistivities must contain one extra basement value")
    any(Δz -> Δz <= 0, t) && error("all layer thicknesses must be positive")
    first_cell > 0 || error("first_cell must be positive")
    growth_factor > 1 || error("growth_factor must be greater than 1")
    padding_cells > 0 || error("padding_cells must be positive")

    Δz = Float64.(t)
    ρv = Float64.(ρ)
    z_interfaces = vcat(0.0, cumsum(Δz))
    z_nodes = Float64[0.0]
    ρfd = Float64[]

    for (i, Δzi) in enumerate(Δz)
        n_cells = max(1, ceil(Int, Δzi / first_cell))
        h = Δzi / n_cells
        for _ in 1:n_cells
            push!(z_nodes, z_nodes[end] + h)
            push!(ρfd, ρv[i])
        end
    end

    h_tail = length(z_nodes) > 1 ? z_nodes[end] - z_nodes[end - 1] : Float64(first_cell)
    for _ in 1:padding_cells
        h_tail *= growth_factor
        push!(z_nodes, z_nodes[end] + h_tail)
        push!(ρfd, ρv[end])
    end

    MT1DMesh(
        interfaces = z_interfaces,
        thicknesses = Δz,
        fd_nodes = z_nodes,
        fd_resistivities = ρfd,
    )
end

"""
    BuildMesh1D(t, ρ; kwargs...)

Inputs:
- `t`: Layer thicknesses in metres.
- `ρ`: Layer resistivities including the basement half-space.

Output:
- `MT1DMesh`: Layered model and finite-difference mesh.

Description:
- Public alias for `build_mt1d_mesh`.
"""
BuildMesh1D(t::AbstractVector{<:Real}, ρ::AbstractVector{<:Real}; kwargs...) =
    build_mt1d_mesh(t, ρ; kwargs...)

"""
    write_mt1d_model(path, mesh; title="MT1D finite-difference model")

Inputs:
- `path`: Output model path.
- `mesh`: 1D mesh to write.
- `title`: Header line written to disk.

Output:
- `String`: Path to the written model file.

Description:
- Writes a 1D finite-difference model file used by the package examples and inversion workflows.
"""
function write_mt1d_model(
    path::AbstractString,
    mesh::MT1DMesh;
    title::AbstractString = "MT1D finite-difference model",
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# $(title)")
        println(io, "NODES $(length(mesh.fd_nodes))")
        for z in mesh.fd_nodes
            println(io, Printf.@sprintf("%.8e", z))
        end
        println(io, "RESISTIVITY $(length(mesh.fd_resistivities))")
        for ρ in mesh.fd_resistivities
            println(io, Printf.@sprintf("%.8e", ρ))
        end
    end
    String(path)
end

"""
    load_mt1d_model(path)

Inputs:
- `path`: Path to a saved 1D model file.

Output:
- `MT1DMesh`: Parsed layered and finite-difference mesh.

Description:
- Reads a 1D model file and reconstructs the layered model implied by the finite-difference cells.
"""
function load_mt1d_model(path::AbstractString)
    isfile(path) || error("1D model file not found: $path")
    lines = readlines(path)

    z_nodes = Float64[]
    ρfd = Float64[]
    mode = :none
    expected = 0

    for raw in lines
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "#") && continue
        if startswith(line, "NODES ")
            expected = parse(Int, split(line)[2])
            empty!(z_nodes)
            mode = :nodes
            continue
        elseif startswith(line, "RESISTIVITY ")
            expected == length(z_nodes) || error("node count mismatch while reading $path")
            expected = parse(Int, split(line)[2])
            empty!(ρfd)
            mode = :resistivity
            continue
        end

        value = parse(Float64, line)
        if mode == :nodes
            push!(z_nodes, value)
        elseif mode == :resistivity
            push!(ρfd, value)
        else
            error("unexpected content in $path: $line")
        end
    end

    length(z_nodes) >= 2 || error("1D model file $path must define at least two nodes")
    length(ρfd) == length(z_nodes) - 1 || error("1D model resistivity count must equal node count minus one")

    Δz_cells = diff(z_nodes)
    z_interfaces = Float64[0.0]
    t = Float64[]
    current_ρ = ρfd[1]
    current_t = Δz_cells[1]

    for i in 2:length(ρfd)
        if isapprox(ρfd[i], current_ρ; rtol = 0.0, atol = 1e-10)
            current_t += Δz_cells[i]
        else
            push!(t, current_t)
            push!(z_interfaces, z_interfaces[end] + current_t)
            current_ρ = ρfd[i]
            current_t = Δz_cells[i]
        end
    end

    MT1DMesh(
        interfaces = z_interfaces,
        thicknesses = t,
        fd_nodes = Float64.(z_nodes),
        fd_resistivities = Float64.(ρfd),
    )
end

"""
    mt1d_layered_model(mesh)

Inputs:
- `mesh`: 1D finite-difference mesh.

Output:
- Named tuple with `thicknesses` and `resistivities`.

Description:
- Collapses equal-resistivity finite-difference cells back into a layered-earth description.
"""
function mt1d_layered_model(mesh::MT1DMesh)
    Δz_cells = diff(mesh.fd_nodes)
    ρlayers = Float64[]
    tlayers = Float64[]
    current_ρ = mesh.fd_resistivities[1]
    current_t = Δz_cells[1]

    for i in 2:length(mesh.fd_resistivities)
        if isapprox(mesh.fd_resistivities[i], current_ρ; rtol = 0.0, atol = 1e-10)
            current_t += Δz_cells[i]
        else
            push!(ρlayers, current_ρ)
            push!(tlayers, current_t)
            current_ρ = mesh.fd_resistivities[i]
            current_t = Δz_cells[i]
        end
    end

    push!(ρlayers, current_ρ)
    t = isempty(tlayers) ? Float64[] : tlayers
    if !isempty(t) && length(ρlayers) != length(t) + 1
        t = t[1:(length(ρlayers) - 1)]
    end
    (thicknesses = t, resistivities = ρlayers)
end

"""
    write_mt1d_data_template(path, f; rho_error_fraction=0.05, phase_error_deg=1.5)

Inputs:
- `path`: Output data-template path.
- `f`: Frequencies in hertz.
- `rho_error_fraction`, `phase_error_deg`: Uncertainty floors.

Output:
- `String`: Path to the written template.

Description:
- Writes a 1D MT data template that contains frequencies and uncertainty columns but no observations.
"""
function write_mt1d_data_template(
    path::AbstractString,
    f::AbstractVector{<:Real};
    rho_error_fraction::Real = 0.05,
    phase_error_deg::Real = 1.5,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "frequency_hz,apparent_resistivity_ohm_m,phase_deg,re_z,im_z,rho_error_fraction,phase_error_deg")
        for fi in Float64.(f)
            Printf.@printf(io, "%.8e,NaN,NaN,NaN,NaN,%.6f,%.6f\n", fi, Float64(rho_error_fraction), Float64(phase_error_deg))
        end
    end
    String(path)
end

"""
    load_mt1d_data_spec(path)

Inputs:
- `path`: Path to a 1D data template or observed-data file.

Output:
- `MT1DDataSpec`: Parsed frequency axis and uncertainty vectors.

Description:
- Reads the uncertainty specification used by the 1D forward and inversion workflows.
"""
function load_mt1d_data_spec(path::AbstractString)
    isfile(path) || error("1D data file not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("1D data file $path is empty")

    f = Float64[]
    σρ = Float64[]
    σϕ = Float64[]

    for raw in lines[2:end]
        line = strip(raw)
        isempty(line) && continue
        tokens = split(line, ",")
        length(tokens) >= 7 || error("1D data row must contain 7 columns in $path")
        push!(f, parse(Float64, tokens[1]))
        push!(σρ, parse(Float64, tokens[6]))
        push!(σϕ, parse(Float64, tokens[7]))
    end

    MT1DDataSpec(
        frequencies = f,
        rho_error_fraction = σρ,
        phase_error_deg = σϕ,
        path = String(path),
    )
end

"""
    write_mt1d_observed_data(path, specification, response)

Inputs:
- `path`: Output observed-data path.
- `specification`: Survey frequencies and uncertainties.
- `response`: 1D MT response sampled on the same frequencies.

Output:
- `String`: Path to the written observed-data file.

Description:
- Writes a complete 1D observed-data file with impedance, apparent resistivity, phase, and errors.
"""
function write_mt1d_observed_data(
    path::AbstractString,
    specification::MT1DDataSpec,
    response::MT1DResponse,
)
    length(specification.frequencies) == length(response.frequencies) || error("specification and response sizes do not match")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "frequency_hz,apparent_resistivity_ohm_m,phase_deg,re_z,im_z,rho_error_fraction,phase_error_deg")
        for i in eachindex(specification.frequencies)
            Printf.@printf(
                io,
                "%.8e,%.8e,%.6f,%.8e,%.8e,%.6f,%.6f\n",
                specification.frequencies[i],
                response.apparent_resistivity[i],
                response.phase[i],
                real(response.impedance[i]),
                imag(response.impedance[i]),
                specification.rho_error_fraction[i],
                specification.phase_error_deg[i],
            )
        end
    end
    String(path)
end

"""
    load_mt1d_observed_data(path)

Inputs:
- `path`: Path to a 1D observed-data file.I sti

Output:
- Named tuple with `specification` and `response`.

Description:
- Reads a 1D observed-data file and returns both the survey specification and observed response.
"""
function load_mt1d_observed_data(path::AbstractString)
    isfile(path) || error("1D observed-data file not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("1D observed-data file $path is empty")

    f = Float64[]
    ρa = Float64[]
    ϕ = Float64[]
    Z = ComplexF64[]
    σρ = Float64[]
    σϕ = Float64[]

    for raw in lines[2:end]
        line = strip(raw)
        isempty(line) && continue
        tokens = split(line, ",")
        length(tokens) >= 7 || error("1D observed-data row must contain 7 columns in $path")
        push!(f, parse(Float64, tokens[1]))
        push!(ρa, parse(Float64, tokens[2]))
        push!(ϕ, parse(Float64, tokens[3]))
        push!(Z, ComplexF64(parse(Float64, tokens[4]), parse(Float64, tokens[5])))
        push!(σρ, parse(Float64, tokens[6]))
        push!(σϕ, parse(Float64, tokens[7]))
    end

    specification = MT1DDataSpec(
        frequencies = f,
        rho_error_fraction = σρ,
        phase_error_deg = σϕ,
        path = String(path),
    )
    response = MT1DResponse(
        method = :observed,
        frequencies = f,
        impedance = Z,
        apparent_resistivity = ρa,
        phase = ϕ,
    )
    (specification = specification, response = response)
end

"""
    build_default_mt1d_mesh(; resistivities=DEFAULT_MT1D_LAYER_RESISTIVITIES)

Inputs:
- Optional layer resistivities including the basement half-space.

Output:
- `MT1DMesh`: Default benchmark mesh populated with the requested resistivities.

Description:
- Builds the standard 1D benchmark discretization used by the examples, benchmarks, and inversion workflow.
"""
function build_default_mt1d_mesh(;
    resistivities::AbstractVector{<:Real} = DEFAULT_MT1D_LAYER_RESISTIVITIES,
)
    build_mt1d_mesh(
        DEFAULT_MT1D_LAYER_THICKNESSES,
        resistivities;
        first_cell = DEFAULT_MT1D_FIRST_CELL_M,
        growth_factor = DEFAULT_MT1D_GROWTH_FACTOR,
        padding_cells = DEFAULT_MT1D_PADDING_CELLS,
    )
end

"""
    MakeMesh1D(; output_dir=..., model_name="Layered1D.true", thicknesses=..., resistivities=..., first_cell=20.0, growth_factor=1.16, padding_cells=30)

Inputs:
- Output location and layered-earth parameters.

Output:
- Named tuple with `mesh` and `model_path`.

Description:
- Builds the default 1D benchmark mesh and writes it to disk.
"""
function MakeMesh1D(;
    output_dir::AbstractString = joinpath(dirname(@__DIR__), "Models"),
    model_name::AbstractString = "Layered1D.true",
    thicknesses::AbstractVector{<:Real} = DEFAULT_MT1D_LAYER_THICKNESSES,
    resistivities::AbstractVector{<:Real} = DEFAULT_MT1D_LAYER_RESISTIVITIES,
    first_cell::Real = DEFAULT_MT1D_FIRST_CELL_M,
    growth_factor::Real = DEFAULT_MT1D_GROWTH_FACTOR,
    padding_cells::Integer = DEFAULT_MT1D_PADDING_CELLS,
)
    mesh = build_mt1d_mesh(
        thicknesses,
        resistivities;
        first_cell = first_cell,
        growth_factor = growth_factor,
        padding_cells = padding_cells,
    )
    model_path = write_mt1d_model(joinpath(output_dir, model_name), mesh)
    (mesh = mesh, model_path = model_path)
end

"""
    _package_mt1d_response(method, f, Z)

Inputs:
- `method`: Solver label.
- `f`: Frequencies in hertz.
- `Z`: Complex impedance values.

Output:
- `MT1DResponse`: Derived MT response.

Description:
- Converts complex impedances into apparent-resistivity and phase curves.
"""
function _package_mt1d_response(
    method::Symbol,
    f::AbstractVector{<:Real},
    Z::Vector{ComplexF64},
)
    fv = Float64.(f)
    ω = 2π .* fv
    MT1DResponse(
        method = method,
        frequencies = fv,
        impedance = Z,
        apparent_resistivity = abs.(Z) .^ 2 ./ (ω .* μ₀_1D),
        phase = rad2deg.(atan.(imag.(Z), real.(Z))),
    )
end

"""
    solve_mt1d_analytical(f, ρ, t)

Inputs:
- `f`: Frequencies in hertz.
- `ρ`: Layer resistivities including the basement half-space.
- `t`: Layer thicknesses in metres.

Output:
- `MT1DResponse`: Analytical 1D MT response.

Description:
- Solves the layered-earth MT problem with the recursive impedance formula.
"""
function solve_mt1d_analytical(
    f::AbstractVector{<:Real},
    ρ::AbstractVector{<:Real},
    t::AbstractVector{<:Real},
)
    fv = Float64.(f)
    ρv = Float64.(ρ)
    tv = Float64.(t)
    σ = 1.0 ./ ρv
    ω = 2π .* fv

    k = sqrt.(-1im .* ω .* μ₀_1D .* σ[end])
    Z = (ω .* μ₀_1D) ./ k

    for i in reverse(eachindex(tv))
        ki = sqrt.(-1im .* ω .* μ₀_1D .* σ[i])
        Zi = (ω .* μ₀_1D) ./ ki
        q = tanh.(1im .* ki .* tv[i])
        Z = Zi .* (Z .+ Zi .* q) ./ (Zi .+ Z .* q)
    end

    _package_mt1d_response(:analytical, fv, Z)
end

"""
    solve_mt1d_fd(f, mesh)

Inputs:
- `f`: Frequencies in hertz.
- `mesh`: 1D finite-difference mesh.

Output:
- `MT1DResponse`: Finite-difference 1D MT response.

Description:
- Solves the 1D MT diffusion equation on the package finite-difference mesh.
"""
function solve_mt1d_fd(
    f::AbstractVector{<:Real},
    mesh::MT1DMesh,
)
    z = mesh.fd_nodes
    n_nodes = length(z)
    n_unknowns = n_nodes - 1
    n_unknowns >= 1 || error("finite-difference mesh needs at least one cell")

    ρnode = Vector{Float64}(undef, n_nodes)
    ρnode[1:end-1] .= mesh.fd_resistivities
    ρnode[end] = mesh.fd_resistivities[end]
    σ = 1.0 ./ ρnode

    Z = Vector{ComplexF64}(undef, length(f))
    for (i, fi) in enumerate(Float64.(f))
        ω = 2π * fi
        k² = -1im * ω * μ₀_1D .* σ
        A = zeros(ComplexF64, n_unknowns, n_unknowns)
        b = zeros(ComplexF64, n_unknowns)

        if n_nodes == 2
            h = z[2] - z[1]
            kbasement = sqrt(-1im * ω * μ₀_1D * σ[end])
            A[1, 1] = 1 / h + 1im * kbasement
            b[1] = -1 / h
        else
            h1 = z[2] - z[1]
            h2 = z[3] - z[2]
            A[1, 1] = -2 / (h1 * h2) + k²[2]
            if n_unknowns > 1
                A[1, 2] = 2 / (h2 * (h1 + h2))
            end
            b[1] = -2 / (h1 * (h1 + h2))

            for node in 3:(n_nodes - 1)
                row = node - 1
                hprev = z[node] - z[node - 1]
                hnext = z[node + 1] - z[node]
                A[row, row - 1] = 2 / (hprev * (hprev + hnext))
                A[row, row] = -2 / (hprev * hnext) + k²[node]
                if row < n_unknowns
                    A[row, row + 1] = 2 / (hnext * (hprev + hnext))
                end
            end

            hlast = z[end] - z[end - 1]
            kbasement = sqrt(-1im * ω * μ₀_1D * σ[end])
            if n_unknowns > 1
                A[end, end - 1] = -1 / hlast
            end
            A[end, end] = 1 / hlast + 1im * kbasement
        end

        E = vcat(1.0 + 0im, A \ b)
        ∂E∂z = (E[2] - E[1]) / (z[2] - z[1])
        Z[i] = -1im * ω * μ₀_1D * E[1] / ∂E∂z
    end

    _package_mt1d_response(:fd, f, Z)
end

"""
    run_mt1d_forward(f, ρ, t; methods=[:analytical, :fd], mesh_kwargs...)

Inputs:
- `f`: Frequencies in hertz.
- `ρ`: Layer resistivities including the basement half-space.
- `t`: Layer thicknesses in metres.
- `methods`: Forward solvers to evaluate.

Output:
- Named tuple with `mesh` and `responses`.

Description:
- Runs one or both 1D forward solvers on a layered model.
"""
function run_mt1d_forward(
    f::AbstractVector{<:Real},
    ρ::AbstractVector{<:Real},
    t::AbstractVector{<:Real};
    methods::AbstractVector{Symbol} = [:analytical, :fd],
    mesh_kwargs...,
)
    mesh = build_mt1d_mesh(t, ρ; mesh_kwargs...)
    responses = Dict{Symbol, MT1DResponse}()

    for method in methods
        if method == :analytical
            responses[method] = solve_mt1d_analytical(f, ρ, t)
        elseif method == :fd
            responses[method] = solve_mt1d_fd(f, mesh)
        else
            error("unknown 1D method: $method")
        end
    end

    (mesh = mesh, responses = responses)
end

"""
    write_mt1d_response_csv(path, responses)

Inputs:
- `path`: Output CSV path.
- `responses`: Response dictionary keyed by method name.

Output:
- `String`: Path to the written CSV file.

Description:
- Writes one or more 1D MT responses to a compact CSV table.
"""
function write_mt1d_response_csv(path::AbstractString, responses::Dict{Symbol, MT1DResponse})
    mkpath(dirname(path))
    ordered = Symbol[]
    for method in (:observed, :true, :recovered, :mean, :median, :analytical, :fd)
        haskey(responses, method) && push!(ordered, method)
    end
    for method in sort!(collect(setdiff(keys(responses), Set(ordered))); by = string)
        push!(ordered, method)
    end

    open(path, "w") do io
        println(io, "frequency_hz,method,apparent_resistivity_ohm_m,phase_deg,re_z,im_z")
        for method in ordered
            response = responses[method]
            for i in eachindex(response.frequencies)
                Printf.@printf(
                    io,
                    "%.8e,%s,%.8e,%.6f,%.8e,%.8e\n",
                    response.frequencies[i],
                    String(method),
                    response.apparent_resistivity[i],
                    response.phase[i],
                    real(response.impedance[i]),
                    imag(response.impedance[i]),
                )
            end
        end
    end

    String(path)
end

"""
    Forward1D(f, ρ, t; kwargs...)

Inputs:
- `f`: Frequencies in hertz.
- `ρ`: Layer resistivities including the basement half-space.
- `t`: Layer thicknesses in metres.

Output:
- Named tuple with `mesh` and `responses`.

Description:
- Public alias for `run_mt1d_forward`.
"""
Forward1D(f::AbstractVector{<:Real}, ρ::AbstractVector{<:Real}, t::AbstractVector{<:Real}; kwargs...) =
    run_mt1d_forward(f, ρ, t; kwargs...)

"""
    RunForward1D(f, ρ, t; kwargs...)

Inputs:
- `f`: Frequencies in hertz.
- `ρ`: Layer resistivities including the basement half-space.
- `t`: Layer thicknesses in metres.

Output:
- Named tuple with `mesh` and `responses`.

Description:
- Compatibility alias for `run_mt1d_forward`.
"""
RunForward1D(f::AbstractVector{<:Real}, ρ::AbstractVector{<:Real}, t::AbstractVector{<:Real}; kwargs...) =
    run_mt1d_forward(f, ρ, t; kwargs...)

"""
    _apply_mt1d_noise(response, specification; rng_seed)

Inputs:
- `response`: Noise-free 1D MT response.
- `specification`: Data errors to apply.
- `rng_seed`: Random-number seed.

Output:
- `MT1DResponse`: Noisy 1D MT response.

Description:
- Perturbs log apparent resistivity and phase using the survey uncertainty specification.
"""
function _apply_mt1d_noise(
    response::MT1DResponse,
    specification::MT1DDataSpec;
    rng_seed::Integer,
)
    rng = Random.MersenneTwister(rng_seed)
    logρ = log10.(response.apparent_resistivity) .+ specification.rho_error_fraction .* randn(rng, length(specification.frequencies))
    ϕ = response.phase .+ specification.phase_error_deg .* randn(rng, length(specification.frequencies))
    Z = similar(response.impedance)

    for i in eachindex(specification.frequencies)
        ω = 2π * specification.frequencies[i]
        magnitude = sqrt((10.0 ^ logρ[i]) * μ₀_1D * ω)
        Z[i] = magnitude * cis(deg2rad(ϕ[i]))
    end

    MT1DResponse(
        method = response.method,
        frequencies = response.frequencies,
        impedance = Z,
        apparent_resistivity = 10 .^ logρ,
        phase = ϕ,
    )
end

"""
    ForwardSolve1D(model_path, data_path; add_noise=false, output_path=nothing, rng_seed=20260308, method=:analytical)

Inputs:
- `model_path`: Saved 1D model file.
- `data_path`: 1D data template or observed-data file providing frequencies and errors.
- `add_noise`, `output_path`, `rng_seed`, `method`: File-driven forward options.

Output:
- `String`: Path to the written predicted or observed data file.

Description:
- Loads a saved 1D model, computes the requested response, optionally adds noise, and writes the result to disk.
"""
function ForwardSolve1D(
    model_path::AbstractString,
    data_path::AbstractString;
    add_noise::Bool = false,
    output_path::Union{Nothing, AbstractString} = nothing,
    rng_seed::Integer = 20260308,
    method::Symbol = :analytical,
)
    mesh = load_mt1d_model(model_path)
    specification = load_mt1d_data_spec(data_path)
    layered = mt1d_layered_model(mesh)
    response = if method == :analytical
        solve_mt1d_analytical(specification.frequencies, layered.resistivities, layered.thicknesses)
    elseif method == :fd
        solve_mt1d_fd(specification.frequencies, mesh)
    else
        error("unsupported 1D forward method: $method")
    end

    observed = add_noise ? _apply_mt1d_noise(response, specification; rng_seed = rng_seed) : response
    destination = something(output_path, joinpath(dirname(data_path), "Data.obs"))
    write_mt1d_observed_data(destination, specification, observed)
end

"""
    _mt1d_style(method)

Inputs:
- Response label symbol.

Output:
- Named tuple with plotting style fields.

Description:
- Returns the plotting style used for a 1D response family.
"""
function _mt1d_style(method::Symbol)
    method == :analytical && return (; color = :navy, linestyle = nothing)
    method == :fd && return (; color = :darkorange, linestyle = :dash)
    (; color = :black, linestyle = nothing)
end

"""
    plot_mt1d_data(responses; output_path)

Inputs:
- Dictionary of 1D responses and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Writes the standard 1D apparent-resistivity and phase plot.
"""
function plot_mt1d_data(
    responses::Dict{Symbol, MT1DResponse};
    output_path::AbstractString,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    figure = Figure(size = (900, 700))
    rho_axis = Axis(
        figure[1, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Apparent resistivity (Ω·m)",
        xscale = log10,
        yscale = log10,
        xreversed = true,
        title = "1D MT response",
    )
    phase_axis = Axis(
        figure[2, 1],
        xlabel = "Frequency (Hz)",
        ylabel = "Phase (deg)",
        xscale = log10,
        xreversed = true,
        title = "Phase",
    )

    methods = collect(keys(responses))
    ordered_methods = Symbol[]
    :observed in methods && push!(ordered_methods, :observed)
    :predicted in methods && push!(ordered_methods, :predicted)
    for method in sort(filter(method -> method ∉ (:observed, :predicted), methods); by = string)
        push!(ordered_methods, method)
    end

    for method in ordered_methods
        response = responses[method]
        if method == :observed
            scatter!(
                rho_axis,
                response.frequencies,
                response.apparent_resistivity;
                color = :white,
                strokecolor = :black,
                strokewidth = 2.5,
                markersize = 13,
                label = "Observed",
            )
            scatter!(
                phase_axis,
                response.frequencies,
                response.phase;
                color = :white,
                strokecolor = :black,
                strokewidth = 2.5,
                markersize = 13,
                label = "Observed",
            )
        elseif method == :predicted
            lines!(rho_axis, response.frequencies, response.apparent_resistivity;
                   color = :firebrick, linewidth = 3, label = "Predicted")
            lines!(phase_axis, response.frequencies, response.phase;
                   color = :firebrick, linewidth = 3, label = "Predicted")
        else
            style = _mt1d_style(method)
            lines!(rho_axis, response.frequencies, response.apparent_resistivity;
                   color = style.color, linestyle = style.linestyle, linewidth = 3, label = String(method))
            lines!(phase_axis, response.frequencies, response.phase;
                   color = style.color, linestyle = style.linestyle, linewidth = 3, label = String(method))
        end
    end

    axislegend(rho_axis, position = :rb)
    axislegend(phase_axis, position = :rb)
    save(output_path, figure)
    String(output_path)
end

"""
    PlotData1D(observed_data_path, predicted_data_path=nothing; output_path)

Inputs:
- 1D observed-data path, optional predicted-data path, and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Loads one or two 1D data files and writes the standard observed-versus-predicted response plot.
"""
function PlotData1D(
    observed_data_path::AbstractString,
    predicted_data_path::Union{Nothing, AbstractString} = nothing;
    output_path::AbstractString,
)
    observed = load_mt1d_observed_data(observed_data_path)
    responses = Dict(:observed => observed.response)
    if predicted_data_path !== nothing
        predicted = load_mt1d_observed_data(predicted_data_path)
        responses[:predicted] = predicted.response
    end
    plot_mt1d_data(responses; output_path = output_path)
end

"""
    plot_mt1d_model(mesh, resistivities; output_path, maximum_depth_m=nothing)

Inputs:
- 1D mesh, layered resistivities, output image path, and optional plotted depth limit.

Output:
- `String`: Path to the written plot.

Description:
- Writes the standard 1D layered-model plot.
"""
function plot_mt1d_model(
    mesh::MT1DMesh,
    resistivities::AbstractVector{<:Real};
    output_path::AbstractString,
    maximum_depth_m::Union{Nothing, Real} = nothing,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    depth_limits = vcat(mesh.interfaces, mesh.fd_nodes[end])
    x_values = Float64[]
    y_values = Float64[]
    for layer_index in eachindex(resistivities)
        push!(x_values, resistivities[layer_index], resistivities[layer_index])
        push!(y_values, depth_limits[layer_index], depth_limits[layer_index + 1])
    end
    max_depth_m = maximum_depth_m === nothing ? depth_limits[end] : min(Float64(maximum_depth_m), depth_limits[end])

    figure = Figure(size = (650, 700))
    axis = Axis(
        figure[1, 1],
        xlabel = "Resistivity (Ω·m)",
        ylabel = "Depth (m)",
        xscale = log10,
        yreversed = true,
        title = "1D layered model (max depth = $(round(max_depth_m, digits = 1)) m)",
    )
    lines!(axis, x_values, y_values, color = :black, linewidth = 3)
    scatter!(axis, resistivities, depth_limits[1:end-1], color = :firebrick, markersize = 10)
    ylims!(axis, max_depth_m, 0.0)

    save(output_path, figure)
    String(output_path)
end

"""
    PlotModel1D(model_path, data_path; output_path, maximum_depth_m=nothing)

Inputs:
- Model path, companion data path, output image path, and optional plotted depth limit.

Output:
- `String`: Path to the written plot.

Description:
- Loads a saved 1D model and writes the standard model plot.
"""
function PlotModel1D(
    model_path::AbstractString,
    data_path::AbstractString;
    output_path::AbstractString,
    maximum_depth_m::Union{Nothing, Real} = nothing,
)
    isfile(data_path) || error("1D data file not found: $data_path")
    mesh = load_mt1d_model(model_path)
    layered = mt1d_layered_model(mesh)
    plot_mt1d_model(mesh, layered.resistivities; output_path = output_path, maximum_depth_m = maximum_depth_m)
end
