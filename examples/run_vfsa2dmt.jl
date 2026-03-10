# This script runs a 2D VFSA inversion and writes all outputs into a
# timestamped result directory beside itself.

using MTGeophysics

result = VFSA2DMT(
    VFSA2DMTParams(
        script_path = @__FILE__,
        start_model_path = normpath(@__DIR__, "0COMEMI2D-I", "Comemi2D1.ini"),
        data_path = normpath(@__DIR__, "0COMEMI2D-I", "Comemi2D1.obs"),
        config = VFSA2DMTConfig(
            n_chains = 2,
            n_ctrl = 400,
            max_iter = 3000, # 3000 for a real run; 300 for smoke test
            n_trials = 4,
            log_bounds = (0.0, 4.0),
            step_scale = 0.11,
            seed = 20260308,
            keep_models = true,
        ),
    ),
)

println("ResultDir = ", result.run_info.run_dir)
println("Summary = ", result.summary_path)
