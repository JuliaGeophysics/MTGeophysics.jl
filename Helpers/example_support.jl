# This helper script provides a shared result-directory helper for repository example scripts.

using Dates

"""
    create_example_result_dir(script_path; root_dir=joinpath(dirname(dirname(script_path)), "Results"), now_value=now())

Inputs:
- Example script path, optional results root, and optional timestamp.

Output:
- Named tuple with `run_dir`, `timestamp_slug`, and `timestamp_label`.

Description:
- Builds the timestamped results directory used by the example scripts.
"""
function create_example_result_dir(
    script_path::AbstractString;
    root_dir::AbstractString = joinpath(dirname(dirname(script_path)), "Results"),
    now_value::DateTime = now(),
)
    script_name = splitext(basename(script_path))[1]
    timestamp_slug = Dates.format(now_value, "yyyymmdd_HHMMSS")
    timestamp_label = Dates.format(now_value, "yyyy-mm-dd HH:MM:SS")
    run_dir = joinpath(root_dir, "$(script_name)_$(timestamp_slug)")
    mkpath(run_dir)
    (run_dir = run_dir, timestamp_slug = timestamp_slug, timestamp_label = timestamp_label)
end
