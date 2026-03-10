# Contributing to MTGeophysics.jl

Thank you for your interest in contributing! This guide covers how to
report issues, suggest improvements, and submit code.

## Reporting bugs

Open a [GitHub issue](../../issues) with:

- A short description of the problem.
- Steps to reproduce it (input files, commands, Julia version).
- The full error message or unexpected output.

## Suggesting features

Open a GitHub issue labelled **enhancement** describing the use case and
expected behaviour.

## Submitting code

1. Fork the repository and create a branch from `main`.
2. Install the project: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
3. Make your changes.
4. Run the tests: `julia --project=. test/runtests.jl`
5. Open a pull request against `main`.

### Code style

- Follow standard Julia conventions (4-space indent, lowercase functions).
- Add docstrings for new public functions.
- Keep commits focused; one logical change per commit.

### Tests

Add or update tests in `test/` for any new functionality. All tests must
pass before a PR will be merged.

## Code of Conduct

Contributors are expected to be respectful and constructive. Harassment
of any kind will not be tolerated.

## License

By contributing you agree that your contributions will be licensed under
the [MIT License](LICENSE.md).
