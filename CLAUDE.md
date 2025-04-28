# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build/Test Commands
- Run all tests: `julia --project=. -e "using Pkg; Pkg.test()"`
- Run single test: `julia --project=. -e "using Pkg; Pkg.test(\"RobustGRAPE\", test_args=[\"path/to/test.jl\"])"`

## Code Style Guidelines
- Use 4-space indentation
- Import order: Julia standard libraries first, then external packages, then internal modules
- Type annotations should be used for function arguments and struct fields when helpful
- Use docstrings with Parameters section for public functions and types
- Error handling: use descriptive error messages
- Follow Julia naming conventions: snake_case for functions and variables, CamelCase for types
- Export explicit symbols rather than `export *` patterns
- Parameters.jl `@with_kw` is used selectively, not required for all structs