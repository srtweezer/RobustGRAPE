#!/usr/bin/env julia

# This script builds the documentation for RobustGRAPE.jl without requiring registration

# Get the package directory
pkg_path = dirname(@__FILE__)

# Add the source directory to the load path
push!(LOAD_PATH, pkg_path)

# Clean up the build directory if it exists
build_dir = joinpath(pkg_path, "docs", "build")
if isdir(build_dir)
    println("Cleaning existing documentation build directory...")
    rm(build_dir, recursive=true, force=true)
    println("Build directory cleared.")
end

# Ensure the docs project environment is activated and packages are installed
import Pkg
Pkg.activate(pkg_path)
Pkg.resolve()
Pkg.instantiate()

# Develop the main package from the current directory
#Pkg.develop(path=pkg_path)

# Try to load the package
try
    @eval using RobustGRAPE
catch e
    @warn "Failed to load RobustGRAPE normally, trying direct include" exception=e
    # If loading fails, try to include the module directly
    include(joinpath(pkg_path, "src", "RobustGRAPE.jl"))
    @eval using Main.RobustGRAPE
end

# Build the documentation
include(joinpath(pkg_path, "docs", "make.jl"))

println("\nDocumentation build attempt completed.")
println("If successful, you can view the documentation by opening:")
println(joinpath(pkg_path, "docs", "build", "index.html"))