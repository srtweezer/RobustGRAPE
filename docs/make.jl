using Documenter

# Add the package directory to the load path
push!(LOAD_PATH, "../src/")

# Load the package - this is needed when the package isn't registered
try
    @eval using RobustGRAPE
catch e
    @warn "Error loading RobustGRAPE" exception=e
    # In case there's an issue, try to load it from the module directly
    include("../src/RobustGRAPE.jl")
    @eval using Main.RobustGRAPE
end

makeinfo = (
    modules = [RobustGRAPE],
    format = Documenter.HTML(
        prettyurls = !isempty(get(ENV, "CI", "")),
        canonical = "https://srtweezer.github.io/RobustGRAPE.jl/stable/",
        assets = [],
        sidebar_sitename = false,
    ),
    sitename = "RobustGRAPE.jl",
    authors = "Endres Lab",
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Introduction" => "guide/introduction.md",
            "Getting Started" => "guide/getting-started.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Unitary Calculations" => "api/unitary.md",
            "Fidelity Calculations" => "api/fidelity.md",
            "Regularization" => "api/regularization.md",
            "Rydberg Tools" => "api/rydberg.md",
        ],
        "Examples" => "examples.md",
    ],
    repo = "https://github.com/srtweezer/RobustGRAPE.jl/blob/{commit}{path}#L{line}",
    warnonly = true
)

Documenter.makedocs(; makeinfo...)

# Uncomment this when ready to deploy to GitHub Pages
Documenter.deploydocs(
    repo = "github.com/srtweezer/RobustGRAPE.jl.git",
    devbranch = "main",
    push_preview = true,
)