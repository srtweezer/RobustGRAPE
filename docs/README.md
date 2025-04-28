# RobustGRAPE.jl Documentation

This directory contains the documentation for the RobustGRAPE.jl package.

## Building the Documentation

To build the documentation, you can use the provided script in the root directory:

```bash
./build_docs.jl
```

Alternatively, you can build the documentation manually with the following steps:

1. Activate the documentation environment:
   ```julia
   julia --project=docs
   ```

2. Install dependencies:
   ```julia
   import Pkg; Pkg.develop(path="."); Pkg.instantiate()
   ```

3. Build the documentation:
   ```julia
   include("docs/make.jl")
   ```

## Viewing the Documentation

After building, you can view the documentation by opening `docs/build/index.html` in your web browser.

## Documentation Structure

- `make.jl`: Documentation build script
- `Project.toml`: Dependencies for the documentation
- `src/`: Source files for the documentation
  - `index.md`: Home page
  - `guide/`: User guide pages
  - `api/`: API reference pages
  - `examples.md`: Example usage

## Editing the Documentation

The documentation is generated from docstrings in the source code and Markdown files in the `src/` directory. To add or modify documentation:

1. **API Reference**: Edit the docstrings in the source code files
2. **User Guide**: Edit the Markdown files in the `src/guide/` directory
3. **Examples**: Edit the `src/examples.md` file or add new example files and update `make.jl`

## Publishing the Documentation

To publish the documentation to GitHub Pages (when the repository is on GitHub), the `deploydocs` function in `make.jl` is already configured. You'll need to:

1. Update the repository URL in `make.jl` to your actual repository URL
2. Set up the GitHub repository with GitHub Actions to build and deploy the documentation

See the [Documenter.jl documentation](https://juliadocs.github.io/Documenter.jl/stable/man/hosting/) for more details.