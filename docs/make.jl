using Documenter, Literate, DocumenterCitations

using QuantumMeasurements

dir = pkgdir(QuantumMeasurements)
repo = "https://github.com/marcsgil/QuantumMeasurements.jl"

for file ∈ readdir(joinpath(dir, "examples"), join=true)
    if endswith(file, ".jl")
        Literate.markdown(file, joinpath(dir, "docs/src"); documenter=true, repo_root_url=joinpath(repo, "tree/master"))
    end
end

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[QuantumMeasurements],
    authors="Marcos Gil",
    sitename="QuantumMeasurements.jl",
    pages=[
        "Home" => "index.md",
        "Quick-Start" => "quick_start.md",
        "Examples" => ["twin_photons_pvm.md", "twin_photons_prop.md", "spatial_structure.md",
            "obstructed_spatial_structure.md"],
        "Theory" => ["random_states.md", "mathematical_foundations.md", "maximum_likelihood.md", "bayesian_inference.md", "proportional_measurements.md"],
        "How-to Guides" => ["choosing_methods.md"],
        "API" => "api.md",
        "References" => "references.md",
    ],
    warnonly=true,
    plugins=[bib],
)

deploydocs(;
    repo="github.com/marcsgil/QuantumMeasurements.jl",
)
