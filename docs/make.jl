using Documenter, DocumenterVitepress, Literate

using QuantumMeasurements

dir = pkgdir(QuantumMeasurements)
repo = "https://github.com/marcsgil/QuantumMeasurements.jl"

for file ∈ readdir(joinpath(dir, "examples"), join=true)
    if endswith(file, ".jl")
        Literate.markdown(file, joinpath(dir, "docs/src"); documenter=true, repo_root_url=joinpath(repo, "tree/master"))
    end
end

makedocs(;
    modules=[QuantumMeasurements],
    authors="Marcos Gil",
    repo,
    sitename="QuantumMeasurements.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/marcsgil/QuantumMeasurements.jl",
    ),
    pages=[
        "QuantumMeasurements.jl" => "index.md",
        "Quick-Start" => "quick_start.md",
        "Explanations" => "explanation.md",
        "Examples" => ["twin_photons_pvm.md"],
        "API" => "api.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/marcsgil/QuantumMeasurements.jl",
    push_preview=true,
)
