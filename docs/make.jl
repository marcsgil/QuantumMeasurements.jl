using Documenter, DocumenterVitepress, Literate

using QuantumMeasurements

dir = pkgdir(QuantumMeasurements)

for file ∈ readdir(joinpath(dir, "examples"), join=true)
    Literate.markdown(file, joinpath(dir, "docs/src"); documenter=true)
end

makedocs(;
    modules=[QuantumMeasurements],
    authors="Marcos Gil",
    repo="https://github.com/marcsgil/QuantumMeasurements.jl",
    sitename="QuantumMeasurements.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/marcsgil/QuantumMeasurements.jl",
        devurl="dev",
        deploy_url="marcsgil.github.io/QuantumMeasurements.jl",
    ),
    pages=[
        "QuantumMeasurements.jl" => "index.md",
        "Quick-Start" => "quick_start.md",
        "Explanations" => "explanation.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/marcsgil/QuantumMeasurements.jl",
    push_preview=true,
)
