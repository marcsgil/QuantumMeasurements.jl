using Documenter, QuantumMeasurements

makedocs(
    sitename="QuantumMeasurements.jl",
    pages=[
        "Introduction" => "index.md",
        "Quick-Start" => "quick_start.md",
        "Explanations" => "explanation.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ]
)