CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using LevelSetSublimation
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath((@__DIR__), "bib", "Zotero_LevelSet.bib");
        style=:authoryear)

@info "Building Documentation"
makedocs(
        sitename="Level Set Sublimation Documentation",
        pages = [
            "Table of Contents" => "index.md",
            "Running a Simulation" => "sim.md",
            "Visualizing Simulation Results" => "plots.md",
            "Implementation and Numerical Details" => [
                "Solved Equations" => "eqs.md",
                "Solving Physical Equations" => "heatmass.md",
                "Level Set" => "levelset.md",
                "Discretizations" => "disc.md",
            ],
            "References" => "refs.md",
            "Everything for Good Measure" => "alldocs.md",
        ]
        ; plugins = [bib], 
        )

if CI
    @info "Deploying Documentation"
    deploydocs(
        # `repo` MUST be set correctly. Once your GitHub name is set
        # the auto-generated documentation will be hosted at:
        # https://PutYourGitHubNameHere.github.io/LyoPronto.jl/dev/
        # (assuming you have enabled `gh-pages` deployment)
        repo = "github.com/ickaser/LevelSetSublimation.jl.git",
        target = "build",
        push_preview = true,
        devbranch = "main",
    )
end

@info "Finished with Documentation"

