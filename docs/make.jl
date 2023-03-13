using DrWatson
@quickactivate :LevelSetSublimation
using Documenter
using DocumenterCitations

bib = CitationBibliography(projectdir("bib", "Zotero_LevelSet.bib"))

makedocs(
        bib,
        sitename="Level Set Sublimation Documentation",
        pages = [
            "Table of Contents" => "index.md",
            "Running a Simulation" => "sim.md",
            "Visualizing Simulation Results" => "plots.md",
            "Implementation and Numerical Details" => [
                "Solved Equations" => "eqs.md",
                "Level Set" => "levelset.md",
                "Heat Equation" => "heat.md",
                "Discretizations" => "disc.md",
            ]
        ] 
        )

