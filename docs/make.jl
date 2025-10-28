using LevelSetSublimation
using Documenter
using DocumenterCitations

bib = CitationBibliography(projectdir("bib", "Zotero_LevelSet.bib");
        style=:authoryear)

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

