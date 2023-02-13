using DrWatson
@quickactivate :LevelSetSublimation
using Documenter

makedocs(
        sitename="Level Set Sublimation Documentation",
        pages = [
            "Table of Contents" => "index.md",
            "Running a Simulation" => "sim.md",
            "Implementation and Details" => [
                "Level Set" => "levelset.md",
                "Heat Equation" => "heat.md",
                "Discretizations" => "disc.md"
            ]
        ] 
        )

