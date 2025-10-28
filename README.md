# LevelSetSublimation

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ickaser.github.io/LevelSetSublimation.jl/dev) 
[![](https://zenodo.org/badge/DOI/10.5281/zenodo.17469616.svg)](https://doi.org/10.5281/zenodo.17469616)

## Installation

From the Julia REPL's Pkg mode (open a REPL and type `]` so that the prompt turns blue), add this package with:
```
add https://github.com/Ickaser/LevelSetSublimation.jl
```
Doing so will install the `main` branch of this repo, which will get you the most current version of the code.

# Documentation
If all is well, the documentation can be found [here](https://ickaser.github.io/LevelSetSublimation.jl/). But you can also set up a local website-like interface from your computer with the following steps:
1. navigate to `LevelSetSublimation/docs`, 
2. run `$ julia --project=. make.jl` in your terminal,
3. navigate to `LevelSetSublimation/docs/build`,
4. run `$ python3 -m http.server --bind localhost`
