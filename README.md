# LevelSetSublimation

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> LevelSetSublimation

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate :LevelSetSublimation
```
which auto-activates the project and enables local path handling from DrWatson.

# Documentation
The documentation is set up to run with Documenter.jl , which is where you can find more extensive instructions on how to navigate it...

but a nice way to set up a website-like interface (running from your computer) is to follow these steps: 
1. navigate to `LevelSetSublimation/docs`, 
2. run `$ julia make.jl` in your terminal,
3. navigate to `LevelSetSublimation/docs/build`,
4. run `$ python3 -m http.server --bind localhost`
