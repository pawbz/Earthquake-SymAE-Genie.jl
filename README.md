## Genie app to viz. earthquake source time functions derived from SymAE

## Installation


Clone the repository and install the dependencies:

First `cd` into the project directory then run:

```bash
$> julia --project -e 'using Pkg; Pkg.instantiate()'
```


Finally, run the app

```bash
$> julia --project
```

```julia
julia> using GenieFramework
julia> Genie.loadapp() # load app
julia> up() # start server
```
Finally, open your browser and navigate to `http://localhost:8000/` to use the app.
