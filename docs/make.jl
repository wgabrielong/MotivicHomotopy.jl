using MotivicHomotopy
using Documenter

DocMeta.setdocmeta!(MotivicHomotopy, :DocTestSetup, :(using MotivicHomotopy); recursive=true)

makedocs(;
    modules=[MotivicHomotopy],
    authors="The MotivicHomotopy.jl contributors",
    sitename="MotivicHomotopy.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl",
        edit_link="main",
        sidebar_sitename=false,
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Working with forms" => "forms.md",
        "Arithmetic methods" => "arithmetic.md",
        "Computing A¹-Brouwer degrees" => "degrees.md",
    ],
    # The examples in docstrings are captured REPL sessions, not run at build
    # time (they require Oscar to be loaded and exact display matching).
    doctest=false,
    # Internal (underscore) functions carry docstrings for readers of the
    # source — the KR algorithm functions are documented line-by-line against
    # the paper — but they are not part of the manual; only exported names
    # must appear in it.
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl",
    devbranch="main",
)
