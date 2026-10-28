using MotivicHomotopy
using Documenter

DocMeta.setdocmeta!(MotivicHomotopy, :DocTestSetup, :(using MotivicHomotopy); recursive=true)

makedocs(;
    modules=[MotivicHomotopy],
    authors="The MotivicHomotopy.jl contributors",
    sitename="MotivicHomotopy.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/stable",
        edit_link="dev",
        sidebar_sitename=false,
        assets=["assets/favicon.ico"],
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

# Each branch's documentation is built from that branch's code: `main`, which
# holds the released code, is published at /stable, and `dev`, the development
# branch, at /dev. Other builds (pull requests, tags, local builds) are not
# deployed.
const DOCS_FOLDERS = Dict("main" => "stable", "dev" => "dev")
const BUILD_BRANCH = replace(get(ENV, "GITHUB_REF", ""), r"^refs/heads/" => "")

if haskey(DOCS_FOLDERS, BUILD_BRANCH)
    deploydocs(;
        repo="github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl",
        devbranch=BUILD_BRANCH,
        devurl=DOCS_FOLDERS[BUILD_BRANCH],
        versions=["stable" => "stable", "dev" => "dev"],
    )
end
