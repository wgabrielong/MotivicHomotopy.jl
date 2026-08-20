```@meta
CurrentModule = MotivicHomotopy
```

# Computing A¹-Brouwer degrees

The local and global ``\mathbb{A}^1``-Brouwer degrees, stable and
unstable, together with the local-algebra basis helper. Exact over
``\mathbb{Q}`` and finite fields of odd characteristic; over
``\mathbb{C}`` the unstable degrees compute numerically from
`HomotopyContinuation` expressions — or exactly from the coefficients read
as exact Gaussian rationals with `certify = true` (see *Certified
computation* on the manual's front page); over
``\mathbb{R}`` compute over ``\mathbb{Q}`` and base-change the result.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["degrees.jl", "unstable_degrees.jl"]
```
