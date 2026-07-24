# Generates oracle.m and ../CORPUS.md (committed alongside this script) — the
# self-contained Magma file described in README.md, plus the corpus-intent
# record.
#
# Run from the repository root, with the number-fields sources present:
#     julia --project=. differential/magma-oracle/generate_oracle.jl
#
# The corpus is defined in corpus.jl; the base test fields are read from
# test/nf_test_fields.jl — the same objects the test suite uses — and the
# per-form relevant-prime lists come from the package's own
# S = 𝔓(q) ∪ {dyadic} data, so the Magma corpus and the Julia corpus cannot
# drift apart. Before emitting anything, the whole corpus is run through the
# instrumented number-field layer (coverage.jl / src/nf_coverage.jl) and the
# MEASURED coverage matrix — not the declared intent — is embedded as a
# comment block in oracle.m and as the table in ../CORPUS.md. Only
# field/form/prime LITERALS are generated; the Magma logic below is a fixed
# template, kept minimal so a Magma user can audit it before running (we
# cannot run Magma ourselves).

include(joinpath(@__DIR__, "corpus.jl"))
include(joinpath(@__DIR__, "coverage.jl"))

println("verifying the field properties the coverage argument rests on...")
_verify_corpus_fields()

println("measuring branch coverage of the corpus (instrumented run)...")
measurements, cov_matrix = measure_corpus(corpus)

# every row must have >= 3 covering records unless CORPUS.md argues the
# exception; a regeneration that silently loses coverage should be loud
for (name, c, _) in cov_matrix
    c < 3 && println("WARNING: branch \"$name\" has only $c covering record(s) — ",
                     "either fix the corpus or argue the exception in corpus.jl")
end

# rational primes at which the Hasse invariant can be nontrivial, from the
# package's fixed-diagonalisation S
function record_primes(K, entries)
    if K === :Q
        beta = diagonal_form(QQ, Tuple(entries))
        return sort(unique(vcat(ZZ.(relevant_primes(beta)), [ZZ(2)])))
    end
    beta = diagonal_form(K, Tuple(entries))
    sort(unique([minimum(P) for P in MH._kr_form(beta).S]))
end

rat(c::QQFieldElem) = isone(denominator(c)) ? string(numerator(c)) :
    string(numerator(c)) * "/" * string(denominator(c))

function magma_elt(a, t::String)
    a isa Integer && return "($a)"
    terms = String[]
    for i in 0:(degree(parent(a)) - 1)
        c = coeff(a, i)
        iszero(c) && continue
        push!(terms, i == 0 ? "($(rat(c)))" : "($(rat(c)))*$t^$i")
    end
    isempty(terms) ? "0" : join(terms, " + ")
end

function magma_poly(K)
    f = defining_polynomial(K)
    join(["($(rat(coeff(f, i))))*x^$i" for i in 0:degree(f) if !iszero(coeff(f, i))],
         " + ")
end

# assign Magma variable names per distinct field, in corpus order
fieldvars = Dict{Any, Tuple{String, String}}()   # K => (Kname, tname)
fieldorder = Any[]
for (_, K, _, _) in corpus
    haskey(fieldvars, K) && continue
    idx = length(fieldorder)
    fieldvars[K] = ("K$idx", "t$idx")
    push!(fieldorder, K)
end

# the measured coverage matrix as plain text lines (shared by oracle.m and
# CORPUS.md)
matrix_lines = String[]
let w = maximum(length(name) for (name, _, _) in cov_matrix)
    for (name, c, idxs) in cov_matrix
        push!(matrix_lines,
              rpad(name, w) * "  " * lpad(string(c), 3) * "  records " *
              (isempty(idxs) ? "-" : join(idxs, ",")))
    end
end

io = IOBuffer()

print(io, """
// ============================================================================
// oracle.m -- differential-test oracle for the MotivicHomotopy.jl number-field
// layer, to be run in Magma. GENERATED FILE: produced by generate_oracle.jl
// (this repository, branch differential-tests); do not edit records by hand.
//
// Run:      magma -b oracle.m > results.txt
// Tier 1 needs plain Magma only. Tier 2 (the anisotropic part itself) needs
// Przemyslaw Koprowski's CQF package (ACM Commun. Comput. Algebra 54(2)
// (2020) 53-56). If you have CQF, attach it before the run by uncommenting
// and adjusting the single line below; otherwise Tier 2 is skipped cleanly.
//
// AttachSpec("/path/to/cqf.spec");   // <-- Tier 2 only; see README.md
//
// Everything printed is a NATIVE Magma value, labelled with the Magma
// function that produced it; the one composite quantity is the explicit
// product prod_{i<j} HilbertSymbol(d_i, d_j, P), printed with its formula in
// the label. Places are identified by DATA, never by position: real places
// by the numerical value of the field generator under the embedding, finite
// places by (p, e, f, generators). Records are ordered cheapest-first and
// flushed as they complete -- a partial run is still useful.
//
// The corpus is COVERAGE-DRIVEN: every record was chosen to force a specific
// internal branch of our implementation (see ../CORPUS.md for the record-by-
// record intent). The matrix below is MEASURED, by running each record
// through the instrumented implementation at generation time -- branch
// against the record ids that exercise it:
//
""")
for line in matrix_lines
    println(io, "// ", line)
end
print(io, """
// ============================================================================

SetColumns(0);
v1, v2, v3 := GetVersion();
printf "MAGMA-ORACLE-FORMAT 1\\n";
printf "magma_version: %o.%o-%o\\n", v1, v2, v3;
printf "run_date: ";
System("date");
CQFAvailable := true;
try
  _ := eval "AnisotropicPart";
catch e
  CQFAvailable := false;
end try;
printf "cqf_available: %o\\n", CQFAvailable;
if not CQFAvailable then
  printf "note: Tier 2 (the anisotropic part) will be skipped -- CQF is not attached. Tier 1 output is complete and useful on its own.\\n";
end if;

// Tier 1: plain Magma. Prints, per form: determinant; per real place the
// numerical embedding of the generator and the sign vector of the diagonal;
// per finite place over the listed rational primes its (p, e, f, gens) and
// the Hilbert-symbol product; global isotropy and Witt index.
procedure Tier1(id, label, K, diag, primes)
  printf "=== RECORD %o (tier 1) ===\\n", id;
  printf "label: %o\\n", label;
  printf "field_poly [Magma: DefiningPolynomial]: %o\\n", DefiningPolynomial(K);
  printf "diag: %o\\n", diag;
  n := #diag;
  try
    M := DiagonalMatrix(K, diag);
    printf "determinant [Magma: Determinant]: %o\\n", Determinant(M);
    for v in RealPlaces(K) do
      printf "real_place theta_approx [Magma: Evaluate(K.1, place)]: %o\\n", Evaluate(K.1, v);
      printf "real_place sign_vector [Magma: Sign(Evaluate(d_i, place))]: %o\\n", [Sign(Evaluate(d, v)) : d in diag];
    end for;
    O := MaximalOrder(K);
    for p in primes do
      for pair in Decomposition(O, p) do
        P := pair[1];
        h := 1;
        for i in [1..n] do
          for j in [i+1..n] do
            h := h * HilbertSymbol(K!diag[i], K!diag[j], P);
          end for;
        end for;
        printf "prime p=%o e=%o f=%o gens=%o hasse=prod_(i<j)HilbertSymbol(d_i,d_j,P) [Magma: HilbertSymbol]: %o\\n",
          p, RamificationIndex(P), InertiaDegree(P), Generators(P), h;
      end for;
    end for;
    V := QuadraticSpace(M);
    printf "is_isotropic [Magma: IsIsotropic]: %o\\n", IsIsotropic(V);
    printf "witt_index [Magma: WittIndex]: %o\\n", WittIndex(V);
  catch e
    printf "RECORD-ERROR (tier 1, record %o): %o\\n", id, e`Object;
  end try;
  printf "=== END RECORD %o (tier 1) ===\\n", id;
end procedure;

// Tier 2: CQF only. The anisotropic part, printed raw. CQF's calling
// convention could not be tested on our side (no Magma licence); if the
// call errors, the message is printed and the run continues -- and we would
// be grateful for a one-line fix matching the CQF documentation.
procedure Tier2(id, label, K, diag, primes)
  printf "=== RECORD %o (tier 2) ===\\n", id;
  try
    M := DiagonalMatrix(K, diag);
    ap := (eval "AnisotropicPart")(M);
    printf "anisotropic_part [CQF: AnisotropicPart(gram matrix)]: %o\\n", ap;
  catch e
    try
      ap2 := (eval "AnisotropicPart")(diag);
      printf "anisotropic_part [CQF: AnisotropicPart(diagonal sequence)]: %o\\n", ap2;
    catch e2
      printf "RECORD-ERROR (tier 2, record %o): %o\\n", id, e2`Object;
    end try;
  end try;
  printf "=== END RECORD %o (tier 2) ===\\n", id;
end procedure;

QQx<x> := PolynomialRing(Rationals());
""")

# field declarations
for K in fieldorder
    Kname, tname = fieldvars[K]
    if K === :Q
        println(io, "$Kname := RationalsAsNumberField();")
    else
        println(io, "$Kname<$tname> := NumberField($(magma_poly(K)));")
    end
end

println(io, "\nrecords := [* *];")
for (i, (label, K, entries, covers)) in enumerate(corpus)
    Kname, tname = fieldvars[K]
    diag = join([magma_elt(a, tname) for a in entries], ", ")
    primes = join(string.(record_primes(K, entries)), ", ")
    println(io, "// record $i chosen for: ", join(covers, "; "))
    println(io, "Append(~records, [* $i, \"$label\", $Kname, [$Kname | $diag], [$primes] *]);")
end

print(io, """

// Tier 1 first, in full; Tier 2 afterwards, only with CQF. Each record is
// flushed as it completes.
for r in records do
  Tier1(r[1], r[2], r[3], r[4], r[5]);
end for;
if CQFAvailable then
  for r in records do
    Tier2(r[1], r[2], r[3], r[4], r[5]);
  end for;
end if;
printf "MAGMA-ORACLE-DONE\\n";
quit;
""")

outpath = joinpath(@__DIR__, "oracle.m")
write(outpath, String(take!(io)))
println("wrote ", outpath)

# ---------------------------------------------------------------------------
# ../CORPUS.md — the record-by-record intent and the measured coverage matrix
# ---------------------------------------------------------------------------

md = IOBuffer()
print(md, """
# The differential corpus: what each record is for

GENERATED FILE — produced by `magma-oracle/generate_oracle.jl` together with
`oracle.m`; edit `magma-oracle/corpus.jl` and regenerate.

The corpus is coverage-driven. Each record was chosen to force a specific
internal branch of the number-field anisotropic-part implementation
([KR21] Algorithms 1–8 and the [KC18] anisotropic dimension), so that the
external oracle (Magma / CQF) witnesses every branch from inputs and outputs
alone. The instrumentation that measured this is `src/nf_coverage.jl`
(off by default) driven by `magma-oracle/coverage.jl`. When a divergence is
reported for record N, the first question — "what was this form for" — is
answered by the table below.

## Fields

| Magma var | Field | Chosen because |
|---|---|---|
""")
fieldwhy = Dict{Any, String}(
    :Q => "calibration only — the Macaulay2-validated rational path",
)
for (K, why) in field_reasons          # defined in corpus.jl
    fieldwhy[K] = why
end
for K in fieldorder
    Kname, _ = fieldvars[K]
    desc = K === :Q ? "Q (as a degree-1 number field)" :
        "Q[x]/(" * string(defining_polynomial(K)) * ")"
    println(md, "| `$Kname` | $desc | $(get(fieldwhy, K, "")) |")
end

print(md, """

## Records

| # | Field | Diagonal | Chosen to exercise |
|---|---|---|---|
""")
for (i, (label, K, entries, covers)) in enumerate(corpus)
    Kname, _ = fieldvars[K]
    diagstr = "⟨" * join(string.(entries), ", ") * "⟩"
    println(md, "| $i | `$Kname` | $diagstr | ", join(covers, "; "), " |")
end

print(md, """

## The measured coverage matrix

Measured at generation time by running every record through the
instrumented implementation (branch → number of covering records → which):

```
""")
for line in matrix_lines
    println(md, line)
end
print(md, """
```

## Coverage limits — what stays below three records, and why

""")
print(md, coverage_limits_md)             # defined in corpus.jl
write(joinpath(@__DIR__, "..", "CORPUS.md"), String(take!(md)))
println("wrote ", joinpath(@__DIR__, "..", "CORPUS.md"))
