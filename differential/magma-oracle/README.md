# Magma oracle run — request for help

**Possibly no licence needed:** Magma's free online calculator
(http://magma.maths.usyd.edu.au/calc/) has historically limited submissions
to about two minutes of computation and cannot load user packages. We could
not verify the current limits ourselves. If those limits still hold, Tier 1
(see below) for most of the records here should fit: paste the whole of
`oracle.m` into the calculator, and if it times out, delete some of the
later `Append(~records, ...)` lines — every record is independent, and the
one expensive record (dimension 10 over a degree-10 field) is deliberately
last. Even a handful of records run this way is valuable to us.

## What we are asking

Run the single file `oracle.m` in Magma and send us the text it prints.
That is the entire request.

## Why

We have written an open-source Julia implementation of the algorithms of
P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic form
over a number field* (arXiv:2109.04172), and we want to cross-validate it
against the reference implementation — Koprowski's CQF package for Magma —
and against Magma's own built-in quadratic-form arithmetic. **We have no
Magma licence and cannot run this script ourselves.** That is why it is
deliberately minimal: it only computes built-in quantities and prints them,
with the name of the Magma function recorded next to every value, so that
you can read it in a few minutes and satisfy yourself about what it does
before running it. All checking happens on our side afterwards.

## The two tiers

- **Tier 1** needs plain Magma only: determinants, signs at the real
  embeddings, Hilbert symbols at the relevant primes, isotropy, Witt
  indices. It runs first and is complete on its own.

  **If you are short of time, Tier 1 is the more valuable half**, not the
  lesser one. Magma's built-in quadratic-form arithmetic is an independent
  implementation of the underlying invariants, while CQF is the reference
  implementation of the *same algorithm* we implemented — so Tier 1 output
  is evidence from an independent code path, and Tier 2 output, valuable as
  it is, could in principle share a mistake with us by construction. Someone
  who runs only Tier 1 has given us the stronger result and should know
  that, rather than assume they have given us the weaker one.
- **Tier 2** computes the anisotropic parts themselves and needs the CQF
  package: P. Koprowski, *CQF Magma package*, ACM Communications in
  Computer Algebra 54(2) (2020) 53–56 — available from the author
  (Przemysław Koprowski, University of Silesia; see the paper). If CQF is
  not attached, the script says so in one line, skips Tier 2 cleanly, and
  the run is still a success.

To enable Tier 2, uncomment and adjust the single `AttachSpec` line at the
top of `oracle.m`. Note that we could not test the exact CQF calling
convention; if the `AnisotropicPart` call errors, the error is printed and
the run continues — a one-line correction from someone with the CQF
documentation in front of them would be very welcome.

## How to run it

```
magma -b oracle.m > results.txt
```

Nothing to configure: no paths, no environment variables, no arguments.
Records are printed one at a time, cheapest first, and flushed as they
complete — **partial output is genuinely welcome. If the run is taking
longer than you want to give it, kill it at any point and send whatever it
printed.** The file it produces is plain text.

## Expected runtime

An estimate, extrapolated from the timings in §7 of the arXiv paper (their
dimension-10 forms over a degree-10 field averaged 2.82 s each in Magma):
almost all records are rank ≤ 8 over quadratic fields and individually
trivial, so Tier 1 over all records should take well under ten minutes on a
modern machine, and Tier 1 plus Tier 2 under twenty; the last record (the
degree-10 field) dominates. If it is far slower than that on your machine,
killing it early is fine — see above.

## Where to send the output

Either of:

- a pull request to this repository adding your `results.txt` under
  `differential/magma-oracle/results/`, or
- a plain email with `results.txt` attached to the maintainer:
  gabriel_ong@ieee.org.

Please include the Magma version if you trimmed the header, and — if you
are willing — tell us how you would like to be credited (name and
affiliation, a link, or not at all); we will acknowledge the contribution
in the package documentation either way unless you prefer otherwise.

**If some value disagrees with ours, that is expected and interesting — it
is the entire point of the exercise, and it does not mean you did something
wrong.** The corpus was built to force every branch of the implementations
through inputs chosen to be hard; we have already found one gap in the
published algorithm this way, so divergences are treated as findings to be
investigated on our side, never as a problem with your run.

## Files here

- `oracle.m` — the file to run. Generated (fields and forms come from this
  repository's test fixtures via `generate_oracle.jl`); do not edit records
  by hand.
- `generate_oracle.jl` — the generator. Julia-side only; you do not need it.
- `corpus.jl` — the corpus definition: every record annotated with the
  internal branch of our implementation it was chosen to exercise. See
  `../CORPUS.md` for the human-readable version and the measured coverage
  matrix.
- `coverage.jl` — the Julia-side coverage measurement the corpus was built
  with. You do not need it.
- `results-stub.txt` — a **hand-computed** example of what output looks
  like, for format reference. It is not Magma output.
