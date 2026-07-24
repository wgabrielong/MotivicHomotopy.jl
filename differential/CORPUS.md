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
| `K0` | Q (as a degree-1 number field) | calibration only — the Macaulay2-validated rational path |
| `K1` | Q[x]/(x^2 - 2) | totally real, 2 ramified, h = 1: the test suite's Algorithm 4 correction field |
| `K2` | Q[x]/(x^2 + 21) | imaginary, Cl = (Z/2)^2, 2 ramified non-principal: lambda_b fires with unit coefficients, is suppressed once S grows |
| `K3` | Q[x]/(x^2 - 34) | totally real, h = 2, dyadic prime PRINCIPAL: the only field here where lambda_b coexists with real places, and the only scanned field with a blind-spot witness class |
| `K4` | Q[x]/(x^2 + 2) | imaginary, 2 ramified, h = 1: non-real control isolating dyadic ramification from formal reality |
| `K5` | Q[x]/(x^2 - 10) | totally real, 2 ramified, h = 2: Algorithm 4 dyadic correction + mixed signs; lambda_b provably never fires here (see the coverage-limits section below) |
| `K6` | Q[x]/(x^2 + 7) | imaginary, 2 split, h = 1: KR's own example field; completions at the dyadic places are Q_2 itself |
| `K7` | Q[x]/(x^2 - 5) | totally real, 2 inert, h = 1 |
| `K8` | Q[x]/(x^2 + 5) | imaginary, h = 2, dyadic class generates Cl: the recorded Algorithm 1 divergence field |
| `K9` | Q[x]/(x^3 - 3*x - 1) | totally real cubic, three real places: genuinely non-constant sign patterns |
| `K10` | Q[x]/(x^3 - 2) | mixed signature (one real, one complex place) |
| `K11` | Q[x]/(x^4 - 4*x^2 + 1) | totally real quartic |
| `K12` | Q[x]/(x^4 - 4*x^2 + 2) | Q(zeta_16)^+: 2 totally ramified, e = 4 — wild dyadic ramification, the deepest prime-2 corner |
| `K13` | Q[x]/(zx^10 - 10*zx^8 + 35*zx^6 + zx^5 - 50*zx^4 - 5*zx^3 + 25*zx^2 + 5*zx - 1) | degree 10, totally real: the FacElem scale stress, one expensive record only |

## Records

| # | Field | Diagonal | Chosen to exercise |
|---|---|---|---|
| 1 | `K0` | ⟨1, 1, -3, -3, 1, -1⟩ | calibration |
| 2 | `K1` | ⟨1, sqrt(2)⟩ | dim_a = 2 over the ramified-dyadic real field; A block nonempty |
| 3 | `K1` | ⟨1, -1⟩ | dim_a = 0 |
| 4 | `K2` | ⟨sqrt(-21), -sqrt(-21)⟩ | dim_a = 0 |
| 5 | `K2` | ⟨1, sqrt(-21)⟩ | dim_a = 2, non-real; A block empty; lambda_b suppressed by S-classes (negative control) |
| 6 | `K2` | ⟨1, 1⟩ | lambda_b in E_S basis; shape: all coefficients units |
| 7 | `K3` | ⟨1, -6*sqrt(34) - 35⟩ | dimension correction decisive; lambda_b in E_S basis; shape: all coefficients units |
| 8 | `K3` | ⟨2, -1⟩ | dimension correction decisive; lambda_b in E_S basis |
| 9 | `K3` | ⟨6*sqrt(34) + 35, -1⟩ | dimension correction decisive; lambda_b in E_S basis |
| 10 | `K3` | ⟨1, -sqrt(34)⟩ | lambda_b AND a nonempty A block |
| 11 | `K3` | ⟨1, 1⟩ | lambda_b AND a nonempty A block; shape: all coefficients units |
| 12 | `K3` | ⟨-1, -6*sqrt(34) - 35⟩ | lambda_b AND a nonempty A block; correction non-decisive contrast |
| 13 | `K4` | ⟨1, sqrt(-2)⟩ | dim_a = 2, ramified dyadic without real places; A block empty |
| 14 | `K5` | ⟨1, -sqrt(10)⟩ | dim_a = 2 over the triple-point field; A block nonempty |
| 15 | `K6` | ⟨1, 2, -2⟩ | dim_a = 1 |
| 16 | `K4` | ⟨3, 1, -3⟩ | dim_a = 1 |
| 17 | `K3` | ⟨sqrt(34), 2, -1⟩ | dim_a = 1; shape: even 2-adic valuation; shape: mixed valuations |
| 18 | `K2` | ⟨1, 3, sqrt(-21)⟩ | dim_a = 1; shape: valuation >= 2 at an odd prime |
| 19 | `K7` | ⟨4, 4*sqrt(5), 12⟩ | shape: even 2-adic valuation; shape: mixed valuations |
| 20 | `K8` | ⟨1, 2, sqrt(-5)⟩ | Algorithm 6; dim_a = 3 |
| 21 | `K6` | ⟨1, 1, 1⟩ | Algorithm 6; dim_a = 3; shape: all coefficients units |
| 22 | `K3` | ⟨1, 1, 1⟩ | Algorithm 7; Algorithm 8 retry; dim_a = 3 |
| 23 | `K7` | ⟨1, 1, 1⟩ | Algorithm 7 over an inert dyadic field |
| 24 | `K1` | ⟨-1, -1, -1⟩ | Algorithm 4 correction load-bearing at a dyadic prime |
| 25 | `K1` | ⟨sqrt(2), sqrt(2), sqrt(2)⟩ | Algorithm 4 correction load-bearing at a dyadic prime; odd dyadic valuations |
| 26 | `K5` | ⟨-1, -1, -1⟩ | Algorithm 4 correction load-bearing at a dyadic prime |
| 27 | `K5` | ⟨sqrt(10), 2, -1⟩ | Algorithm 4 correction load-bearing at a dyadic prime; mixed sign pattern (r = 2); Algorithm 8 retry |
| 28 | `K5` | ⟨1, -1, 1, -1⟩ | dim_a = 0 |
| 29 | `K1` | ⟨1, sqrt(2), 1, -1⟩ | Algorithm 8 entered with w = 1 mod 4; dim_a = 2 |
| 30 | `K2` | ⟨2, sqrt(-21), -1, 3⟩ | dim_a = 2; shape: even 2-adic valuation; shape: valuation >= 2 at an odd prime; shape: mixed valuations |
| 31 | `K6` | ⟨1, 1, 1, 1⟩ | Algorithm 5, non-real branch; dim_a = 4; Algorithm 8 entered with w = 2 mod 4; padding changed a dyadic Hasse invariant |
| 32 | `K2` | ⟨1, -2, -11, 22⟩ | Algorithm 5, non-real branch; dim_a = 4 |
| 33 | `K4` | ⟨1, -2, -11, 22⟩ | Algorithm 5, non-real branch; dim_a = 4 |
| 34 | `K8` | ⟨1, -2, -3, 6⟩ | Algorithm 5, non-real branch; dim_a = 4 |
| 35 | `K1` | ⟨1, 1, 1, 1⟩ | Algorithm 5, real branch; dim_a = 4; Algorithm 8 retry |
| 36 | `K2` | ⟨1, 1, 1, -1, 3, -3⟩ | Algorithm 8 entered with w = 2 mod 4; lambda_b in E_S basis |
| 37 | `K6` | ⟨1, sqrt(-7), 1, -1, 3, -3⟩ | Algorithm 8 entered with w = 2 mod 4; padding changed a dyadic Hasse invariant |
| 38 | `K1` | ⟨1, sqrt(2), 1, -1, 2, -2⟩ | Algorithm 8 entered with w = 2 mod 4 |
| 39 | `K1` | ⟨1, sqrt(2), 1, -1, 2, -2, 3, -3⟩ | Algorithm 8 entered with w = 3 mod 4; padding changed a dyadic Hasse invariant |
| 40 | `K1` | ⟨1, 1, 1, 1, 1⟩ | dim_a = 5; Algorithm 5 iterated; Algorithm 8 entered with w = 3 mod 4 |
| 41 | `K5` | ⟨-1, -1, -1, -1, -1⟩ | dim_a = 5; Algorithm 5 iterated; Algorithm 4 correction load-bearing at a dyadic prime; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4 |
| 42 | `K3` | ⟨1, 1, 1, 1, 1⟩ | dim_a = 5; Algorithm 5 iterated; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4 |
| 43 | `K1` | ⟨1, 1, 1, 1, 1, 1⟩ | dim_a = 6; Algorithm 5 iterated |
| 44 | `K3` | ⟨1, 1, 1, 1, 1, 1⟩ | dim_a = 6; Algorithm 5 iterated; Algorithm 8 retry |
| 45 | `K9` | ⟨-1, -1, -1⟩ | Algorithm 7; Algorithm 4 correction load-bearing at a dyadic prime; shape: all coefficients units |
| 46 | `K9` | ⟨a, a, a⟩ | mixed sign pattern (r = 3); Algorithm 4 correction load-bearing at a dyadic prime |
| 47 | `K9` | ⟨1, a, a^2 - 2⟩ | mixed sign pattern (r = 3) |
| 48 | `K10` | ⟨1, 1, 1⟩ | Algorithm 7; dim_a = 3 |
| 49 | `K9` | ⟨1, 1, 1, 1, 1, 1⟩ | dim_a = 6; Algorithm 5 iterated; A block with three rows |
| 50 | `K11` | ⟨1, -c, 2⟩ | Algorithm 7; dim_a = 3 |
| 51 | `K12` | ⟨-1, -1, -1⟩ | Algorithm 4 correction load-bearing at a dyadic prime; wild ramification (e = 4) |
| 52 | `K12` | ⟨w, 1, 1⟩ | mixed sign pattern (r = 4); wild ramification (e = 4) |
| 53 | `K12` | ⟨1, 1, 1⟩ | Algorithm 7; wild ramification (e = 4) |
| 54 | `K12` | ⟨1, -w⟩ | dim_a = 2; wild ramification (e = 4); A block nonempty |
| 55 | `K6` | ⟨-9*sqrt(-7) - 3, -1, -6*sqrt(-7) - 2, -sqrt(-7) + 1, 4*sqrt(-7) - 6, 2*sqrt(-7) - 3, -4*sqrt(-7) + 4⟩ | the paper's own worked example; Algorithm 6; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4 |
| 56 | `K13` | ⟨1, c10, c10^2 - 2, 3, -c10 - 1, 2*c10, c10^3 - 3*c10, 5, c10 - 2, -7⟩ | scale: degree 10 (FacElem stress); dim_a = 4; Algorithm 5, real branch |

## The measured coverage matrix

Measured at generation time by running every record through the
instrumented implementation (branch → number of covering records → which):

```
dim_a = 0                                                3  records 3,4,28
dim_a = 1                                                4  records 15,16,17,18
dim_a = 2                                               18  records 2,5,6,7,8,9,10,11,12,13,14,29,30,36,37,38,39,54
dim_a = 3                                               18  records 19,20,21,22,23,24,25,26,27,45,46,47,48,50,51,52,53,55
dim_a = 4                                                6  records 31,32,33,34,35,56
dim_a = 5                                                3  records 40,41,42
dim_a = 6                                                3  records 43,44,49
Algorithm 5, real branch                                 8  records 35,40,41,42,43,44,49,56
Algorithm 5, non-real branch                             4  records 31,32,33,34
Algorithm 5 iterated (>= 2 runs)                         6  records 40,41,42,43,44,49
Algorithm 6 (dim 3, non-real)                            7  records 20,21,31,32,33,34,55
Algorithm 7 (dim 3, real)                               23  records 19,22,23,24,25,26,27,35,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56
Algorithm 4 called                                      23  records 19,22,23,24,25,26,27,35,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56
Algorithm 4 correction load-bearing at a dyadic prime    9  records 24,25,26,27,41,45,46,51,56
Algorithm 4 correction load-bearing somewhere in S       9  records 24,25,26,27,41,45,46,51,56
Algorithm 8 reached                                     48  records 2,5,6,7,8,9,10,11,12,13,14,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56
Algorithm 8: lambda_b in the E_S basis                   8  records 6,7,8,9,10,11,12,36
Algorithm 8: lambda_b AND a nonempty A block             3  records 10,11,12
Algorithm 8: A block empty (real half inert)            17  records 5,6,7,8,9,13,20,21,27,30,31,32,33,34,36,37,55
Algorithm 8: A block nonempty                           31  records 2,10,11,12,14,19,22,23,24,25,26,29,35,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,56
Algorithm 8: retry loop enlarged S                       9  records 22,27,35,40,41,42,43,44,55
Algorithm 8 entered with w = 0 mod 4 (no padding)       15  records 2,5,6,7,8,9,10,11,12,13,14,43,44,49,54
Algorithm 8 entered with w = 1 mod 4                    20  records 19,20,21,22,23,24,25,26,27,29,30,45,46,47,48,50,51,52,53,56
Algorithm 8 entered with w = 2 mod 4                     8  records 31,32,33,34,35,36,37,38
Algorithm 8 entered with w = 3 mod 4                     5  records 39,40,41,42,55
Algorithm 8 padding changed a dyadic Hasse invariant    18  records 19,20,22,24,25,26,27,29,31,37,39,40,41,42,46,50,55,56
dimension correction decisive (blind-spot witness)       3  records 7,8,9
shape: all coefficients units                           26  records 3,6,7,9,11,12,21,22,23,24,26,28,31,35,40,41,42,43,44,45,46,47,48,49,51,53
shape: valuation >= 2 at an odd prime                    4  records 18,30,36,56
shape: even 2-adic valuation >= 2                       13  records 8,17,19,20,27,30,32,33,34,38,39,50,55
shape: mixed odd/dyadic valuations                      13  records 10,14,17,19,20,27,30,32,33,34,39,55,56
```

## Coverage limits — what stays below three records, and why

- **"dimension correction decisive" is single-field** (three forms, all over
  Q(sqrt34)). Witness existence is a property of the field, not the form: a
  witness is a totally positive non-square with even valuation everywhere
  that is a local square at every dyadic prime. A scan over 17 fields (the
  corpus fields plus real quadratic Q(sqrt d) for
  d = 3, 6, 7, 11, 14, 15, 19, 21, 22, 23, 33, 46, searching products of
  unit-group generators and 2·units) found witnesses over Q(sqrt34) ONLY.
  For h = 1 real quadratic fields with N(epsilon) = -1 — Q(sqrt2), Q(sqrt5)
  among them — non-existence is provable: every totally positive unit is a
  square. Over Q(sqrt34) the witness class is unique
  (2·(35+6·sqrt34) = (6+sqrt34)^2, so 2 and the fundamental unit are the
  same class); the three records use different literal coefficients of that
  one class.
- **"lambda_b AND a nonempty A block" is single-field** (three forms, all
  over Q(sqrt34)). This needs a formally real field whose class group keeps
  2-torsion after quotienting by the classes of S ∋ dyadic primes. Q(sqrt10)
  was proposed for this but VERIFIED NOT TO WORK: its dyadic prime is
  non-principal and generates Cl = Z/2, so C_S is trivial for every
  admissible S — the same mechanism as the recorded Q(sqrt-5) divergence.
  Q(sqrt34) works because its dyadic prime is principal (6+sqrt34 has
  norm 2). Q(sqrt65) (2 split) was also checked: both dyadic classes are
  non-principal and generate Cl.
- **Algorithm 5's non-real branch never iterates**: over a non-real field
  dim_a <= 4, so the branch runs at most once per input — a structural
  fact, not a coverage gap. Its records all enter at dim_a = 4 exactly.
- **"Algorithm 4 correction load-bearing at an ODD prime only"** has no
  record: in every measured instance the non-unit valuation of alpha_1
  appeared at dyadic primes (the t^2 integralisation scaling concentrates
  there). The dyadic case is the prime-2 hazard the correction exists for;
  an odd-prime-only instance was not constructed.
- **The calibration record** covers nothing in the matrix by design: it
  never enters the number-field layer.
