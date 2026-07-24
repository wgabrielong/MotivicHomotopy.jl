# The differential corpus, shared by generate_oracle.jl (which turns it into
# oracle.m and ../CORPUS.md) and coverage.jl (which measures which internal
# branches each record exercises). Each record is
#     (label, K, entries, covers)
# where `covers` names the coverage-matrix rows the record was CHOSEN to
# exercise — the intent, recorded so that when a divergence is reported the
# question "what was this form for" has an answer without re-derivation. The
# matrix embedded in oracle.m / CORPUS.md is the MEASURED one (coverage.jl),
# not this declared one; where measurement contradicted an earlier intent the
# annotation below was corrected to the measurement.
#
# Ordered cheapest-first for the Magma run: the Q calibration record, then
# quadratic fields (rank 2, then 3, then >= 4), cubics, quartics, the KR
# section-6 fixture, and the degree-10 field last — a killed run still sends
# useful data. Record 2 is kept in position 2 so that results-stub.txt stays
# valid.

using Oscar
using MotivicHomotopy
const MH = MotivicHomotopy

include(joinpath(@__DIR__, "..", "..", "test", "nf_test_fields.jl"))

# The corpus fields come from test/nf_test_fields.jl (the same objects the
# test suite uses); _verify_corpus_fields below asserts every property the
# coverage argument rests on. Only the degree-10 scale field is defined
# here, matching the suite's scale probe.
Zx, zx = polynomial_ring(ZZ, "zx")
K10, c10 = number_field(map_coefficients(QQ, cos_minpoly(25, zx)), "c10")

eps34 = 35 + 6 * s34    # fundamental unit of Q(sqrt34): totally positive,
                        # NOT a square, and a local square at the dyadic
                        # prime — the anisotropic-dimension blind-spot
                        # witness. Note 2 * eps34 = (6 + s34)^2, so 2 and
                        # eps34 are the SAME square class: Q(sqrt34) has
                        # exactly one witness class.

# verified computationally at generation time — these are the properties the
# coverage argument rests on, so they are asserted, not trusted
function _verify_corpus_fields()
    for (K, h, e2, dyadic_principal) in
        [(Ks10, 2, 2, false), (K162, 1, 4, true), (K2i, 1, 2, true), (K34, 2, 2, true)]
        OK = maximal_order(K)
        C, mC = class_group(OK)
        order(C) == h || error("class number of $K is not $h")
        dec = prime_decomposition(OK, 2)
        maximum(ramification_index(P) for (P, _) in dec) == e2 ||
            error("dyadic ramification of $K is not e = $e2")
        all(is_principal(P) for (P, _) in dec) == dyadic_principal ||
            error("dyadic principality of $K is not $dyadic_principal")
    end
    # the K162 field really is the maximal real subfield of Q(zeta_16)
    # (coefficient-wise: the two polynomials live in different rings)
    let f = defining_polynomial(K162), g = cos_minpoly(16, zx)
        (degree(f) == degree(g) &&
         all(coeff(f, i) == coeff(g, i) for i in 0:degree(g))) ||
            error("K162 is not Q(zeta_16)^+")
    end
    # the blind-spot witness class over Q(sqrt34): eps34 ~ 2 mod squares
    is_square(2 * eps34) || error("2 and eps34 are not the same square class")
    is_square(eps34) && error("eps34 is a square")
    nothing
end

corpus = Any[
    # --- calibration -------------------------------------------------------
    ("CALIBRATION over Q: <1,1,-3,-3,1,-1>, Macaulay2-validated. Not here for coverage - a canary: disagreement on THIS record means conventions or parsing, not mathematics.",
     :Q, Any[1, 1, -3, -3, 1, -1],
     ["calibration"]),
    # --- quadratic fields, rank 2 ------------------------------------------
    ("Q(sqrt2), rank 2: <1, sqrt2>", K2, Any[K2(1), s2],
     ["dim_a = 2 over the ramified-dyadic real field", "A block nonempty"]),
    ("Q(sqrt2), rank 2 hyperbolic: <1, -1>", K2, Any[K2(1), K2(-1)],
     ["dim_a = 0"]),
    ("Q(sqrt-21), rank 2 hyperbolic: <s, -s>", K21, Any[s21, -s21],
     ["dim_a = 0"]),
    ("Q(sqrt-21), rank 2: <1, s>. S = {P2,P3,P7} generates Cl = (Z/2)^2, so lambda_b does NOT fire here - contrast with <1,1>",
     K21, Any[K21(1), s21],
     ["dim_a = 2, non-real", "A block empty", "lambda_b suppressed by S-classes (negative control)"]),
    ("Q(sqrt-21), rank 2: <1, 1> - unit coefficients keep S dyadic-only, lambda_b fires",
     K21, Any[K21(1), K21(1)],
     ["lambda_b in E_S basis", "shape: all coefficients units"]),
    ("Q(sqrt34), rank 2: the blind-spot form <1, -(35+6*sqrt34)>", K34, Any[K34(1), -eps34],
     ["dimension correction decisive", "lambda_b in E_S basis", "shape: all coefficients units"]),
    ("Q(sqrt34), rank 2: <2, -1> - second blind-spot witness (2 ~ eps34 mod squares)",
     K34, Any[K34(2), K34(-1)],
     ["dimension correction decisive", "lambda_b in E_S basis"]),
    ("Q(sqrt34), rank 2: <35+6*sqrt34, -1> - third blind-spot witness",
     K34, Any[eps34, K34(-1)],
     ["dimension correction decisive", "lambda_b in E_S basis"]),
    ("Q(sqrt34), rank 2: <1, -sqrt34> - lambda_b together with a real-place condition",
     K34, Any[K34(1), -s34],
     ["lambda_b AND a nonempty A block"]),
    ("Q(sqrt34), rank 2: <1, 1> - lambda_b, A block from both real places",
     K34, Any[K34(1), K34(1)],
     ["lambda_b AND a nonempty A block", "shape: all coefficients units"]),
    ("Q(sqrt34), rank 2: <-1, -(35+6*sqrt34)> - witness class but detected at the real places (correction NOT decisive)",
     K34, Any[K34(-1), -eps34],
     ["lambda_b AND a nonempty A block", "correction non-decisive contrast"]),
    ("Q(sqrt-2), rank 2: <1, sqrt-2> - ramified dyadic, non-real control", K2i, Any[K2i(1), s2i],
     ["dim_a = 2, ramified dyadic without real places", "A block empty"]),
    ("Q(sqrt10), rank 2: <1, -sqrt10>", Ks10, Any[Ks10(1), -s10],
     ["dim_a = 2 over the triple-point field", "A block nonempty"]),
    # --- quadratic fields, rank 3 ------------------------------------------
    ("Q(sqrt-7), rank 3: <1, 2, -2> - hyperbolic pair plus a unit", K7, Any[K7(1), K7(2), K7(-2)],
     ["dim_a = 1"]),
    ("Q(sqrt-2), rank 3: <3, 1, -3>", K2i, Any[K2i(3), K2i(1), K2i(-3)],
     ["dim_a = 1"]),
    ("Q(sqrt34), rank 3", K34, Any[s34, K34(2), K34(-1)],
     ["dim_a = 1", "shape: even 2-adic valuation", "shape: mixed valuations"]),
    ("Q(sqrt-21), rank 3", K21, Any[K21(1), K21(3), s21],
     ["dim_a = 1", "shape: valuation >= 2 at an odd prime"]),
    ("Q(sqrt5), rank 3, even dyadic valuations: <4, 4*sqrt5, 12>", Kr5, Any[Kr5(4), 4 * sr5, Kr5(12)],
     ["shape: even 2-adic valuation", "shape: mixed valuations"]),
    ("Q(sqrt-5), rank 3, class number 2", K5, Any[K5(1), K5(2), s5],
     ["Algorithm 6", "dim_a = 3"]),
    ("Q(sqrt-7), rank 3, definite entries (anisotropic only at the dyadic places)", K7, Any[K7(1), K7(1), K7(1)],
     ["Algorithm 6", "dim_a = 3", "shape: all coefficients units"]),
    ("Q(sqrt34), rank 3, definite", K34, Any[K34(1), K34(1), K34(1)],
     ["Algorithm 7", "Algorithm 8 retry", "dim_a = 3"]),
    ("Q(sqrt5), rank 3, definite - inert dyadic prime", Kr5, Any[Kr5(1), Kr5(1), Kr5(1)],
     ["Algorithm 7 over an inert dyadic field"]),
    ("Q(sqrt2), rank 3, negative definite: <-1, -1, -1>", K2, Any[K2(-1), K2(-1), K2(-1)],
     ["Algorithm 4 correction load-bearing at a dyadic prime"]),
    ("Q(sqrt2), rank 3: <sqrt2, sqrt2, sqrt2>", K2, Any[s2, s2, s2],
     ["Algorithm 4 correction load-bearing at a dyadic prime", "odd dyadic valuations"]),
    ("Q(sqrt10), rank 3, negative definite: <-1, -1, -1>", Ks10, Any[Ks10(-1), Ks10(-1), Ks10(-1)],
     ["Algorithm 4 correction load-bearing at a dyadic prime"]),
    ("Q(sqrt10), rank 3, mixed signature: <sqrt10, 2, -1>", Ks10, Any[s10, Ks10(2), Ks10(-1)],
     ["Algorithm 4 correction load-bearing at a dyadic prime", "mixed sign pattern (r = 2)", "Algorithm 8 retry"]),
    # --- quadratic fields, rank >= 4 ---------------------------------------
    ("Q(sqrt10), rank 4 hyperbolic: <1, -1, 1, -1>", Ks10, Any[Ks10(1), Ks10(-1), Ks10(1), Ks10(-1)],
     ["dim_a = 0"]),
    ("Q(sqrt2), rank 4, Witt-padding case: <1, sqrt2, 1, -1>", K2, Any[K2(1), s2, K2(1), K2(-1)],
     ["Algorithm 8 entered with w = 1 mod 4", "dim_a = 2"]),
    ("Q(sqrt-21), rank 4", K21, Any[K21(2), s21, K21(-1), K21(3)],
     ["dim_a = 2", "shape: even 2-adic valuation", "shape: valuation >= 2 at an odd prime", "shape: mixed valuations"]),
    ("Q(sqrt-7), rank 4, definite entries: <1, 1, 1, 1> - anisotropic at the split dyadic primes",
     K7, Any[K7(1), K7(1), K7(1), K7(1)],
     ["Algorithm 5, non-real branch", "dim_a = 4", "Algorithm 8 entered with w = 2 mod 4", "padding changed a dyadic Hasse invariant"]),
    ("Q(sqrt-21), rank 4: <1, -2, -11, 22> - quaternion norm form, anisotropic at a place over the split prime 11",
     K21, Any[K21(1), K21(-2), K21(-11), K21(22)],
     ["Algorithm 5, non-real branch", "dim_a = 4"]),
    ("Q(sqrt-2), rank 4: <1, -2, -11, 22> - quaternion norm form, anisotropic at a place over the split prime 11",
     K2i, Any[K2i(1), K2i(-2), K2i(-11), K2i(22)],
     ["Algorithm 5, non-real branch", "dim_a = 4"]),
    ("Q(sqrt-5), rank 4: <1, -2, -3, 6> - quaternion norm form, anisotropic at a place over the split prime 3",
     K5, Any[K5(1), K5(-2), K5(-3), K5(6)],
     ["Algorithm 5, non-real branch", "dim_a = 4"]),
    ("Q(sqrt2), rank 4, definite: <1, 1, 1, 1>", K2, Any[K2(1), K2(1), K2(1), K2(1)],
     ["Algorithm 5, real branch", "dim_a = 4", "Algorithm 8 retry"]),
    ("Q(sqrt-21), rank 6: <1, 1> + 2H - lambda_b meeting w = 2 mod 4",
     K21, Any[K21(1), K21(1), K21(1), K21(-1), K21(3), K21(-3)],
     ["Algorithm 8 entered with w = 2 mod 4", "lambda_b in E_S basis"]),
    ("Q(sqrt-7), rank 6: <1, sqrt-7> + 2H", K7, Any[K7(1), s7, K7(1), K7(-1), K7(3), K7(-3)],
     ["Algorithm 8 entered with w = 2 mod 4", "padding changed a dyadic Hasse invariant"]),
    ("Q(sqrt2), rank 6: <1, sqrt2> + 2H", K2, Any[K2(1), s2, K2(1), K2(-1), K2(2), K2(-2)],
     ["Algorithm 8 entered with w = 2 mod 4"]),
    ("Q(sqrt2), rank 8: <1, sqrt2> + 3H", K2,
     Any[K2(1), s2, K2(1), K2(-1), K2(2), K2(-2), K2(3), K2(-3)],
     ["Algorithm 8 entered with w = 3 mod 4", "padding changed a dyadic Hasse invariant"]),
    ("Q(sqrt2), rank 5, definite: <1, 1, 1, 1, 1>", K2, Any[fill(K2(1), 5)...],
     ["dim_a = 5", "Algorithm 5 iterated", "Algorithm 8 entered with w = 3 mod 4"]),
    ("Q(sqrt10), rank 5, negative definite: <-1, ..., -1>", Ks10, Any[fill(Ks10(-1), 5)...],
     ["dim_a = 5", "Algorithm 5 iterated", "Algorithm 4 correction load-bearing at a dyadic prime", "Algorithm 8 retry", "Algorithm 8 entered with w = 3 mod 4"]),
    ("Q(sqrt34), rank 5, definite: <1, 1, 1, 1, 1>", K34, Any[fill(K34(1), 5)...],
     ["dim_a = 5", "Algorithm 5 iterated", "Algorithm 8 retry", "Algorithm 8 entered with w = 3 mod 4"]),
    ("Q(sqrt2), rank 6, definite: <1, ..., 1>", K2, Any[fill(K2(1), 6)...],
     ["dim_a = 6", "Algorithm 5 iterated"]),
    ("Q(sqrt34), rank 6, definite: <1, ..., 1>", K34, Any[fill(K34(1), 6)...],
     ["dim_a = 6", "Algorithm 5 iterated", "Algorithm 8 retry"]),
    # --- cubic fields -------------------------------------------------------
    ("Q[x]/(x^3-3x-1), rank 3, negative definite", Kc, Any[Kc(-1), Kc(-1), Kc(-1)],
     ["Algorithm 7", "Algorithm 4 correction load-bearing at a dyadic prime", "shape: all coefficients units"]),
    ("Q[x]/(x^3-3x-1), rank 3: <a, a, a> - sign pattern (+,-,-) across the three real places",
     Kc, Any[ac, ac, ac],
     ["mixed sign pattern (r = 3)", "Algorithm 4 correction load-bearing at a dyadic prime"]),
    ("Q[x]/(x^3-3x-1), rank 3: <1, a, a^2-2> - sign patterns differing at every real place",
     Kc, Any[Kc(1), ac, ac^2 - 2],
     ["mixed sign pattern (r = 3)"]),
    ("Q[x]/(x^3-2), rank 3 (mixed-signature field)", K32, Any[K32(1), K32(1), K32(1)],
     ["Algorithm 7", "dim_a = 3"]),
    ("Q[x]/(x^3-3x-1), rank 6, definite: <1, ..., 1>", Kc, Any[fill(Kc(1), 6)...],
     ["dim_a = 6", "Algorithm 5 iterated", "A block with three rows"]),
    # --- quartic fields -----------------------------------------------------
    ("quartic Q[x]/(x^4-4x^2+1), rank 3", Kq, Any[Kq(1), -aq, Kq(2)],
     ["Algorithm 7", "dim_a = 3"]),
    ("Q(zeta16)+, rank 3, negative definite: <-1, -1, -1> - wild dyadic ramification e = 4",
     K162, Any[K162(-1), K162(-1), K162(-1)],
     ["Algorithm 4 correction load-bearing at a dyadic prime", "wild ramification (e = 4)"]),
    ("Q(zeta16)+, rank 3: <w, 1, 1> - mixed sign pattern over the wild dyadic field",
     K162, Any[w16, K162(1), K162(1)],
     ["mixed sign pattern (r = 4)", "wild ramification (e = 4)"]),
    ("Q(zeta16)+, rank 3, definite", K162, Any[K162(1), K162(1), K162(1)],
     ["Algorithm 7", "wild ramification (e = 4)"]),
    ("Q(zeta16)+, rank 2: <1, -w>", K162, Any[K162(1), -w16],
     ["dim_a = 2", "wild ramification (e = 4)", "A block nonempty"]),
    # --- the expensive tail -------------------------------------------------
    ("KR21 section-6 fixture, dimension 7 over Q(sqrt-7)", K7,
     Any[-3 - 9s7, K7(-1), -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7],
     ["the paper's own worked example", "Algorithm 6", "Algorithm 8 retry", "Algorithm 8 entered with w = 3 mod 4"]),
    ("dimension 10 over the degree-10 totally real field Q(2cos(2pi/25)) - the expensive record, deliberately last",
     K10, Any[K10(1), c10, c10^2 - 2, K10(3), -c10 - 1, 2c10, c10^3 - 3c10, K10(5), c10 - 2, K10(-7)],
     ["scale: degree 10 (FacElem stress)", "dim_a = 4", "Algorithm 5, real branch"]),
]

# why each field is in the corpus — consumed by generate_oracle.jl for the
# CORPUS.md field table
field_reasons = Any[
    K2 => "totally real, 2 ramified, h = 1: the test suite's Algorithm 4 correction field",
    K21 => "imaginary, Cl = (Z/2)^2, 2 ramified non-principal: lambda_b fires with unit coefficients, is suppressed once S grows",
    K34 => "totally real, h = 2, dyadic prime PRINCIPAL: the only field here where lambda_b coexists with real places, and the only scanned field with a blind-spot witness class",
    K2i => "imaginary, 2 ramified, h = 1: non-real control isolating dyadic ramification from formal reality",
    Ks10 => "totally real, 2 ramified, h = 2: Algorithm 4 dyadic correction + mixed signs; lambda_b provably never fires here (see the coverage-limits section below)",
    K7 => "imaginary, 2 split, h = 1: KR's own example field; completions at the dyadic places are Q_2 itself",
    Kr5 => "totally real, 2 inert, h = 1",
    K5 => "imaginary, h = 2, dyadic class generates Cl: the recorded Algorithm 1 divergence field",
    Kc => "totally real cubic, three real places: genuinely non-constant sign patterns",
    K32 => "mixed signature (one real, one complex place)",
    Kq => "totally real quartic",
    K162 => "Q(zeta_16)^+: 2 totally ramified, e = 4 — wild dyadic ramification, the deepest prime-2 corner",
    K10 => "degree 10, totally real: the FacElem scale stress, one expensive record only",
]

# the explicitly-argued exceptions to the >= 3 covering-forms rule — consumed
# by generate_oracle.jl for CORPUS.md
coverage_limits_md = """
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
"""
