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
// dim_a = 0                                                3  records 3,4,28
// dim_a = 1                                                4  records 15,16,17,18
// dim_a = 2                                               18  records 2,5,6,7,8,9,10,11,12,13,14,29,30,36,37,38,39,54
// dim_a = 3                                               18  records 19,20,21,22,23,24,25,26,27,45,46,47,48,50,51,52,53,55
// dim_a = 4                                                6  records 31,32,33,34,35,56
// dim_a = 5                                                3  records 40,41,42
// dim_a = 6                                                3  records 43,44,49
// Algorithm 5, real branch                                 8  records 35,40,41,42,43,44,49,56
// Algorithm 5, non-real branch                             4  records 31,32,33,34
// Algorithm 5 iterated (>= 2 runs)                         6  records 40,41,42,43,44,49
// Algorithm 6 (dim 3, non-real)                            7  records 20,21,31,32,33,34,55
// Algorithm 7 (dim 3, real)                               23  records 19,22,23,24,25,26,27,35,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56
// Algorithm 4 called                                      23  records 19,22,23,24,25,26,27,35,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56
// Algorithm 4 correction load-bearing at a dyadic prime    9  records 24,25,26,27,41,45,46,51,56
// Algorithm 4 correction load-bearing somewhere in S       9  records 24,25,26,27,41,45,46,51,56
// Algorithm 8 reached                                     48  records 2,5,6,7,8,9,10,11,12,13,14,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56
// Algorithm 8: lambda_b in the E_S basis                   8  records 6,7,8,9,10,11,12,36
// Algorithm 8: lambda_b AND a nonempty A block             3  records 10,11,12
// Algorithm 8: A block empty (real half inert)            17  records 5,6,7,8,9,13,20,21,27,30,31,32,33,34,36,37,55
// Algorithm 8: A block nonempty                           31  records 2,10,11,12,14,19,22,23,24,25,26,29,35,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,56
// Algorithm 8: retry loop enlarged S                       9  records 22,27,35,40,41,42,43,44,55
// Algorithm 8 entered with w = 0 mod 4 (no padding)       15  records 2,5,6,7,8,9,10,11,12,13,14,43,44,49,54
// Algorithm 8 entered with w = 1 mod 4                    20  records 19,20,21,22,23,24,25,26,27,29,30,45,46,47,48,50,51,52,53,56
// Algorithm 8 entered with w = 2 mod 4                     8  records 31,32,33,34,35,36,37,38
// Algorithm 8 entered with w = 3 mod 4                     5  records 39,40,41,42,55
// Algorithm 8 padding changed a dyadic Hasse invariant    18  records 19,20,22,24,25,26,27,29,31,37,39,40,41,42,46,50,55,56
// dimension correction decisive (blind-spot witness)       3  records 7,8,9
// shape: all coefficients units                           26  records 3,6,7,9,11,12,21,22,23,24,26,28,31,35,40,41,42,43,44,45,46,47,48,49,51,53
// shape: valuation >= 2 at an odd prime                    4  records 18,30,36,56
// shape: even 2-adic valuation >= 2                       13  records 8,17,19,20,27,30,32,33,34,38,39,50,55
// shape: mixed odd/dyadic valuations                      13  records 10,14,17,19,20,27,30,32,33,34,39,55,56
// ============================================================================

SetColumns(0);
v1, v2, v3 := GetVersion();
printf "MAGMA-ORACLE-FORMAT 1\n";
printf "magma_version: %o.%o-%o\n", v1, v2, v3;
printf "run_date: ";
System("date");
CQFAvailable := true;
try
  _ := eval "AnisotropicPart";
catch e
  CQFAvailable := false;
end try;
printf "cqf_available: %o\n", CQFAvailable;
if not CQFAvailable then
  printf "note: Tier 2 (the anisotropic part) will be skipped -- CQF is not attached. Tier 1 output is complete and useful on its own.\n";
end if;

// Tier 1: plain Magma. Prints, per form: determinant; per real place the
// numerical embedding of the generator and the sign vector of the diagonal;
// per finite place over the listed rational primes its (p, e, f, gens) and
// the Hilbert-symbol product; global isotropy and Witt index.
procedure Tier1(id, label, K, diag, primes)
  printf "=== RECORD %o (tier 1) ===\n", id;
  printf "label: %o\n", label;
  printf "field_poly [Magma: DefiningPolynomial]: %o\n", DefiningPolynomial(K);
  printf "diag: %o\n", diag;
  n := #diag;
  try
    M := DiagonalMatrix(K, diag);
    printf "determinant [Magma: Determinant]: %o\n", Determinant(M);
    for v in RealPlaces(K) do
      printf "real_place theta_approx [Magma: Evaluate(K.1, place)]: %o\n", Evaluate(K.1, v);
      printf "real_place sign_vector [Magma: Sign(Evaluate(d_i, place))]: %o\n", [Sign(Evaluate(d, v)) : d in diag];
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
        printf "prime p=%o e=%o f=%o gens=%o hasse=prod_(i<j)HilbertSymbol(d_i,d_j,P) [Magma: HilbertSymbol]: %o\n",
          p, RamificationIndex(P), InertiaDegree(P), Generators(P), h;
      end for;
    end for;
    V := QuadraticSpace(M);
    printf "is_isotropic [Magma: IsIsotropic]: %o\n", IsIsotropic(V);
    printf "witt_index [Magma: WittIndex]: %o\n", WittIndex(V);
  catch e
    printf "RECORD-ERROR (tier 1, record %o): %o\n", id, e`Object;
  end try;
  printf "=== END RECORD %o (tier 1) ===\n", id;
end procedure;

// Tier 2: CQF only. The anisotropic part, printed raw. CQF's calling
// convention could not be tested on our side (no Magma licence); if the
// call errors, the message is printed and the run continues -- and we would
// be grateful for a one-line fix matching the CQF documentation.
procedure Tier2(id, label, K, diag, primes)
  printf "=== RECORD %o (tier 2) ===\n", id;
  try
    M := DiagonalMatrix(K, diag);
    ap := (eval "AnisotropicPart")(M);
    printf "anisotropic_part [CQF: AnisotropicPart(gram matrix)]: %o\n", ap;
  catch e
    try
      ap2 := (eval "AnisotropicPart")(diag);
      printf "anisotropic_part [CQF: AnisotropicPart(diagonal sequence)]: %o\n", ap2;
    catch e2
      printf "RECORD-ERROR (tier 2, record %o): %o\n", id, e2`Object;
    end try;
  end try;
  printf "=== END RECORD %o (tier 2) ===\n", id;
end procedure;

QQx<x> := PolynomialRing(Rationals());
K0 := RationalsAsNumberField();
K1<t1> := NumberField((-2)*x^0 + (1)*x^2);
K2<t2> := NumberField((21)*x^0 + (1)*x^2);
K3<t3> := NumberField((-34)*x^0 + (1)*x^2);
K4<t4> := NumberField((2)*x^0 + (1)*x^2);
K5<t5> := NumberField((-10)*x^0 + (1)*x^2);
K6<t6> := NumberField((7)*x^0 + (1)*x^2);
K7<t7> := NumberField((-5)*x^0 + (1)*x^2);
K8<t8> := NumberField((5)*x^0 + (1)*x^2);
K9<t9> := NumberField((-1)*x^0 + (-3)*x^1 + (1)*x^3);
K10<t10> := NumberField((-2)*x^0 + (1)*x^3);
K11<t11> := NumberField((1)*x^0 + (-4)*x^2 + (1)*x^4);
K12<t12> := NumberField((2)*x^0 + (-4)*x^2 + (1)*x^4);
K13<t13> := NumberField((-1)*x^0 + (5)*x^1 + (25)*x^2 + (-5)*x^3 + (-50)*x^4 + (1)*x^5 + (35)*x^6 + (-10)*x^8 + (1)*x^10);

records := [* *];
// record 1 chosen for: calibration
Append(~records, [* 1, "CALIBRATION over Q: <1,1,-3,-3,1,-1>, Macaulay2-validated. Not here for coverage - a canary: disagreement on THIS record means conventions or parsing, not mathematics.", K0, [K0 | (1), (1), (-3), (-3), (1), (-1)], [2, 3] *]);
// record 2 chosen for: dim_a = 2 over the ramified-dyadic real field; A block nonempty
Append(~records, [* 2, "Q(sqrt2), rank 2: <1, sqrt2>", K1, [K1 | (1), (1)*t1^1], [2] *]);
// record 3 chosen for: dim_a = 0
Append(~records, [* 3, "Q(sqrt2), rank 2 hyperbolic: <1, -1>", K1, [K1 | (1), (-1)], [2] *]);
// record 4 chosen for: dim_a = 0
Append(~records, [* 4, "Q(sqrt-21), rank 2 hyperbolic: <s, -s>", K2, [K2 | (1)*t2^1, (-1)*t2^1], [2, 3, 7] *]);
// record 5 chosen for: dim_a = 2, non-real; A block empty; lambda_b suppressed by S-classes (negative control)
Append(~records, [* 5, "Q(sqrt-21), rank 2: <1, s>. S = {P2,P3,P7} generates Cl = (Z/2)^2, so lambda_b does NOT fire here - contrast with <1,1>", K2, [K2 | (1), (1)*t2^1], [2, 3, 7] *]);
// record 6 chosen for: lambda_b in E_S basis; shape: all coefficients units
Append(~records, [* 6, "Q(sqrt-21), rank 2: <1, 1> - unit coefficients keep S dyadic-only, lambda_b fires", K2, [K2 | (1), (1)], [2] *]);
// record 7 chosen for: dimension correction decisive; lambda_b in E_S basis; shape: all coefficients units
Append(~records, [* 7, "Q(sqrt34), rank 2: the blind-spot form <1, -(35+6*sqrt34)>", K3, [K3 | (1), (-35) + (-6)*t3^1], [2] *]);
// record 8 chosen for: dimension correction decisive; lambda_b in E_S basis
Append(~records, [* 8, "Q(sqrt34), rank 2: <2, -1> - second blind-spot witness (2 ~ eps34 mod squares)", K3, [K3 | (2), (-1)], [2] *]);
// record 9 chosen for: dimension correction decisive; lambda_b in E_S basis
Append(~records, [* 9, "Q(sqrt34), rank 2: <35+6*sqrt34, -1> - third blind-spot witness", K3, [K3 | (35) + (6)*t3^1, (-1)], [2] *]);
// record 10 chosen for: lambda_b AND a nonempty A block
Append(~records, [* 10, "Q(sqrt34), rank 2: <1, -sqrt34> - lambda_b together with a real-place condition", K3, [K3 | (1), (-1)*t3^1], [2, 17] *]);
// record 11 chosen for: lambda_b AND a nonempty A block; shape: all coefficients units
Append(~records, [* 11, "Q(sqrt34), rank 2: <1, 1> - lambda_b, A block from both real places", K3, [K3 | (1), (1)], [2] *]);
// record 12 chosen for: lambda_b AND a nonempty A block; correction non-decisive contrast
Append(~records, [* 12, "Q(sqrt34), rank 2: <-1, -(35+6*sqrt34)> - witness class but detected at the real places (correction NOT decisive)", K3, [K3 | (-1), (-35) + (-6)*t3^1], [2] *]);
// record 13 chosen for: dim_a = 2, ramified dyadic without real places; A block empty
Append(~records, [* 13, "Q(sqrt-2), rank 2: <1, sqrt-2> - ramified dyadic, non-real control", K4, [K4 | (1), (1)*t4^1], [2] *]);
// record 14 chosen for: dim_a = 2 over the triple-point field; A block nonempty
Append(~records, [* 14, "Q(sqrt10), rank 2: <1, -sqrt10>", K5, [K5 | (1), (-1)*t5^1], [2, 5] *]);
// record 15 chosen for: dim_a = 1
Append(~records, [* 15, "Q(sqrt-7), rank 3: <1, 2, -2> - hyperbolic pair plus a unit", K6, [K6 | (1), (2), (-2)], [2] *]);
// record 16 chosen for: dim_a = 1
Append(~records, [* 16, "Q(sqrt-2), rank 3: <3, 1, -3>", K4, [K4 | (3), (1), (-3)], [2, 3] *]);
// record 17 chosen for: dim_a = 1; shape: even 2-adic valuation; shape: mixed valuations
Append(~records, [* 17, "Q(sqrt34), rank 3", K3, [K3 | (1)*t3^1, (2), (-1)], [2, 17] *]);
// record 18 chosen for: dim_a = 1; shape: valuation >= 2 at an odd prime
Append(~records, [* 18, "Q(sqrt-21), rank 3", K2, [K2 | (1), (3), (1)*t2^1], [2, 3, 7] *]);
// record 19 chosen for: shape: even 2-adic valuation; shape: mixed valuations
Append(~records, [* 19, "Q(sqrt5), rank 3, even dyadic valuations: <4, 4*sqrt5, 12>", K7, [K7 | (4), (4)*t7^1, (12)], [2, 3, 5] *]);
// record 20 chosen for: Algorithm 6; dim_a = 3
Append(~records, [* 20, "Q(sqrt-5), rank 3, class number 2", K8, [K8 | (1), (2), (1)*t8^1], [2, 5] *]);
// record 21 chosen for: Algorithm 6; dim_a = 3; shape: all coefficients units
Append(~records, [* 21, "Q(sqrt-7), rank 3, definite entries (anisotropic only at the dyadic places)", K6, [K6 | (1), (1), (1)], [2] *]);
// record 22 chosen for: Algorithm 7; Algorithm 8 retry; dim_a = 3
Append(~records, [* 22, "Q(sqrt34), rank 3, definite", K3, [K3 | (1), (1), (1)], [2] *]);
// record 23 chosen for: Algorithm 7 over an inert dyadic field
Append(~records, [* 23, "Q(sqrt5), rank 3, definite - inert dyadic prime", K7, [K7 | (1), (1), (1)], [2] *]);
// record 24 chosen for: Algorithm 4 correction load-bearing at a dyadic prime
Append(~records, [* 24, "Q(sqrt2), rank 3, negative definite: <-1, -1, -1>", K1, [K1 | (-1), (-1), (-1)], [2] *]);
// record 25 chosen for: Algorithm 4 correction load-bearing at a dyadic prime; odd dyadic valuations
Append(~records, [* 25, "Q(sqrt2), rank 3: <sqrt2, sqrt2, sqrt2>", K1, [K1 | (1)*t1^1, (1)*t1^1, (1)*t1^1], [2] *]);
// record 26 chosen for: Algorithm 4 correction load-bearing at a dyadic prime
Append(~records, [* 26, "Q(sqrt10), rank 3, negative definite: <-1, -1, -1>", K5, [K5 | (-1), (-1), (-1)], [2] *]);
// record 27 chosen for: Algorithm 4 correction load-bearing at a dyadic prime; mixed sign pattern (r = 2); Algorithm 8 retry
Append(~records, [* 27, "Q(sqrt10), rank 3, mixed signature: <sqrt10, 2, -1>", K5, [K5 | (1)*t5^1, (2), (-1)], [2, 5] *]);
// record 28 chosen for: dim_a = 0
Append(~records, [* 28, "Q(sqrt10), rank 4 hyperbolic: <1, -1, 1, -1>", K5, [K5 | (1), (-1), (1), (-1)], [2] *]);
// record 29 chosen for: Algorithm 8 entered with w = 1 mod 4; dim_a = 2
Append(~records, [* 29, "Q(sqrt2), rank 4, Witt-padding case: <1, sqrt2, 1, -1>", K1, [K1 | (1), (1)*t1^1, (1), (-1)], [2] *]);
// record 30 chosen for: dim_a = 2; shape: even 2-adic valuation; shape: valuation >= 2 at an odd prime; shape: mixed valuations
Append(~records, [* 30, "Q(sqrt-21), rank 4", K2, [K2 | (2), (1)*t2^1, (-1), (3)], [2, 3, 7] *]);
// record 31 chosen for: Algorithm 5, non-real branch; dim_a = 4; Algorithm 8 entered with w = 2 mod 4; padding changed a dyadic Hasse invariant
Append(~records, [* 31, "Q(sqrt-7), rank 4, definite entries: <1, 1, 1, 1> - anisotropic at the split dyadic primes", K6, [K6 | (1), (1), (1), (1)], [2] *]);
// record 32 chosen for: Algorithm 5, non-real branch; dim_a = 4
Append(~records, [* 32, "Q(sqrt-21), rank 4: <1, -2, -11, 22> - quaternion norm form, anisotropic at a place over the split prime 11", K2, [K2 | (1), (-2), (-11), (22)], [2, 11] *]);
// record 33 chosen for: Algorithm 5, non-real branch; dim_a = 4
Append(~records, [* 33, "Q(sqrt-2), rank 4: <1, -2, -11, 22> - quaternion norm form, anisotropic at a place over the split prime 11", K4, [K4 | (1), (-2), (-11), (22)], [2, 11] *]);
// record 34 chosen for: Algorithm 5, non-real branch; dim_a = 4
Append(~records, [* 34, "Q(sqrt-5), rank 4: <1, -2, -3, 6> - quaternion norm form, anisotropic at a place over the split prime 3", K8, [K8 | (1), (-2), (-3), (6)], [2, 3] *]);
// record 35 chosen for: Algorithm 5, real branch; dim_a = 4; Algorithm 8 retry
Append(~records, [* 35, "Q(sqrt2), rank 4, definite: <1, 1, 1, 1>", K1, [K1 | (1), (1), (1), (1)], [2] *]);
// record 36 chosen for: Algorithm 8 entered with w = 2 mod 4; lambda_b in E_S basis
Append(~records, [* 36, "Q(sqrt-21), rank 6: <1, 1> + 2H - lambda_b meeting w = 2 mod 4", K2, [K2 | (1), (1), (1), (-1), (3), (-3)], [2] *]);
// record 37 chosen for: Algorithm 8 entered with w = 2 mod 4; padding changed a dyadic Hasse invariant
Append(~records, [* 37, "Q(sqrt-7), rank 6: <1, sqrt-7> + 2H", K6, [K6 | (1), (1)*t6^1, (1), (-1), (3), (-3)], [2, 3, 7] *]);
// record 38 chosen for: Algorithm 8 entered with w = 2 mod 4
Append(~records, [* 38, "Q(sqrt2), rank 6: <1, sqrt2> + 2H", K1, [K1 | (1), (1)*t1^1, (1), (-1), (2), (-2)], [2] *]);
// record 39 chosen for: Algorithm 8 entered with w = 3 mod 4; padding changed a dyadic Hasse invariant
Append(~records, [* 39, "Q(sqrt2), rank 8: <1, sqrt2> + 3H", K1, [K1 | (1), (1)*t1^1, (1), (-1), (2), (-2), (3), (-3)], [2, 3] *]);
// record 40 chosen for: dim_a = 5; Algorithm 5 iterated; Algorithm 8 entered with w = 3 mod 4
Append(~records, [* 40, "Q(sqrt2), rank 5, definite: <1, 1, 1, 1, 1>", K1, [K1 | (1), (1), (1), (1), (1)], [2] *]);
// record 41 chosen for: dim_a = 5; Algorithm 5 iterated; Algorithm 4 correction load-bearing at a dyadic prime; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4
Append(~records, [* 41, "Q(sqrt10), rank 5, negative definite: <-1, ..., -1>", K5, [K5 | (-1), (-1), (-1), (-1), (-1)], [2] *]);
// record 42 chosen for: dim_a = 5; Algorithm 5 iterated; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4
Append(~records, [* 42, "Q(sqrt34), rank 5, definite: <1, 1, 1, 1, 1>", K3, [K3 | (1), (1), (1), (1), (1)], [2] *]);
// record 43 chosen for: dim_a = 6; Algorithm 5 iterated
Append(~records, [* 43, "Q(sqrt2), rank 6, definite: <1, ..., 1>", K1, [K1 | (1), (1), (1), (1), (1), (1)], [2] *]);
// record 44 chosen for: dim_a = 6; Algorithm 5 iterated; Algorithm 8 retry
Append(~records, [* 44, "Q(sqrt34), rank 6, definite: <1, ..., 1>", K3, [K3 | (1), (1), (1), (1), (1), (1)], [2] *]);
// record 45 chosen for: Algorithm 7; Algorithm 4 correction load-bearing at a dyadic prime; shape: all coefficients units
Append(~records, [* 45, "Q[x]/(x^3-3x-1), rank 3, negative definite", K9, [K9 | (-1), (-1), (-1)], [2] *]);
// record 46 chosen for: mixed sign pattern (r = 3); Algorithm 4 correction load-bearing at a dyadic prime
Append(~records, [* 46, "Q[x]/(x^3-3x-1), rank 3: <a, a, a> - sign pattern (+,-,-) across the three real places", K9, [K9 | (1)*t9^1, (1)*t9^1, (1)*t9^1], [2] *]);
// record 47 chosen for: mixed sign pattern (r = 3)
Append(~records, [* 47, "Q[x]/(x^3-3x-1), rank 3: <1, a, a^2-2> - sign patterns differing at every real place", K9, [K9 | (1), (1)*t9^1, (-2) + (1)*t9^2], [2] *]);
// record 48 chosen for: Algorithm 7; dim_a = 3
Append(~records, [* 48, "Q[x]/(x^3-2), rank 3 (mixed-signature field)", K10, [K10 | (1), (1), (1)], [2] *]);
// record 49 chosen for: dim_a = 6; Algorithm 5 iterated; A block with three rows
Append(~records, [* 49, "Q[x]/(x^3-3x-1), rank 6, definite: <1, ..., 1>", K9, [K9 | (1), (1), (1), (1), (1), (1)], [2] *]);
// record 50 chosen for: Algorithm 7; dim_a = 3
Append(~records, [* 50, "quartic Q[x]/(x^4-4x^2+1), rank 3", K11, [K11 | (1), (-1)*t11^1, (2)], [2] *]);
// record 51 chosen for: Algorithm 4 correction load-bearing at a dyadic prime; wild ramification (e = 4)
Append(~records, [* 51, "Q(zeta16)+, rank 3, negative definite: <-1, -1, -1> - wild dyadic ramification e = 4", K12, [K12 | (-1), (-1), (-1)], [2] *]);
// record 52 chosen for: mixed sign pattern (r = 4); wild ramification (e = 4)
Append(~records, [* 52, "Q(zeta16)+, rank 3: <w, 1, 1> - mixed sign pattern over the wild dyadic field", K12, [K12 | (1)*t12^1, (1), (1)], [2] *]);
// record 53 chosen for: Algorithm 7; wild ramification (e = 4)
Append(~records, [* 53, "Q(zeta16)+, rank 3, definite", K12, [K12 | (1), (1), (1)], [2] *]);
// record 54 chosen for: dim_a = 2; wild ramification (e = 4); A block nonempty
Append(~records, [* 54, "Q(zeta16)+, rank 2: <1, -w>", K12, [K12 | (1), (-1)*t12^1], [2] *]);
// record 55 chosen for: the paper's own worked example; Algorithm 6; Algorithm 8 retry; Algorithm 8 entered with w = 3 mod 4
Append(~records, [* 55, "KR21 section-6 fixture, dimension 7 over Q(sqrt-7)", K6, [K6 | (-3) + (-9)*t6^1, (-1), (-2) + (-6)*t6^1, (1) + (-1)*t6^1, (-6) + (4)*t6^1, (-3) + (2)*t6^1, (4) + (-4)*t6^1], [2, 3, 37] *]);
// record 56 chosen for: scale: degree 10 (FacElem stress); dim_a = 4; Algorithm 5, real branch
Append(~records, [* 56, "dimension 10 over the degree-10 totally real field Q(2cos(2pi/25)) - the expensive record, deliberately last", K13, [K13 | (1), (1)*t13^1, (-2) + (1)*t13^2, (3), (-1) + (-1)*t13^1, (2)*t13^1, (-3)*t13^1 + (1)*t13^3, (5), (-2) + (1)*t13^1, (-7)], [2, 3, 5, 7] *]);

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
printf "MAGMA-ORACLE-DONE\n";
quit;
