# Convention gates: every Hecke intrinsic used as an oracle is pinned on
# hand-known values before any suite runs; a gate failure aborts the run.
# Each expected value is derivable by hand, with the derivation stated at
# the gate.

function run_gates()
    QQx, x = polynomial_ring(QQ, "x")
    KQ, _ = Hecke.rationals_as_number_field()
    K2, s2 = quadratic_field(2)
    Kr5, _ = quadratic_field(5)
    K7, _ = quadratic_field(-7)
    K34, s34 = quadratic_field(34)

    gate(ok, what) = ok || error("Hecke convention gate failed: $what")

    # --- is_isometric, global ------------------------------------------------
    # <1,1> ~ <2,2> over Q: 2 = 1^2 + 1^2 is represented, same rank/disc.
    gate(oracle_isometric(KQ, [KQ(1), KQ(1)], [KQ(2), KQ(2)]),
         "is_isometric: <1,1> ~ <2,2> over Q")
    # <1,1> is definite at both real places of Q(sqrt2); <1,-1> is not.
    gate(!oracle_isometric(K2, [K2(1), K2(1)], [K2(1), K2(-1)]),
         "is_isometric: <1,1> vs <1,-1> over Q(sqrt2)")
    # 2 = (sqrt2)^2 is a square in Q(sqrt2), so <2,2> ~ <1,1> there.
    gate(oracle_isometric(K2, [K2(1), K2(1)], [K2(2), K2(2)]),
         "is_isometric: square-class invariance over Q(sqrt2)")
    # The blind-spot binary: eps = 35 + 6 sqrt34 is a totally positive
    # nonsquare, so <1, -eps> is anisotropic while <1, -1> is hyperbolic.
    eps34 = K34(35 + 6 * s34)
    gate(!oracle_isometric(K34, [K34(1), -eps34], [K34(1), K34(-1)]),
         "is_isometric: blind-spot binary vs hyperbolic over Q(sqrt34)")

    # --- hasse_invariant convention -------------------------------------------
    # The calibration form <1,1,-3,-3,1,-1> over Q has Hasse invariant -1 at 3
    # in Lam's prod_{i<j} convention (Macaulay2-validated; also reproduced by
    # Magma in the July 2026 calculator runs).
    P3 = prime_decomposition(maximal_order(KQ), 3)[1][1]
    gate(oracle_hasse(KQ, KQ.([1, 1, -3, -3, 1, -1]), P3) == -1,
         "hasse_invariant: calibration form at 3 (Lam prod_{i<j})")
    # (2,3)_2 = -1 over Q_2 by the classical formula
    # (2, b)_2 = (-1)^((b^2-1)/8) with b = 3: (9-1)/8 = 1.
    P2Q = prime_decomposition(maximal_order(KQ), 2)[1][1]
    gate(oracle_hasse(KQ, KQ.([2, 3]), P2Q) == -1,
         "hasse_invariant: (2,3) at the dyadic place of Q")

    # --- is_isotropic, LOCAL --------------------------------------------------
    # Over K7 = Q(sqrt-7) the prime 2 splits, so the completion is Q_2 itself:
    # <1,1,1,1> is the norm form of the division quaternion (-1,-1)_{Q_2},
    # hence anisotropic; rank 5 is always isotropic over a local field.
    P2K7 = prime_decomposition(maximal_order(K7), 2)[1][1]
    gate(!oracle_isotropic_at(K7, K7.([1, 1, 1, 1]), P2K7),
         "is_isotropic local: <1,1,1,1> anisotropic over Q_2")
    gate(oracle_isotropic_at(K7, K7.([1, 1, 1, 1, 1]), P2K7),
         "is_isotropic local: rank 5 isotropic over Q_2")
    # Over the RAMIFIED dyadic completion Q_2(sqrt2) and the INERT one
    # Q_2(sqrt5), every quadratic extension of Q_2 splits the division
    # quaternion, so <1,1,1,1> becomes isotropic.
    P2K2 = prime_decomposition(maximal_order(K2), 2)[1][1]
    gate(oracle_isotropic_at(K2, K2.([1, 1, 1, 1]), P2K2),
         "is_isotropic local: <1,1,1,1> isotropic over ramified Q_2(sqrt2)")
    P2K5 = prime_decomposition(maximal_order(Kr5), 2)[1][1]
    gate(oracle_isotropic_at(Kr5, Kr5.([1, 1, 1, 1]), P2K5),
         "is_isotropic local: <1,1,1,1> isotropic over inert Q_2(sqrt5)")
    # Odd place, binary unit form: <1,-3> at 5 is anisotropic (3 is a
    # nonsquare mod 5); <1,-4> is isotropic (4 = 2^2).
    P5 = prime_decomposition(maximal_order(KQ), 5)[1][1]
    gate(!oracle_isotropic_at(KQ, KQ.([1, -3]), P5),
         "is_isotropic local: <1,-3> anisotropic at 5")
    gate(oracle_isotropic_at(KQ, KQ.([1, -4]), P5),
         "is_isotropic local: <1,-4> isotropic at 5")
    # Real places: <1,1> definite, <1,-s2> indefinite at exactly one of the
    # two real places of Q(sqrt2) (the embeddings send s2 to +-sqrt2).
    rp = real_places(K2)
    gate(!oracle_isotropic_at(K2, K2.([1, 1]), rp[1]) &&
         !oracle_isotropic_at(K2, K2.([1, 1]), rp[2]),
         "is_isotropic local: <1,1> definite at both real places")
    mixed = [oracle_isotropic_at(K2, [K2(1), -s2], v) for v in rp]
    gate(sort(mixed) == [false, true],
         "is_isotropic local: <1,-sqrt2> indefinite at exactly one real place")

    # --- the known GLOBAL is_isotropic defect ---------------------------------
    # Hecke's global is_isotropic reports the definite <1,1,1> over Q(sqrt2)
    # isotropic, which is false. The oracle layer never calls it; this gate
    # exists so a Hecke release that FIXES the defect is noticed immediately.
    hecke_global = is_isotropic(oracle_space(K2, K2.([1, 1, 1])))
    hecke_global &&
        @info "Hecke's global is_isotropic still has the definite-form defect (expected; the oracle layer avoids it)"
    hecke_global ||
        @warn "Hecke's global is_isotropic now answers the definite gate correctly — the defect may be fixed upstream; oracle_isotropic_global could be simplified"

    # --- discriminant convention ----------------------------------------------
    # Hecke's discriminant(V) is the SIGNED discriminant: <1,1> -> -1 (n = 2,
    # sign (-1)^1), <1,1,1> -> -1 (n = 3, sign (-1)^3).
    gate(oracle_signed_disc(KQ, KQ.([1, 1])) == KQ(-1),
         "discriminant: signed convention on <1,1>")
    gate(oracle_signed_disc(KQ, KQ.([1, 1, 1])) == KQ(-1),
         "discriminant: signed convention on <1,1,1>")

    # --- rank-2 exactness of the global oracle ---------------------------------
    # The blind-spot binary <1,-eps> over Q(sqrt34) is locally isotropic at
    # every real and dyadic place yet globally anisotropic; the rank-2 branch
    # of oracle_isotropic_global must catch it via is_square.
    gate(!oracle_isotropic_global(K34, [K34(1), -eps34]),
         "oracle_isotropic_global: blind-spot binary is anisotropic")
    gate(oracle_isotropic_global(K34, [K34(1), K34(-1)]),
         "oracle_isotropic_global: hyperbolic plane is isotropic")

    println("all Hecke convention gates passed")
    nothing
end
