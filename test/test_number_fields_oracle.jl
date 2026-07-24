# Specialisation oracle (memo §3): forms over QQ base-changed to
# Hecke.rationals_as_number_field() must reproduce every invariant of the
# QQ implementation. Any divergence is a bug in one of the two paths.

using Test
using Oscar
using MotivicHomotopy

const MH = MotivicHomotopy

# The shared number-field test matrix (see nf_test_fields.jl); KQ comes from it.
include("nf_test_fields.jl")
OKQ = maximal_order(KQ)
toKQ(M) = map_entries(x -> KQ(x), M)
P_over(p) = prime_decomposition(OKQ, p)[1][1]

@testset "specialisation: QQ invariants vs rationals_as_number_field" begin
    mats = [
        QQ[2 1; 1 3],
        QQ[1 4 7; 4 3 -1; 7 -1 5],
        QQ[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5],
        gw_matrix(diagonal_form(QQ, (3, -3, 2, 5, 1, -9))),
        gw_matrix(diagonal_form(QQ, (1, -1))),
        gw_matrix(diagonal_form(QQ, (6, 7, 22, -5, 11))),
        gw_matrix(diagonal_form(QQ, (2, 3))),
        gw_matrix(diagonal_form(QQ, (1, 1, 1))),
        gw_matrix(diagonal_form(QQ, (-1, -1, -1))),
        gw_matrix(diagonal_form(QQ, 7)),
        gw_matrix(hyperbolic_form(QQ, 4)),
        gw_matrix(pfister_form(QQ, (2, 6))),
    ]
    classesQ = [GWClass(M) for M in mats]
    classesK = [GWClass(toKQ(M)) for M in mats]
    sigmaKQ = real_places(KQ)[1]

    for (bq, bk) in zip(classesQ, classesK)
        @test form_rank(bq) == form_rank(bk)
        @test form_signature(bq) == form_signature(bk, sigmaKQ)
        @test form_signatures(bk) == [form_signature(bq)]
        # discriminant as a square class
        @test is_square(KQ(QQ(signed_discriminant(bq))) * signed_discriminant(bk))
        # 𝔓(q) matches as rational primes (same deterministic diagonalisation,
        # squarefree part over QQ vs odd-valuation primes over K)
        @test sort([minimum(P) for P in relevant_primes(bk)]) ==
              sort(relevant_primes(bq))
        # Hasse at every relevant prime and at 2 (the unique prime above p)
        for p in unique(vcat(relevant_primes(bq), ZZRingElem[2]))
            @test hasse_witt_invariant(bq, p) == hasse_witt_invariant(bk, P_over(p))
        end
        @test anisotropic_dimension(bq) == anisotropic_dimension(bk)
        @test gw_witt_index(bq) == gw_witt_index(bk)
        @test is_isotropic_form(bq) == is_isotropic_form(bk)
        @test is_anisotropic_form(bq) == is_anisotropic_form(bk)
    end

    # pairwise isometry decisions agree
    for i in eachindex(classesQ), j in i:length(classesQ)
        @test is_isomorphic_form(classesQ[i], classesQ[j]) ==
              is_isomorphic_form(classesK[i], classesK[j])
    end
end

@testset "specialisation: anisotropic_part vs the QQ path" begin
    # The number-field driver, run through rationals_as_number_field, must be
    # isometric to what the QQ path (validated against Macaulay2) returns.
    # Compared by mapping back through the degree-1 isomorphism and calling
    # the QQ predicate — is_isomorphic_form errors on mixed base fields by
    # design and must stay that way. Divergences are reported individually
    # (form, both answers, invariants) and never reconciled.
    #
    # The corpus deliberately includes every negative-signature dimension-3
    # shape from the decomposition suite — the reach of the published M2
    # alpha_1 defect. The QQ port repaired that defect by verifying BOTH
    # candidate signs of its CRT representative (decomposition.jl,
    # _reduce_aniso_dim3), so no divergence is expected; any that appears
    # belongs in the erratum.
    toQQmat(M) = map_entries(x -> coeff(x, 0), M)
    corpus = Any[
        (3, -3, 2, 5, 1, -9), (1, -1), (6, 7, 22, -5, 11), (2, 3), (1, 1, 1),
        (-1, -1, -1), (1, 1, 1, 1, 1, -1), (1, 1, -3, -3, 1, -1),
        (1, 1, 1, 1, -1), (3, 5, 1, -1), (2, 1, -1), (-1, -1, -1, -1, -1, 1),
        (-1, -1, -1, -1, 1), (-1, -2, -5), (-3, -5, -30), (7,), (-2, -3),
        (1, -1, 1, -1),
    ]
    grams = [gw_matrix(diagonal_form(QQ, e)) for e in corpus]
    push!(grams, QQ[2 1; 1 3], QQ[1 4 7; 4 3 -1; 7 -1 5],
          QQ[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5])
    for M in grams
        bq = GWClass(M)
        bk = GWClass(toKQ(M))
        qaQ = anisotropic_part(bq)
        qaK = anisotropic_part(bk)
        agree = form_rank(qaQ) == form_rank(qaK) &&
                (form_rank(qaQ) == 0 ||
                 is_isomorphic_form(qaQ, GWClass(toQQmat(gw_matrix(qaK)))))
        if !agree
            @info "anisotropic_part DIVERGENCE (report, do not reconcile)" gram = M qaQ_entries =
                diagonal_entries(qaQ) qaK_entries = diagonal_entries(qaK) rankQ =
                form_rank(qaQ) rankK = form_rank(qaK) sigQ = form_signature(bq) discQ =
                integral_discriminant(bq)
        end
        @test agree
    end
end
