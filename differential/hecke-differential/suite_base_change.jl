# gw_base_change GW(Q) -> GW(K): rank preserved, every K-signature equals
# the Q-signature, the discriminant maps to the image square class, sums and
# tensors commute with the change of base, and isotropy is monotone.

# Rational diagonals: large square content, denominators, high rank,
# definite and indefinite blocks.
const BASE_CHANGE_QQ_DIAGONALS = [
    [QQ(1), QQ(-1)],
    [QQ(1), QQ(1), QQ(1)],
    [QQ(2), QQ(-3), QQ(5)],
    [QQ(1024), QQ(-9), QQ(49)],
    [QQ(1 // 2), QQ(-3 // 4), QQ(2 // 9)],
    [QQ(7), QQ(7), QQ(-7), QQ(-7)],
    [QQ(1), QQ(2), QQ(-2), QQ(3), QQ(-6), QQ(30)],
    [QQ(-15)],
    [QQ(1000003), QQ(-999983)],
    [QQ(2), QQ(2), QQ(2), QQ(2)],
    [QQ(1), QQ(-1), QQ(1), QQ(-1), QQ(5), QQ(-5), QQ(7), QQ(-30)],
    [QQ(-1 // 1024), QQ(343)],
]

function suite_base_change!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        for Dq in BASE_CHANGE_QQ_DIAGONALS
            q = diagonal_form(QQ, Tuple(Dq))
            bc = gw_base_change(K, q)
            DK = [K(a) for a in Dq]

            check!(res, "gw_base_change", gw_matrix(bc) == oracle_gram(K, DK),
                   "entry-wise image: field=$(f.name) diagQ=$(Dq)")
            check!(res, "gw_base_change", form_rank(bc) == length(Dq),
                   "rank preserved: field=$(f.name) diagQ=$(Dq)")

            sigQ = form_signature(q)
            for sigma in real_places(K)
                check!(res, "gw_base_change", form_signature(bc, sigma) == sigQ,
                       "signature at every real place = Q-signature: field=$(f.name) diagQ=$(Dq)")
            end

            dQ = signed_discriminant(q)
            dK = signed_discriminant(bc)
            check!(res, "gw_base_change", is_square(K(dQ) * K(dK)),
                   "discriminant maps to the image class: field=$(f.name) diagQ=$(Dq)")

            # Isotropy is monotone under field extension.
            if is_isotropic_form(q)
                check!(res, "gw_base_change", oracle_isotropic_global(K, DK),
                       "isotropy survives base change (oracle): field=$(f.name) diagQ=$(Dq)")
            end
        end

        # Base change is a semiring map: it commutes with ⟂ and ⊗ literally.
        for _ in 1:3
            D1 = BASE_CHANGE_QQ_DIAGONALS[det_int!(rng, length(BASE_CHANGE_QQ_DIAGONALS))]
            D2 = BASE_CHANGE_QQ_DIAGONALS[det_int!(rng, length(BASE_CHANGE_QQ_DIAGONALS))]
            q1, q2 = diagonal_form(QQ, Tuple(D1)), diagonal_form(QQ, Tuple(D2))
            check!(res, "gw_base_change",
                   gw_matrix(gw_base_change(K, gw_direct_sum(q1, q2))) ==
                       gw_matrix(gw_direct_sum(gw_base_change(K, q1), gw_base_change(K, q2))),
                   "commutes with direct sum: field=$(f.name)")
            if length(D1) * length(D2) <= 12
                check!(res, "gw_base_change",
                       gw_matrix(gw_base_change(K, gw_tensor_product(q1, q2))) ==
                           gw_matrix(gw_tensor_product(gw_base_change(K, q1),
                                                       gw_base_change(K, q2))),
                       "commutes with tensor product: field=$(f.name)")
            end
        end
    end
end
