# transfer_gw over number fields: a class over a number field L goes to its
# transfer along Tr_{L/k}, a class over the base field k, which transfer_gw
# computes with Hecke's trace on Hecke's basis of L. Each value is checked
# against routes that share no code with that: for an absolute field, Hecke's
# restriction of scalars of the quadratic space; for every field, the étale
# algebra k[t]/(f) presenting L, whose transfer is built from multiplication
# matrices. Over each field of the matrix, its own transfer to QQ and that of a
# relative quadratic extension K(sqrt m) to K are checked, and one non-simple
# field closes the suite.

# L = k(a) as the étale algebra k[t]/(f), with f the defining polynomial of a
# and a sent to t, so an element goes to the polynomial in t whose coefficients
# are its coordinates in the power basis.
function _transfer_presentation(L)
    k = base_field(L)
    S, (t,) = polynomial_ring(k, [:t])
    A, _ = quo(S, ideal(S, [defining_polynomial(L)(t)]))
    into_A(e) = (coords = coordinates(e);
                 A(sum(coords[i] * t^(i - 1) for i in eachindex(coords))))
    return into_A
end

# A symmetric n x n matrix over L with nonzero determinant: each entry is a
# combination of `spanning` with coefficients in -3..3 drawn from rng.
function _transfer_gram(L, spanning::Vector, n::Int, rng::DetRNG)
    while true
        B = zero_matrix(L, n, n)
        for i in 1:n, j in i:n
            B[i, j] = B[j, i] = sum((det_int!(rng, 7) - 4) * s for s in spanning)
        end
        iszero(det(B)) || return B
    end
end

_hecke_restriction(B) =
    map_entries(QQ, gram_matrix(restrict_scalars(quadratic_space(base_ring(B), B), QQ)[1]))

function suite_transfer!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K

        # K itself, over QQ: ranks 1 to 3, two forms each
        into_A = _transfer_presentation(K)
        for n in (1, 1, 2, 2, 3, 3)
            B = _transfer_gram(K, [one(K), f.gen], n, rng)
            G = gw_matrix(transfer_gw(GWClass(B)))
            check!(res, "transfer_gw", G == gw_matrix(transfer_gw(GWClass(map(into_A, B)))),
                   "$(f.name), rank $n: the étale algebra presenting K, B = $B")
            check!(res, "transfer_gw", G == _hecke_restriction(B),
                   "$(f.name), rank $n: Hecke's restriction of scalars, B = $B")
        end

        # a relative quadratic extension L = K(sqrt m) over K, transferring to K
        m = first(m for m in (2, 3, 5, 6, 7, -1, -3) if !is_square(K(m)))
        Kz, z = K["z"]
        L, r = number_field(z^2 - m, "r")
        into_AL = _transfer_presentation(L)
        for n in (1, 2, 2)
            B = _transfer_gram(L, [one(L), r, L(f.gen)], n, rng)
            G = gw_matrix(transfer_gw(GWClass(B)))
            check!(res, "transfer_gw",
                   base_ring(G) === K &&
                   G == gw_matrix(transfer_gw(GWClass(map(into_AL, B)))),
                   "$(f.name)(sqrt $m) over $(f.name), rank $n: the étale algebra " *
                   "presenting L, B = $B")
        end
    end

    # A non-simple field: QQ(sqrt2, sqrt3) on Hecke's basis 1, u, v, uv, against
    # the étale algebra QQ[u, v]/(u^2 - 2, v^2 - 3) on its standard monomials
    # 1, u, v, uv, an element going to the same coordinates.
    Qx, x = QQ["x"]
    N, (u, v) = number_field([x^2 - 2, x^2 - 3], ["u", "v"])
    S, (U, V) = polynomial_ring(QQ, ["u", "v"])
    biquadratic, _ = quo(S, ideal(S, [U^2 - 2, V^2 - 3]))
    monomials = reverse(monomial_basis(biquadratic))
    into_biquadratic(e) = (coords = coordinates(e);
                           biquadratic(sum(coords[i] * monomials[i] for i in 1:4)))
    for n in (1, 2, 2)
        B = _transfer_gram(N, [one(N), u, v], n, rng)
        check!(res, "transfer_gw",
               gw_matrix(transfer_gw(GWClass(B))) ==
               gw_matrix(transfer_gw(GWClass(map(into_biquadratic, B)))),
               "QQ(sqrt2, sqrt3) as a non-simple field, rank $n, B = $B")
    end
end
