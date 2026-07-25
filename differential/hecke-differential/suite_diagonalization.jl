# diagonalize_via_congruence, diagonal_entries, and GWClass construction on
# dense symmetric inputs (T' D T): Hecke isometry between input and output
# plus square-class preservation of the determinant.

function suite_diagonalization!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        for A in dense_pool(f, rng)
            n = nrows(A)
            iszero(det(A)) && continue

            D = diagonalize_via_congruence(A)
            offdiag_zero = all(iszero(D[i, j]) for i in 1:n for j in 1:n if i != j)
            check!(res, "diagonalize_via_congruence", offdiag_zero,
                   "output is diagonal: field=$(f.name) A=$(A)")
            check!(res, "diagonalize_via_congruence",
                   oracle_isometric_gram(A, D),
                   "congruent to the input (Hecke isometry): field=$(f.name) A=$(A)")
            check!(res, "diagonalize_via_congruence",
                   is_square(det(A) * det(D)),
                   "determinant square class preserved: field=$(f.name) A=$(A)")

            # diagonal_entries on the dense class: entries multiply to the
            # determinant's square class, and rebuilding a diagonal form from
            # them recovers the class (Hecke-checked).
            beta = GWClass(A)
            ents = [K(a) for a in diagonal_entries(beta)]
            check!(res, "diagonal_entries",
                   length(ents) == n && is_square(det(A) * prod(ents; init = one(K))),
                   "entries multiply to det's class: field=$(f.name) A=$(A)")
            check!(res, "diagonal_entries", oracle_isometric_gram(A, oracle_gram(K, ents)),
                   "entries rebuild the class: field=$(f.name) A=$(A)")

            # GWClass construction: accepts the symmetric nondegenerate A,
            # rejects a non-symmetric perturbation and a visibly singular one.
            check!(res, "GWClass", gw_matrix(beta) == A,
                   "constructor stores the Gram matrix: field=$(f.name)")
            asym = deepcopy(A)
            asym[1, min(2, n)] += 1     # n >= 2 for every dense-pool matrix
            check!(res, "GWClass",
                   !_constructs(asym),
                   "rejects non-symmetric input: field=$(f.name)")
            singular = zero_matrix(K, n, n)
            check!(res, "GWClass",
                   !_constructs(singular),
                   "rejects the zero matrix: field=$(f.name)")
        end
    end
end

"Does GWClass accept the matrix without throwing?"
function _constructs(M::MatElem)::Bool
    try
        GWClass(M)
        true
    catch
        false
    end
end
