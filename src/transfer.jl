# Port of M2 Code/GWTransfer.m2: transferGW — the transfer GW(L) -> GW(k) along a
# finite étale algebra L/k, induced by the trace Tr_{L/k}: L -> k. The transferred
# Gram matrix is built from the Gram matrix over L itself, not from a
# diagonalization, so nothing has to be inverted in L, which need not be a field.
# (The GWClass constructor has already checked that L is finite étale.) A number
# field over its base field is a finite étale algebra too, and is transferred
# with Hecke's trace and basis.

"""
    transfer_gw(beta::GWClass)

The image of a Grothendieck–Witt class over a finite étale algebra ``L/k``, or
over a number field ``L`` with base field ``k``, under the transfer map
``\\text{GW}(L) → \\text{GW}(k)`` induced by the trace
``\\text{Tr}_{L/k} : L → k``. A form ``β`` on ``L^n`` is sent to the form
``(v, w) \\mapsto \\text{Tr}_{L/k}(β(v, w))`` on ``L^n`` regarded as a
``k``-vector space, so the rank is multiplied by ``[L : k]``. For example,
``⟨a⟩`` is sent to the scaled trace form ``(x, y) \\mapsto \\text{Tr}_{L/k}(axy)``,
of rank ``[L : k]``.

The Gram matrix is written in the ``k``-basis ``e_i b_s`` of ``L^n``, with ``b``
running fastest, where ``e_1, …, e_n`` is the standard basis of ``L^n`` and
``b_1, …, b_d`` are the standard monomials of ``L`` in the order used by
[`multiplication_matrix`](@ref): its ``(i, j)`` block is
``(\\text{Tr}_{L/k}(β_{ij} b_s b_t))_{s,t}``. Its determinant is
``\\det(T)^n N_{L/k}(\\det β)``, where ``T = (\\text{Tr}_{L/k}(b_s b_t))_{s,t}`` is
the trace form of ``L/k``.

For a number field ``L`` the basis ``b`` is Hecke's `basis(L)` over the base
field ``k``, which is ``\\mathbb{Q}`` for an absolute field. For a simple field
``L = k(a)`` it is the power basis ``1, a, …, a^{d-1}``, so the result is the same
matrix as for the algebra ``k[t]/(f)`` presenting ``L``, with ``f`` the defining
polynomial of ``a``. Any number field is accepted, relative or non-simple
included; the other functions of the package reject non-simple fields because
their computations need a simple one, and the transfer does not.

!!! note
    The transferred form is nondegenerate exactly when ``\\det β`` is a unit of
    ``L``, which always holds when ``L`` is a field. When ``L`` is a product of
    fields and ``\\det β`` is a zero divisor, ``β`` is degenerate on a factor of
    ``L``, and `transfer_gw` errors.

# Examples

Over ``L = \\mathbb{Q}[t]/(t^2 - 1) \\cong \\mathbb{Q} \\times \\mathbb{Q}`` the
transfer is the sum of the specializations at ``t = 1`` and ``t = -1``:

```julia-repl
julia> S, (t,) = polynomial_ring(QQ, ["t"]);

julia> A, _ = quo(S, ideal(S, [t^2 - 1]));

julia> beta = GWClass(A[A(1) A(2); A(2) A(t)]);

julia> transfer_gw(beta)
[2   0   4   0]
[0   2   0   4]
[4   0   0   2]
[0   4   2   0]

julia> is_isomorphic_form(transfer_gw(beta), gw_direct_sum(GWClass(QQ[1 2; 2 1]), GWClass(QQ[1 2; 2 -1])))
true
```

Over a number field the transfer goes to the base field. On ``\\mathbb{Q}(\\sqrt 2)``
with basis ``1, \\sqrt 2`` the trace form is ``\\mathrm{diag}(2, 4)``, and
``⟨3 + 5\\sqrt 2⟩`` goes to that form times multiplication by ``3 + 5\\sqrt 2``:

```julia-repl
julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^2 - 2, "a");

julia> transfer_gw(GWClass(matrix(K, 1, 1, [3 + 5a])))
[ 6   20]
[20   12]
```

See also [`algebra_trace`](@ref), [`algebra_norm`](@ref), [`GWClass`](@ref).
"""
function transfer_gw(alpha::GWClass)
    Alg = gw_algebra(alpha)
    Alg isa MPolyQuoRing ||
        error("transferGW requires a class over a finite étale algebra or a number field")
    kk = coefficient_ring(base_ring(Alg))
    M = gw_matrix(alpha)
    n = nrows(M)
    # det M is a unit exactly when multiplication by it is invertible
    iszero(algebra_norm(Alg, det(M))) &&
        error("the determinant of the form is a zero divisor in the algebra, " *
              "so its transfer is degenerate")
    mons = reverse(monomial_basis(Alg))   # the basis of multiplication_matrix
    d = length(mons)
    T = matrix(kk, [algebra_trace(Alg, mons[s] * mons[t]) for s in 1:d, t in 1:d])
    G = zero_matrix(kk, n * d, n * d)
    for j in 1:n, i in 1:j
        # Tr(M[i, j] b_s b_t) is entry (s, t) of T times multiplication by M[i, j]
        TM = T * multiplication_matrix(Alg, M[i, j])
        for s in 1:d, t in 1:d
            G[(i - 1) * d + s, (j - 1) * d + t] = TM[s, t]
            G[(j - 1) * d + t, (i - 1) * d + s] = TM[s, t]
        end
    end
    GWClass(G)
end

# A number field L over its base field k is a finite étale k-algebra as well.
# Hecke computes Tr_{L/k} and supplies the basis b of L over k, which for a
# simple field is the power basis, so L = k(a) gives the matrix of k[t]/(f).
# L is a field, so det(beta) is a unit and the transfer is nondegenerate.
function transfer_gw(alpha::GWClass{<:MatElem{<:NumFieldElem}})
    M = gw_matrix(alpha)
    L = base_ring(M)
    kk = base_field(L)
    n = nrows(M)
    b = basis(L)
    d = length(b)
    G = zero_matrix(kk, n * d, n * d)
    for j in 1:n, i in 1:j, s in 1:d, t in 1:d
        # entry (s, t) of block (i, j) is Tr(M[i, j] b_s b_t)
        entry = tr(M[i, j] * b[s] * b[t])
        G[(i - 1) * d + s, (j - 1) * d + t] = entry
        G[(j - 1) * d + t, (i - 1) * d + s] = entry
    end
    GWClass(G)
end
