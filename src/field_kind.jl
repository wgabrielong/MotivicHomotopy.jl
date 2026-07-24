# Field-kind trait: the per-field-species dispatch for the form layer.
#
# The per-field branches used to be `isa QQField` / `isa FinField` chains at
# every entry point; the trait names the field species once and lets each
# entry point split into per-kind methods. Formally-real is a *predicate*,
# not a kind: the number-field algorithms of [KR21] branch on it inside their
# NumberFieldKind methods (real vs non-real changes steps within one
# algorithm, not which function runs).
#
# The ℝ/ℂ float backends (Matrix{Float64} / Matrix{ComplexF64}) never route
# through the trait — they are dispatch-separated by matrix type — so there
# are deliberately no kinds and no `field_kind` methods for them.
#
# [KR21] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic
# form over a number field*, arXiv:2109.04172.

abstract type FieldKind end
struct RationalKind    <: FieldKind end   # QQ
struct FiniteKind      <: FieldKind end   # 𝔽_q, q odd
struct NumberFieldKind <: FieldKind end   # absolute number field K
struct UnsupportedKind <: FieldKind end

field_kind(::QQField) = RationalKind()
field_kind(R::FinField) = characteristic(R) == 2 ? UnsupportedKind() : FiniteKind()
field_kind(::AbsSimpleNumField) = NumberFieldKind()
field_kind(::Any) = UnsupportedKind()

# Umbrella guard for the entry points that support exactly QQ and odd 𝔽_q
# today. Number fields deliberately do NOT pass this guard yet; each entry
# point gains a NumberFieldKind method when its mathematics is ported.
_exact_supported(R) = field_kind(R) isa Union{RationalKind, FiniteKind}

# Only absolute simple number fields are supported; a user who builds
# Q(sqrt2, sqrt3) via number_field([x^2-2, x^2-3]) holds a non-simple field
# and needs the conversion named, not a generic unsupported-field error.
const _NONSIMPLE_FIELD_ERR = "non-simple number fields (AbsNonSimpleNumField) are " *
    "not supported; convert to a simple extension first with simple_extension(K)"

_check_nonsimple(R) =
    R isa AbsNonSimpleNumField ? error(_NONSIMPLE_FIELD_ERR) : nothing

function _require_exact_supported(R, err::AbstractString)
    _exact_supported(R) && return nothing
    _check_nonsimple(R)
    error(err)
end

# Guard for the entry points whose mathematics now covers number fields too
# (invariants, anisotropic dimension, isometry).
_form_supported(R) = field_kind(R) isa Union{RationalKind, FiniteKind, NumberFieldKind}

function _require_form_supported(R, err::AbstractString)
    _form_supported(R) && return nothing
    _check_nonsimple(R)
    error(err)
end

# --- formally real / totally real ---------------------------------------------------
# A field is formally real iff -1 is not a sum of squares; for a number field
# that is equivalent to having at least one real embedding (r₁ > 0), and
# totally real means all embeddings are real (r₂ = 0). Kept internal
# (underscore) so a future Hecke/Oscar export of the same name can never
# create a `using Oscar` ambiguity.

_is_formally_real(::QQField) = true
_is_formally_real(::FinField) = false
_is_formally_real(K::AbsSimpleNumField) = signature(K)[1] > 0

_is_totally_real(::QQField) = true
_is_totally_real(::FinField) = false
_is_totally_real(K::AbsSimpleNumField) = signature(K)[2] == 0
