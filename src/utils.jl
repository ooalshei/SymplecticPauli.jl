"""
    _check_string_length(string, Q)

Throw unless `string` fits in the `2Q` bits a `Q`-qubit Pauli string occupies.

The low `Q` bits of a symplectic string hold the X part and the high `Q` bits the Z part, so
a valid string occupies exactly `2Q` bits. Testing `string >> 2Q` instead of comparing
against `4^Q - 1` keeps the check to a single shift and stays correct when `2Q` reaches (or
exceeds) the width of the integer type, where `4^Q` would silently overflow.
"""
@inline _check_string_length(string::Unsigned, Q::Integer) =
    if iszero(string >> (2 * Q))
        nothing
    else
        throw(ArgumentError("String must not exceed $(big(4)^Q - 1)."))
    end

"""
    _check_type(T, Q)

Throw unless the unsigned type `T` is wide enough for a `Q`-qubit string, i.e. has at least
`2Q` bits. Types without a `typemax` are unbounded and always pass.
"""
function _check_type(::Type{T}, s::Integer) where {T<:Unsigned}
    Base.hastypemax(T) &&
        8 * sizeof(T) < 2 * s &&
        throw(
            ArgumentError(
                "String length cannot exceed $(4 * sizeof(T)). " *
                "Consider using a larger unsigned integer type.",
            ),
        )
    nothing
end

# The mask below is built in the string's own type: `2^Q - 1` is `Int` arithmetic and
# overflows to `-1` at `Q = 64`, which would mask in the Z half as well.
"""
    countx(p)
    countx(string::Unsigned, Q)

Number of qubits of a Pauli string whose `X` bit is set.

A `Y` sets both bits of its qubit, so it is counted here as well as by [`countz`](@ref) and
[`county`](@ref): `countx` is the number of `X` *and* `Y` factors. [`counti`](@ref) folds the
two halves together instead, and counts a `Y` once.

# Examples
```jldoctest
julia> p = UPauli("XY-Z");

julia> countx(p), county(p), countz(p), counti(p)
(2, 1, 2, 1)

julia> countx(p.string, p.qubits)
2
```
"""
function countx(string::Unsigned, Q::Integer)
    _check_string_length(string, Q)
    return count_ones(string & ((one(string) << Q) - one(string)))
end
countx(p::AbstractPauli) = countx(p.string, p.qubits)

"""
    county(p)
    county(string::Unsigned, Q)

Number of `Y` factors of a Pauli string — the qubits with both bits set.

This is the exponent of the ``i^{\\#Y}`` phase that separates the bare symplectic string
from the Hermitian Pauli operator, so it turns up wherever a phase does: [`Pauli`](@ref),
[`PauliSentence`](@ref) coefficients, [`tostring`](@ref).

# Examples
```jldoctest
julia> county(UPauli("XY-Z")), county(UPauli("YYY-"))
(1, 3)

julia> Pauli(UPauli("YYY-")).sign == im^county(UPauli("YYY-"))
true
```
"""
county(string::Unsigned, Q::Integer) =
    (_check_string_length(string, Q); count_ones(string & (string >> Q)))
county(p::AbstractPauli) = county(p.string, p.qubits)

"""
    countz(p)
    countz(string::Unsigned, Q)

Number of qubits of a Pauli string whose `Z` bit is set — the number of `Z` *and* `Y`
factors, for the same reason as [`countx`](@ref).

# Examples
```jldoctest
julia> countz(UPauli("XY-Z"))
2
```
"""
countz(string::Unsigned, Q::Integer) =
    (_check_string_length(string, Q); count_ones(string >> Q))
countz(p::AbstractPauli) = countz(p.string, p.qubits)

"""
    counti(p)
    counti(string::Unsigned, Q)

Number of identity factors of a Pauli string — the qubits with neither bit set.

`Q` minus this is the *weight* of the string, the number of qubits it acts on, so sorting by
`counti` orders strings from the identity outwards. Unlike the other three counts it is not a
popcount of one half: a `Y` is one non-identity qubit even though it sets two bits, so the
two halves are folded together before counting.

# Examples
```jldoctest
julia> counti(UPauli("XZ--"))
2

julia> counti(UPauli("XY-Z"))       # a Y counts once, like any other factor
1

julia> counti(UPauli("----")), counti(UPauli("YYYY"))
(4, 0)
```
"""
function counti(string::Unsigned, Q::Integer)
    _check_string_length(string, Q)
    mask = (one(string) << Q) - one(string)
    # A qubit is an identity when neither of its bits is set, so OR the two halves together
    # and count what is left: a `Y` is one non-identity qubit, not two or three.
    return Q - count_ones((string | (string >> Q)) & mask)
end
counti(p::AbstractPauli) = counti(p.string, p.qubits)

"""
    tostring(p::AbstractPauli)
    tostring(v::PauliList)
    tostring(s::PauliSentence)

Human-readable form of a Pauli string, list or sentence: `X`, `Y`, `Z` and `-`, leftmost
qubit first.

For a [`Pauli`](@ref) the phase is prefixed, e.g. `"(i)Z-"`. For a [`PauliSentence`](@ref)
the result is a `Dict` from string to coefficient, with the ``i^{\\#Y}`` convention phase
divided out so the coefficients read as they were entered.

# Examples
```jldoctest
julia> tostring(UPauli("XY-Z"))
"XY-Z"

julia> tostring(PauliList(["ZZ", "-X"]))
2-element Vector{String}:
 "ZZ"
 "-X"

julia> tostring(PauliSentence(PauliList(["ZZ"]), [0.5]))
Dict{String, ComplexF64} with 1 entry:
  "ZZ" => 0.5+0.0im
```
"""
function tostring(p::UPauli)::String
    result = ""
    string = digits(p.string, base=2, pad=2 * p.qubits)
    for i in 1:p.qubits
        if string[i] == string[i+p.qubits] == 1
            result *= "Y"
        elseif string[i] == 1
            result *= "X"
        elseif string[i+p.qubits] == 1
            result *= "Z"
        else
            result *= "-"
        end
    end
    return result
end
function tostring(p::Pauli)::String
    sign = p.sign * C8(-im)^county(p)
    sign == 1 && return "(+)" * tostring(UPauli(p))
    sign == -1 && return "(-)" * tostring(UPauli(p))
    sign == im && return "(i)" * tostring(UPauli(p))
    sign == -im && return "(-i)" * tostring(UPauli(p))
end
function tostring(s::PauliSentence)
    result = Dict{String,valtype(s)}()
    for (key, value) in pairs(s)
        result[tostring(UPauli(UInt(key), s.qubits))] = (-im)^county(key, s.qubits) * value
    end
    return result
end
tostring(v::PauliList) = tostring.(UPauli.(v.strings, v.qubits))

@doc raw"""
    tomatrix(p::AbstractPauli)
    tomatrix(s::PauliSentence)
    tomatrix(string::Unsigned, Q)

Dense ``2^Q \times 2^Q`` matrix of a Pauli string or sentence.

`tomatrix(string, Q)` gives the *bare* embedding — `X` ↦ ``σ_1``, `Z` ↦ ``σ_3``, a `Y`
position ↦ the real antisymmetric ``σ_2`` — which is what [`PauliSentence`](@ref)
coefficients multiply. `tomatrix` of a [`Pauli`](@ref) or [`UPauli`](@ref) includes the
phase that makes it the Hermitian Pauli operator.

Exponential in `Q`: for checking small cases, not for computing with.

# Examples
```jldoctest
julia> tomatrix(Pauli("X"))
2×2 Matrix{Complex{Int8}}:
 0+0im  1+0im
 1+0im  0+0im

julia> tomatrix(PauliSentence(PauliList(["Z"]), [2.0]))
2×2 Matrix{ComplexF64}:
 2.0+0.0im   0.0+0.0im
 0.0+0.0im  -2.0+0.0im
```
"""
function tomatrix(string::Unsigned, Q::Integer)
    _check_string_length(string, Q)
    result = I(1)
    string = digits(string, base=2, pad=2 * Q)
    for i in 1:Q
        if string[i] == string[i+Q] == 1
            result = result ⊗ σ₂real  # σy
        elseif string[i] == 1
            result = result ⊗ σ₁  # σx
        elseif string[i+Q] == 1
            result = result ⊗ σ₃  # σz
        else
            result = result ⊗ I(2)  # I
        end
    end
    return result
end
tomatrix(p::AbstractPauli) = Pauli(p).sign * tomatrix(p.string, p.qubits)

function tomatrix(p::PauliSentence)
    result = zeros(ComplexF64, 2^p.qubits, 2^p.qubits)
    for (key, value) in p
        result .+= value .* tomatrix(key, p.qubits)
    end
    return result
end
