# The low `Q` bits of a symplectic string hold the X part and the high `Q` bits the Z part,
# so a valid string occupies exactly `2Q` bits. Testing `string >> 2Q` instead of comparing
# against `4^Q - 1` keeps the check to a single shift and stays correct when `2Q` reaches
# (or exceeds) the width of the integer type, where `4^Q` would silently overflow.
@inline _check_string_length(string::Unsigned, Q::Integer) = if iszero(string >> (2 * Q))
    nothing
else
    throw(ArgumentError("String must not exceed $(big(4)^Q - 1)."))
end

function _check_type(::Type{T}, s::Integer) where {T<:Unsigned}
    Base.hastypemax(T) &&
        8 * sizeof(T) < 2 * s &&
        throw(ArgumentError("String length cannot exceed $(4 * sizeof(T)). Consider using
                            a larger unsigned integer type."))
    nothing
end

"""
    countx(p), county(p), countz(p), counti(p)

Bit counts of a Pauli string: how many qubits carry the `X` bit, the `Y` (both bits), the
`Z` bit, and `Q` minus the three.

A `Y` sets both bits of its qubit, so it is counted by `countx` *and* `countz` as well as by
`county`, and `counti` therefore subtracts it three times. `counti` is used as a weight
ordering — larger means closer to the identity — rather than as a literal count of identity
factors.

# Examples
```jldoctest
julia> p = UPauli("XY-Z");

julia> countx(p), county(p), countz(p), counti(p)
(2, 1, 2, -1)

julia> counti(UPauli("XZ--"))       # no Y: the plain identity count
2
```
"""
countx(string::Unsigned, Q::Integer) =
    (_check_string_length(string, Q); count_ones(string & (2^Q - 1)))
countx(p::AbstractPauli) = countx(p.string, p.qubits)

county(string::Unsigned, Q::Integer) =
    (_check_string_length(string, Q); count_ones(string & (string >> Q)))
county(p::AbstractPauli) = county(p.string, p.qubits)

countz(string::Unsigned, Q::Integer) =
    (_check_string_length(string, Q); count_ones(string >> Q))
countz(p::AbstractPauli) = countz(p.string, p.qubits)

counti(string::Unsigned, Q::Integer) =
    Q - countx(string, Q) - county(string, Q) - countz(string, Q)
counti(p::AbstractPauli) = counti(p.string, p.qubits)

"""
    tostring(x)

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
    tomatrix(x)

Dense ``2^Q \times 2^Q`` matrix of a Pauli string or sentence.

`tomatrix(string, Q)` gives the *bare* embedding — `X` ↦ ``σ_1``, `Z` ↦ ``σ_3``, a `Y`
position ↦ the real antisymmetric ``σ_2`` — which is what [`PauliSentence`](@ref)
coefficients multiply. `tomatrix` of a [`Pauli`](@ref) or [`UPauli`](@ref) includes the phase
that makes it the Hermitian Pauli operator.

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
