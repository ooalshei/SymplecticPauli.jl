@doc raw"""
    AbstractPauli{T<:Unsigned,Q}

A single Pauli string on `Q` qubits, stored in `T` as `2Q` symplectic bits.

The low `Q` bits mark `X` support and the high `Q` bits mark `Z` support, so a `Y` sets both
bits of its qubit. Multiplying two strings is `⊻` on the bits; the accompanying phase is
what separates the two concrete types, [`UPauli`](@ref) and [`Pauli`](@ref).

Both parameters are part of the type: `T` is the unsigned integer the bits live in — `UInt8`
holds 4 qubits, `UInt` 32, `UInt128` 64 — and `Q` is the number of qubits, so strings of
different widths cannot be multiplied by accident. Every subtype is a bits type, free to
build and to store in an array unboxed.

Subtypes carry a `string` field (the bits) and a `qubits` field (`Q` again, as an `Int`).

# Examples
```jldoctest
julia> p = UPauli{UInt16}("XY-Z")
Pauli{UInt16}(XY-Z)

julia> p isa AbstractPauli{UInt16,4}
true

julia> p.string, p.qubits
(0x00a3, 4)
```
"""
abstract type AbstractPauli{T<:Unsigned,Q} end
Base.copy(p::AbstractPauli) = p
Base.show(io::IO, p::AbstractPauli{T}) where {T} = print(io, "Pauli{$T}(", tostring(p), ")")
# toint(p::AbstractPauli) = p.string
# toint(::Type{T}, p::AbstractPauli) where {T<:Integer} = T(p.string)

@doc raw"""
    UPauli(string::Unsigned, Q)
    UPauli(text::AbstractString)
    UPauli(indices::AbstractVector{<:Integer})
    UPauli(p::AbstractPauli)

An *unsigned* Pauli string: the bit pattern only, with no phase attached.

The text form uses `X`, `Y`, `Z` and `-` (or `I`) for the identity, leftmost character
first; the index form uses `1`–`4` for `I`, `X`, `Y`, `Z`. A `UPauli` is a bits type, so it
is free to build and pass around.

Passing an [`AbstractPauli`](@ref) keeps the bits and drops the phase. Prefixing any
constructor with an integer type — `UPauli{UInt16}("XY")` — chooses the storage; without one
it is `UInt`, which covers 32 qubits. Building from raw bits also takes the qubit count,
since the bits alone do not say how many trailing identities there are.

The product of two `UPauli`s is a [`Pauli`](@ref), because it generally carries a phase; the
phaseless product — the bit pattern alone — is `p ⊻ q`.

# Examples
```jldoctest
julia> p = UPauli("XY-Z")
Pauli{UInt64}(XY-Z)

julia> p.string, p.qubits
(0x00000000000000a3, 4)

julia> tostring(UPauli([2, 3, 1, 4]))
"XY-Z"

julia> UPauli(UInt(0b1001), 2)       # bits given directly: X on qubit 1, Z on qubit 2
Pauli{UInt64}(XZ)

julia> tostring(p ⊻ UPauli("Z-ZZ"))  # product without its phase
"YYZ-"
```
"""
struct UPauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    qubits::Int
    UPauli{T,Q}(string) where {T,Q} = (_check_string_length(string, Q); new{T,Q}(string, Q))
end
UPauli(string::T, Q::Integer) where {T<:Unsigned} = UPauli{T,Q}(string)
# UPauli{T,Q}(p::AbstractPauli) where {T,Q} = UPauli{T,Q}(p.string)
UPauli{T}(p::AbstractPauli) where {T} = UPauli{T,p.qubits}(p.string)
UPauli(p::AbstractPauli{T,Q}) where {T,Q} = UPauli{T,Q}(p.string)
function UPauli{T}(p::AbstractString) where {T}
    unique(p) ⊆ ['X', 'Y', 'Z', 'I', '-'] ||
        throw(ArgumentError("String must contain only 'X', 'Y', 'Z', 'I', or '-'."))
    Q = length(p)
    _check_type(T, Q)
    number = T(0)
    # `one(T) << Q`, not `1 << Q`: an `Int` literal would promote the accumulator out of `T`
    # and overflow the shift for every type narrower than `Int`.
    for char in Iterators.reverse(p)
        number <<= 1
        if char == 'X'
            number |= one(T)
        elseif char == 'Z'
            number |= (one(T) << Q)
        elseif char == 'Y'
            number |= (one(T) << Q)
            number |= one(T)
        end
    end
    return UPauli{T,Q}(number)
end
function UPauli{T}(p::AbstractVector{<:Integer}) where {T}
    unique(p) ⊆ [1, 2, 3, 4] ||
        throw(ArgumentError("Array must contain only 1, 2, 3, or 4."))
    Q = length(p)
    _check_type(T, Q)
    number = T(0)
    for ind in Iterators.reverse(p)
        number <<= 1
        if ind == 2
            number |= one(T)
        elseif ind == 4
            number |= (one(T) << Q)
        elseif ind == 3
            number |= (one(T) << Q)
            number |= one(T)
        end
    end
    UPauli{T,Q}(number)
end
UPauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = UPauli{UInt}(p)
# Base.:(==)(p::UPauli, q::UPauli) = p.string == q.string

@doc raw"""
    Pauli(string::Unsigned, sign, Q)
    Pauli(p::UPauli, [sign])
    Pauli(text::AbstractString, [sign])
    Pauli(indices::AbstractVector{<:Integer}, [sign])

A Pauli string carrying a phase, one of ``1, -1, i, -i``.

`Pauli(p::UPauli)` attaches the phase ``i^{\#Y}`` that turns the bare bit pattern into the
Hermitian Pauli operator, which is the convention [`tomatrix`](@ref) and
[`PauliSentence`](@ref) coefficients follow. Products of `Pauli`s track their phase exactly.
Give a `sign` explicitly to override that convention; anything other than a fourth root of
unity is an `ArgumentError`.

The phase is stored as `Complex{Int8}`, so it is exact and multiplying by it never rounds.

# Examples
```jldoctest
julia> Pauli("Y").sign
0 + 1im

julia> tomatrix(Pauli("Y"))
2×2 Matrix{Complex{Int8}}:
 0+0im  0-1im
 0+1im  0+0im

julia> UPauli("X") * UPauli("Y")        # XY = iZ
Pauli{UInt64}((i)Z)

julia> Pauli("Z-", -im) * Pauli("ZZ")   # phases multiply through
Pauli{UInt64}((-i)-Z)

julia> Pauli(UInt(1), 2, 1)
ERROR: ArgumentError: Sign must be 1, -1, im, or -im
```
"""
struct Pauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    sign::C8
    qubits::Int
    function Pauli{T,Q}(string::Integer, sign::Number) where {T,Q}
        _check_string_length(string, Q)
        if (sign == 1) | (sign == -1) | (sign == im) | (sign == -im)
            new{T,Q}(string, sign, Q)
        else
            throw(ArgumentError("Sign must be 1, -1, im, or -im"))
        end
    end
end
Pauli(string::T, sign::Number, Q::Integer) where {T<:Unsigned} = Pauli{T,Q}(string, sign)
Pauli{T,Q}(string::Unsigned) where {T,Q} = Pauli{T,Q}(string, 1)
Pauli(string::T, Q::Integer) where {T<:Unsigned} = Pauli{T,Q}(string, 1)
Pauli{T}(p::UPauli) where {T} = Pauli{T,p.qubits}(p.string, (im)^county(p))
Pauli(p::UPauli{T,Q}) where {T,Q} = Pauli{T}(p)
Pauli{T}(p::UPauli, sign::Number) where {T} = Pauli{T,p.qubits}(p.string, sign)
Pauli(p::UPauli{T,Q}, sign::Number) where {T,Q} = Pauli{T,Q}(p.string, sign)
Pauli{T}(p::Pauli) where {T} = Pauli{T,p.qubits}(p.string, p.sign)
Pauli(p::Pauli) = p
Pauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}) where {T} =
    Pauli(UPauli{T}(p), (im)^county(UPauli(p)))
Pauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = Pauli{UInt}(p)
Pauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) where {T} =
    Pauli(UPauli{T}(p), sign)
Pauli(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) =
    Pauli(UPauli(p), sign)
# Base.:(==)(p1::Pauli, p2::Pauli) = p1.string == p2.string & p1.sign == p2.sign
