abstract type AbstractPauli{T<:Unsigned,Q} end
Base.copy(p::AbstractPauli) = p
Base.show(io::IO, p::AbstractPauli{T}) where {T} = print(io, "Pauli{$T}(", tostring(p), ")")
# toint(p::AbstractPauli) = p.string
# toint(::Type{T}, p::AbstractPauli) where {T<:Integer} = T(p.string)

struct UPauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    qubits::Integer
    UPauli{T,Q}(string) where {T,Q} = (_check_string_length(string, Q); new{T,Q}(string, Q))
end
UPauli(string::T, Q::Integer) where {T<:Unsigned} = UPauli{T,Q}(string)
# UPauli{T,Q}(p::AbstractPauli) where {T,Q} = UPauli{T,Q}(p.string)
UPauli{T}(p::AbstractPauli) where {T} = UPauli{T,p.qubits}(p.string)
UPauli(p::AbstractPauli{T,Q}) where {T,Q} = UPauli{T,Q}(p.string)
function UPauli{T}(p::AbstractString) where {T}
    unique(p) ⊆ ['X', 'Y', 'Z', 'I', '-'] || throw(ArgumentError("String must contain only 'X', 'Y', 'Z', 'I', or '-'."))
    Q = length(p)
    _check_type(T, Q)
    number = T(0)
    for char in Iterators.reverse(p)
        number <<= 1
        if char == 'X'
            number |= 1
        elseif char == 'Z'
            number |= (1 << Q)
        elseif char == 'Y'
            number |= (1 << Q)
            number |= 1
        end
    end
    return UPauli{T,Q}(number)
end
function UPauli{T}(p::AbstractVector{<:Integer}) where {T}
    unique(p) ⊆ [1, 2, 3, 4] || throw(ArgumentError("Array must contain only 1, 2, 3, or 4."))
    Q = length(p)
    _check_type(T, Q)
    number = T(0)
    for ind in Iterators.reverse(p)
        number <<= 1
        if ind == 2
            number |= 1
        elseif ind == 4
            number |= (1 << Q)
        elseif ind == 3
            number |= (1 << Q)
            number |= 1
        end
    end
    UPauli{T,Q}(number)
end
UPauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = UPauli{UInt}(p)
# Base.:(==)(p::UPauli, q::UPauli) = p.string == q.string

struct Pauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    sign::C8
    qubits::Integer
    function Pauli{T,Q}(string::Integer, sign::Number) where {T,Q}
        _check_string_length(string, Q)
        sign in Set([1, -1, im, -im]) ? new{T,Q}(string, sign, Q) : throw(ArgumentError("Sign must be 1, -1, im, or -im"))
    end
end
Pauli(string::T, sign::Number, Q::Integer) where {T<:Unsigned} = Pauli{T,Q}(string, C8(sign))
Pauli{T,Q}(string::Unsigned) where {T,Q} = Pauli{T,Q}(string, 1)
Pauli(string::T, Q::Integer) where {T<:Unsigned} = Pauli{T,Q}(string, 1)
Pauli{T}(p::UPauli) where {T} = Pauli{T,p.qubits}(p.string, (im)^county(p))
Pauli(p::UPauli{T,Q}) where {T,Q} = Pauli{T}(p)
Pauli{T}(p::UPauli, sign::Number) where {T} = Pauli{T,p.qubits}(p.string, C8(sign))
Pauli(p::UPauli{T,Q}, sign::Number) where {T,Q} = Pauli{T,Q}(p.string, C8(sign))
Pauli{T}(p::Pauli) where {T} = Pauli{T,p.qubits}(p.string, p.sign)
Pauli(p::Pauli) = p
Pauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}) where {T} = Pauli(UPauli{T}(p), (im)^county(UPauli(p)))
Pauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = Pauli{UInt}(p)
Pauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) where {T} = Pauli(UPauli{T}(p), sign)
Pauli(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) = Pauli(UPauli(p), sign)
# Base.:(==)(p1::Pauli, p2::Pauli) = p1.string == p2.string & p1.sign == p2.sign
