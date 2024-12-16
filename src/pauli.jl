abstract type AbstractPauli{T<:Unsigned,Q} end
Base.copy(p::AbstractPauli) = p
Base.show(io::IO, p::AbstractPauli{T}) where {T} = print(io, "Pauli{$T}(", tostring(p), ")")

struct Pauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    qubits::Integer
    Pauli{T,Q}(string) where {T,Q} = string > (4^Q - 1) ? throw(ArgumentError("String must not exceed $(4^Q - 1).")) : new{T,Q}(string, Q)
end
Pauli(string::T, Q::Integer) where {T<:Unsigned} = Pauli{T,Q}(string)
# Pauli{T,Q}(p::AbstractPauli) where {T,Q} = Pauli{T,Q}(p.string)
Pauli{T}(p::AbstractPauli) where {T} = Pauli{T,p.qubits}(p.string)
Pauli(p::AbstractPauli{T,Q}) where {T,Q} = Pauli{T,Q}(p.string)
function Pauli{T}(p::AbstractString) where {T}
    unique(p) ⊆ ['X', 'Y', 'Z', 'I', '-'] || throw(ArgumentError("String must contain only 'X', 'Y', 'Z', 'I', or '-'."))
    Q = length(p)
    Base.hastypemax(T) && typemax(T) < 4^Q && throw(ArgumentError("String length cannot exceed $(count_ones(typemax(T)) ÷ 2)). Consider using a larger unsigned integer type."))
    number = T(0)
    for char in p
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
    return Pauli{T,Q}(number)
end
function Pauli{T}(p::AbstractVector{<:Integer}) where {T}
    unique(p) ⊆ [1, 2, 3, 4] || throw(ArgumentError("Array must contain only 1, 2, 3, or 4."))
    Q = length(p)
    Base.hastypemax(T) && typemax(T) < 4^Q && throw(ArgumentError("Array length cannot exceed $(count_ones(typemax(T)) ÷ 2)). Consider using a larger unsigned integer type."))
    number = T(0)
    for ind in p
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
    Pauli{T,Q}(number)
end
Pauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = Pauli{UInt}(p)
# Base.:(==)(p::Pauli, q::Pauli) = p.string == q.string

function tostring(p::Pauli)::String
    result = ""
    string = digits(p.string, base=2, pad=2*p.qubits)
    for i in 1:p.qubits
        if string[i] == string[i + p.qubits] == 1
            result = "Y" * result
        elseif string[i] == 1
            result = "X" * result
        elseif string[i + p.qubits] == 1
            result = "Z" * result
        else
            result = "-" * result
        end
    end
    return result
end

struct SignedPauli{T<:Unsigned,Q} <: AbstractPauli{T,Q}
    string::T
    sign::C8
    qubits::Integer
    function SignedPauli{T,Q}(string, sign) where {T,Q}
        if string > (4^Q - 1)
            throw(ArgumentError("String must be of length $(4^Q - 1) or less."))
        elseif !(sign in Set([1, -1, im, -im]))
            throw(ArgumentError("Sign must be 1, -1, im, or -im"))
        else
            new{T,Q}(string, sign, Q)
        end
    end
end
SignedPauli(string::T, sign::Number, Q::Integer) where {T<:Unsigned} = SignedPauli{T,Q}(string, C8(sign))
SignedPauli{T,Q}(string::Unsigned) where {T,Q} = SignedPauli{T,Q}(string, C8(1))
SignedPauli(string::T, Q::Integer) where {T<:Unsigned} = SignedPauli{T,Q}(string, C8(1))
SignedPauli{T}(p::Pauli) where {T} = SignedPauli{T,p.qubits}(p.string, C8(1))
SignedPauli(p::Pauli{T,Q}) where {T,Q} = SignedPauli{T,Q}(p.string, C8(1))
SignedPauli{T}(p::Pauli, sign::Number) where {T} = SignedPauli{T,p.qubits}(p.string, C8(sign))
SignedPauli(p::Pauli{T,Q}, sign::Number) where {T,Q} = SignedPauli{T,Q}(p.string, C8(sign))
SignedPauli{T}(p::SignedPauli) where {T} = SignedPauli{T,p.qubits}(p.string, p.sign)
SignedPauli(p::SignedPauli) = p
SignedPauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}) where {T} = SignedPauli(Pauli{T}(p))
SignedPauli(p::Union{AbstractString,AbstractVector{<:Integer}}) = SignedPauli(Pauli(p))
SignedPauli{T}(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) where {T} = SignedPauli(Pauli{T}(p), sign)
SignedPauli(p::Union{AbstractString,AbstractVector{<:Integer}}, sign::Number) = SignedPauli(Pauli(p), sign)
# Base.:(==)(p1::SignedPauli, p2::SignedPauli) = p1.string == p2.string & p1.sign == p2.sign

function tostring(p::SignedPauli)::String
    p.sign == 1 && return "(+)" * tostring(Pauli(p))
    p.sign == -1 && return "(-)" * tostring(Pauli(p))
    p.sign == 1im && return "(i)" * tostring(Pauli(p))
    p.sign == -1im && return "(-i)" * tostring(Pauli(p))
end