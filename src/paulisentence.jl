struct PauliSentence{T<:Number} <: AbstractVector{T}
    sentence::Vector{<:Union{T,Missing}}
    PauliSentence{T}(sentence) where {T} =  new{T}(copy(sentence))
end
Base.axes(s::PauliSentence) = (ZeroTo(length(s)),)
Base.size(s::PauliSentence) = size(s.sentence)
Base.getindex(s::PauliSentence, i::Integer) = s.sentence[i+1]
Base.getindex(s::PauliSentence, p::Pauli) = s.sentence[p.string+1]
# Base.unsafe_getindex(s::PauliSentence, i::Integer) = Base.unsafe_getindex(s.sentence, i + 1)
# Base.unsafe_getindex(s::PauliSentence, p::Pauli) = Base.unsafe_getindex(s.sentence, p.string + 1)
Base.setindex!(s::PauliSentence, value, i::Integer) = (s.sentence[i+1] = value)
Base.setindex!(s::PauliSentence, value, p::Pauli) = (s.sentence[p.string+1] = value)
Base.keys(s::PauliSentence) = ZeroTo(length(s))
Base.show(io::IO, s::PauliSentence) = print(io, tostring(s))
Base.similar(s::PauliSentence) = PauliSentence(similar(s.sentence))
# Base.similar(s::PauliSentence, T::Type) = PauliSentence(similar(s.sentence, T))
Base.similar(s::PauliSentence, m::Int) = PauliSentence(similar(s.sentence, m))
Base.similar(s::PauliSentence, r::ZeroTo) = PauliSentence(similar(s.sentence, r))
Base.similar(s::PauliSentence, ::Type{T}) where {T} = PauliSentence(similar(s.sentence, T))
Base.similar(s::PauliSentence, ::Type{T}, shape::Tuple{ZeroTo,Vararg{ZeroTo}}) where {T} = PauliSentence(similar(s.sentence, T, shape))
Base.similar(s::PauliSentence, ::Type{T}, m::Int) where {T} = PauliSentence(similar(s.sentence, T, m))
# Base.similar(::Type{T}, shape::Tuple{ZeroTo,Vararg{ZeroTo}}) where {T<:PauliSentence} = similar(T, Base.to_shape(shape))

PauliSentence(coeffs::AbstractVector{<:Union{Missing,T}}) where {T<:Number} = PauliSentence{T}(Vector{Union{Missing,T}}(coeffs))
function PauliSentence{T}(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{<:Number}, qubits::Integer) where {T}
    length(paulis) == length(coeffs) || throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    max(paulis) < 4^qubits || throw(ArgumentError("Pauli string must not exceed $(4^qubits - 1)."))
    sentence = Vector{<:Union{T,Missing}}(missing, 4^length(qubits))
    for (p, coeff) in zip(paulis, coeffs)
        sentence[p+1] = coeff
    end
    return PauliSentence{T}(sentence)
end
PauliSentence(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{T}, qubits::Integer) where {T<:Number} = PauliSentence{T}(paulis, coeffs, qubits)
function PauliSentence{T}(paulis::AbstractVector{<:Pauli}, coeffs::AbstractVector{<:Number}) where {T}
    length(paulis) == length(coeffs) || throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    sentence = Vector{<:Union{T,Missing}}(missing, 4^length(paulis[1].qubits))
    for (p, coeff) in zip(paulis, coeffs)
        sentence[p.string+1] = coeff
    end
    return PauliSentence{T}(sentence)
end
PauliSentence(paulis::AbstractVector{<:Pauli}, coeffs::AbstractVector{T}) where {T<:Number} = PauliSentence{T}(paulis, coeffs)
function PauliSentence{T}(paulis::AbstractVector{<:Union{<:AbstractString,<:AbstractVector{<:Integer}}}, coeffs::AbstractVector{<:Number}) where {T}
    length(paulis) == length(coeffs) || throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    sentence = Vector{<:Union{T,Missing}}(missing, 4^length(paulis[1]))
    for (p, coeff) in zip(paulis, coeffs)
        sentence[Pauli(p).string+1] = coeff
    end
    return PauliSentence{T}(sentence)
end
PauliSentence(paulis::AbstractVector{<:Union{<:AbstractString,<:AbstractVector{<:Integer}}}, coeffs::AbstractVector{T}) where {T<:Number} = PauliSentence{T}(paulis, coeffs)
PauliSentence{T}(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{<:Number}) where {T} = PauliSentence{T}(eachcol(paulis), coeffs)
PauliSentence(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{T}) where {T<:Number} = PauliSentence{T}(eachcol(paulis), coeffs)
PauliSentence{T}(paulis::AbstractDict{<:Union{<:Unsigned,<:Pauli,<:AbstractString,<:AbstractVector{<:Integer}},<:Number}) where {T} = PauliSentence{T}(keys(paulis), values(paulis))
PauliSentence(paulis::AbstractDict{<:Union{<:Unsigned,<:Pauli,<:AbstractString,<:AbstractVector{<:Integer}},T}) where {T<:Number} = PauliSentence{T}(keys(paulis), values(paulis))
PauliSentence{T}(s::PauliSentence) where {T} = PauliSentence{T}(s.sentence)
PauliSentence(s::PauliSentence{T}) where {T} = copy(s)

function tostring(s::PauliSentence{T})::Dict{String,T} where {T}
    qubits = Integer(log(4, length(s)))
    result = Dict{String,T}()
    for (key, value) in pairs(skipmissing(s))
        result[tostring(Pauli(UInt(key), qubits))] = value
    end
    return result
end
