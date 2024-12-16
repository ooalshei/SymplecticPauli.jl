const NumOrMiss = Union{Number,Missing}

struct PauliSentence{T<:NumOrMiss} <: AbstractVector{T}
    sentence::Vector{T}
    PauliSentence{T}(sentence) where {T} = isinteger(log2(length(sentence)) / 2) ? new{T}(copy(sentence)) : throw(ArgumentError("Length of PauliSentence must be a power of 4."))
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

PauliSentence(coeffs::AbstractVector{T}) where {T<:NumOrMiss} = PauliSentence{T}(Vector(coeffs))
function PauliSentence{T}(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{<:NumOrMiss}, qubits::Integer) where {T}
    # length(paulis) == length(coeffs) || throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    maximum(paulis) < 4^qubits || throw(ArgumentError("Pauli string must not exceed $(4^qubits - 1)."))
    sentence = Vector{Union{T,Missing}}(missing, 4^qubits)
    sentence[paulis.+1] = coeffs
    return PauliSentence{T}(sentence)
end
PauliSentence(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{T}, qubits::Integer) where {T<:NumOrMiss} = PauliSentence{T}(paulis, coeffs, qubits)
PauliSentence{T}(paulis::AbstractVector{<:Pauli{<:Unsigned,Q}}, coeffs::AbstractVector{<:NumOrMiss}) where {T,Q} = PauliSentence{T}([p.string for p in paulis], coeffs, Q)
PauliSentence(paulis::AbstractVector{<:Pauli}, coeffs::AbstractVector{T}) where {T<:NumOrMiss} = PauliSentence{T}(paulis, coeffs)
PauliSentence{T}(paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}, coeffs::AbstractVector{<:NumOrMiss}) where {T} = PauliSentence{T}(Pauli.(paulis), coeffs)
PauliSentence(paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}, coeffs::AbstractVector{T}) where {T<:NumOrMiss} = PauliSentence{T}(paulis, coeffs)
PauliSentence{T}(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{<:NumOrMiss}) where {T} = PauliSentence{T}(eachcol(paulis), coeffs)
PauliSentence(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{T}) where {T<:NumOrMiss} = PauliSentence{T}(eachcol(paulis), coeffs)
PauliSentence{T}(paulis::AbstractDict{<:Union{Unsigned,Pauli,AbstractString,AbstractVector{<:Integer}},<:NumOrMiss}) where {T} = PauliSentence{T}(keys(paulis), values(paulis))
PauliSentence(paulis::AbstractDict{<:Union{Unsigned,Pauli,AbstractString,AbstractVector{<:Integer}},T}) where {T<:NumOrMiss} = PauliSentence{T}(keys(paulis), values(paulis))
PauliSentence{T}(s::PauliSentence) where {T} = PauliSentence{T}(s.sentence)
PauliSentence(s::PauliSentence{T}) where {T} = copy(s)

function tostring(s::PauliSentence{T})::Dict{String,T} where {T}
    qubits = Int(log2(length(s)) / 2)
    result = Dict{String,T}()
    for (key, value) in pairs(skipmissing(s))
        result[tostring(Pauli(UInt(key), qubits))] = value
    end
    return result
end