struct PauliSentence{T<:Unsigned,N<:Number,Q} <: AbstractDict{T,N}
    sentence::Dict{T,N}
    qubits::Integer
    PauliSentence{T,N,Q}(sentence) where {T,N,Q} = (isempty(sentence) || maximum(keys(sentence)) < 4^Q) ? new{T,N,Q}(copy(sentence), Q) : throw(ArgumentError("String must not exceed $(4^Q - 1)."))
end

Base.show(io::IO, s::PauliSentence) = print(io, tostring(s))
Base.iterate(s::PauliSentence, i=1) = iterate(s.sentence, i)
Base.length(s::PauliSentence) = length(s.sentence)
Base.get(s::PauliSentence, key, default) = get(s.sentence, key, default)
Base.setindex!(s::PauliSentence, value, key) = setindex!(s.sentence, value, key)
Base.empty(::PauliSentence{T,N,Q}) where {T,N,Q} = PauliSentence{T,N,Q}(Dict{T,N}())
Base.delete!(s::PauliSentence, key) = delete!(s.sentence, key)

PauliSentence{T,N}(s::AbstractDict{<:Unsigned,<:Number}, Q::Integer) where {T,N} = PauliSentence{T,N,Q}(s)
PauliSentence(s::AbstractDict{T,N}, Q::Integer) where {T<:Unsigned,N<:Number} = PauliSentence{T,N,Q}(s)
function PauliSentence{T,N,Q}(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{<:Number}) where {T,N,Q}
    length(paulis) == length(coeffs) || throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    maximum(paulis) < 4^Q || throw(ArgumentError("Pauli string must not exceed $(4^Q - 1)."))
    # sentence = Dict(zip(paulis, coeffs))
    return PauliSentence{T,N,Q}(Dict(Pair.(paulis, coeffs)))
end
PauliSentence{T,N}(paulis::AbstractVector{<:Unsigned}, coeffs::AbstractVector{<:Number}, Q::Integer) where {T,N} = PauliSentence{T,N,Q}(paulis, coeffs)
PauliSentence(paulis::AbstractVector{T}, coeffs::AbstractVector{N}, Q::Integer) where {T<:Unsigned,N<:Number} = PauliSentence{T,N,Q}(paulis, coeffs)
PauliSentence{T,N}(paulis::PauliList, coeffs::AbstractVector{<:Number}) where {T,N} = PauliSentence{T,N}(UPauli.(paulis, paulis.qubits), coeffs)
PauliSentence(paulis::PauliList{T,Q}, coeffs::AbstractVector{N}) where {T,Q,N<:Number} = PauliSentence{T,promote_type(C8, N)}(paulis, coeffs)
PauliSentence{T,N}(paulis::AbstractVector{<:UPauli{<:Unsigned,Q}}, coeffs::AbstractVector{<:Number}) where {T,N,Q} = PauliSentence{T,N,Q}(map(x -> x.string, paulis), im .^ (county.(paulis)) .* coeffs)
PauliSentence(paulis::AbstractVector{UPauli{T,Q}}, coeffs::AbstractVector{N}) where {T,Q,N<:Number} = any(p -> isodd(county(p)), paulis) ? PauliSentence{T,promote_type(C8, N)}(paulis, coeffs) : PauliSentence{T,N}(paulis, coeffs)
function PauliSentence{T,N}(paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}, coeffs::AbstractVector{<:Number}) where {T,N}
    ps = Pauli.(paulis)
    c = copy(coeffs)
    for (i, p) in pairs(ps)
        c[i] *= p.sign
    end
    return PauliSentence{T,N,length(paulis[1])}(map(x -> x.string, ps), c)
end
PauliSentence(paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}, coeffs::AbstractVector{N}) where {N<:Number} = any(p -> isodd(county(p)), Pauli.(paulis)) ? PauliSentence{UInt,promote_type(C8, N)}(paulis, coeffs) : PauliSentence{UInt,N}(paulis, coeffs)
PauliSentence{T,N}(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{<:Number}) where {T,N} = PauliSentence{T,N}(eachcol(paulis), coeffs)
PauliSentence(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{N}) where {N<:Number} = PauliSentence{UInt,N}(eachcol(paulis), coeffs)
# PauliSentence{T,N}(sentence::AbstractDict{<:UPauli,<:Number}) where {T,N} = PauliSentence{T,N}(collect(keys(sentence)), collect(values(sentence)))
PauliSentence{T,N}(sentence::AbstractDict{<:Union{UPauli,AbstractString,AbstractVector{<:Integer}},<:Number}) where {T,N} = PauliSentence{T,N}(collect(keys(sentence)), collect(values(sentence)))
PauliSentence(sentence::AbstractDict{<:UPauli{T,Q},N}) where {T,N<:Number,Q} = PauliSentence{T,N}(sentence)
PauliSentence(sentence::AbstractDict{<:Union{AbstractString,AbstractVector{<:Integer}},N}) where {N<:Number} = PauliSentence{UInt,N}(sentence)
function PauliSentence{T,N}(m::AbstractMatrix{<:Number}) where {T,N}
    (x, y) = size(m)
    Q = log2(x)
    (isequal(x, y) & isinteger(Q)) || throw(ArgumentError("Matrix must be square and of size 2^Q x 2^Q for some integer Q."))
    sentence = PauliSentence(Dict{T,N}(), Int(Q))
    for i in 0:x^2-1
        c = conj(tr(tomatrix(T(i), Int(Q)) * m') / x)
        abs(c) > eps(Float64) && (sentence[T(i)] = c)
    end
    return sentence
end
PauliSentence(m::AbstractMatrix{<:Number}) = PauliSentence{UInt,ComplexF64}(m)
PauliSentence{T,N}(s::PauliSentence) where {T,N} = PauliSentence{T,N,s.qubits}(s.sentence)
PauliSentence(s::PauliSentence) = copy(s)
