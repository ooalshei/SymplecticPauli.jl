struct PauliList{T<:Unsigned,Q} <: AbstractVector{T}
    strings::Vector{T}
    qubits::Integer
    PauliList{T,Q}(vector) where {T,Q} = any(i -> i >= 4^Q, vector) ? throw(ArgumentError("Elements of the vector must be less than $(4^Q).")) : new{T,Q}(copy(vector), Q)
end

Base.size(v::PauliList) = (length(v.strings),)
Base.getindex(v::PauliList, i::Integer) = v.strings[i]
Base.getindex(v::PauliList, inds::UnitRange{Int}) = PauliList(v.strings[inds], v.qubits)
Base.setindex!(v::PauliList, value::T, i::Integer) where {T<:Unsigned} = setindex!(v.strings, value, i)
Base.resize!(v::PauliList, n::Integer) = resize!(v.strings, n)
Base.deleteat!(v::PauliList, i::Integer) = deleteat!(v.strings, i)
Base.deleteat!(v::PauliList, inds::UnitRange{Int}) = deleteat!(v.strings, inds)
# Base.similar(v::PauliList{T,Q}) where {T,Q} = PauliList{T,Q}(similar(v.strings))
# Base.similar(v::PauliList{T,Q}, m::Int) where {T,Q} = PauliList{T,Q}(similar(v.strings, m))
# Base.similar(v::PauliList{T,Q}, r::Base.OneTo) where {T,Q} = PauliList{T,Q}(similar(v.strings, r))
# Base.similar(v::PauliList{<:Unsigned,Q}, ::Type{T}) where {T<:Unsigned,Q} = PauliList{T,Q}(similar(v.strings, T))
# Base.similar(v::PauliList{<:Unsigned,Q}, ::Type{T}, m::Int) where {T<:Unsigned,Q} = PauliList{T,Q}(similar(v.strings, T, m))
Base.copy(v::PauliList) = PauliList(copy(v.strings), v.qubits)

toint(v::PauliList) = v.strings
PauliList(v::AbstractVector{T}, Q::Integer) where {T<:Unsigned} = PauliList{T,Q}(v)
PauliList{T}(v::AbstractVector{<:UPauli}) where {T,Q} = PauliList{T,Q}(toint.(v), v[1].qubits)
PauliList(v::AbstractVector{<:UPauli{T,Q}}) where {T,Q} = PauliList{T}(v)
PauliList{T}(v::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}) where {T} = PauliList{T,length(v[1])}(toint.(UPauli.(v)))
PauliList(v::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}) = PauliList{UInt}(v)
PauliList{T}(v::AbstractMatrix{<:Integer}) where {T} = PauliList{T}(eachcol(v))
PauliList(v::AbstractMatrix{<:Integer}) = PauliList{UInt}(v)
PauliList{T,Q}() where {T,Q} = PauliList{T,Q}(T[])
