struct PauliList{T<:Unsigned,Q} <: AbstractVector{T}
    strings::Vector{T}
    qubits::Integer
    function PauliList{T,Q}(vector; iscopy=true, check=true) where {T,Q}
        if check && any(i -> i >= 4^Q, vector)
            throw(ArgumentError("Elements of the vector must be less than $(4^Q)."))
        elseif iscopy
            new{T,Q}(copy(vector), Q)
        else
            new{T,Q}(vector, Q)
        end
    end
end
Base.show(io::IO, v::PauliList) = print(io, tostring(v))
Base.size(v::PauliList) = (length(v.strings),)
Base.getindex(v::PauliList, i::Integer) = v.strings[i]
Base.getindex(v::PauliList, inds::UnitRange{Int}) =
    PauliList(v.strings[inds], v.qubits, iscopy=false, check=false)

Base.setindex!(v::PauliList, value::T, i::Integer) where {T<:Unsigned} =
    setindex!(v.strings, value, i)

Base.resize!(v::PauliList, n::Integer) =
    PauliList(resize!(v.strings, n), v.qubits, iscopy=false, check=false)

Base.unique(v::PauliList) =
    PauliList(unique(v.strings), v.qubits, iscopy=false, check=false)

Base.deleteat!(v::PauliList, i) =
    PauliList(deleteat!(v.strings, i), v.qubits, iscopy=false, check=false)

# Base.deleteat!(v::PauliList, inds) = deleteat!(v.strings, inds)
Base.popat!(v::PauliList, i) = popat!(v.strings, i)
Base.similar(v::PauliList{T,Q}) where {T,Q} =
    PauliList{T,Q}(similar(v.strings), iscopy=false, check=false)

Base.similar(v::PauliList{T,Q}, m::Int) where {T,Q} =
    PauliList{T,Q}(similar(v.strings, m), iscopy=false, check=false)

Base.similar(v::PauliList{T,Q}, r::Base.OneTo) where {T,Q} =
    PauliList{T,Q}(similar(v.strings, r), iscopy=false, check=false)

Base.similar(v::PauliList{<:Unsigned,Q}, ::Type{T}) where {T<:Unsigned,Q} =
    PauliList{T,Q}(similar(v.strings, T), iscopy=false, check=false)

Base.similar(v::PauliList{<:Unsigned,Q}, ::Type{T}, m::Int) where {T<:Unsigned,Q} =
    PauliList{T,Q}(similar(v.strings, T, m), iscopy=false, check=false)
# Base.copy(v::PauliList) = PauliList(v.strings, v.qubits, check=false)

# toint(v::PauliList) = v.strings
PauliList(v::AbstractVector{T}, Q::Integer; iscopy=true, check=true) where {T<:Unsigned} =
    PauliList{T,Q}(v, iscopy=iscopy, check=check)

PauliList{T}(v::AbstractVector{<:UPauli}) where {T} =
    PauliList{T,v[1].qubits}(map(x -> x.string, v), iscopy=false, check=false)

PauliList(v::AbstractVector{<:UPauli{T,Q}}) where {T,Q} = PauliList{T}(v)
PauliList{T}(
    v::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}},
) where {T} = PauliList{T,length(v[1])}(
    map(x -> (x |> UPauli |> (y -> y.string)), v),
    iscopy=false,
    check=false,
)

PauliList(v::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}}) =
    PauliList{UInt}(v)

PauliList{T}(v::AbstractMatrix{<:Integer}) where {T} = PauliList{T}(eachcol(v))
PauliList(v::AbstractMatrix{<:Integer}) = PauliList{UInt}(v)
PauliList{T,Q}() where {T,Q} = PauliList{T,Q}(T[], iscopy=false, check=false)
PauliList{T,Q}(::UndefInitializer, n::Integer) where {T,Q} =
    PauliList{T,Q}(Vector{T}(undef, n), iscopy=false, check=false)

PauliList{T}(::UndefInitializer, n::Integer, Q::Integer) where {T} =
    PauliList{T,Q}(Vector{T}(undef, n), iscopy=false, check=false)

PauliList(::UndefInitializer, n::Integer, Q::Integer) =
    PauliList{UInt,Q}(Vector{UInt}(undef, n), iscopy=false, check=false)
