@doc raw"""
    PauliList(strings::AbstractVector{<:Unsigned}, Q)
    PauliList(text::AbstractVector{<:AbstractString})
    PauliList(indices::AbstractVector{<:AbstractVector{<:Integer}})
    PauliList(indices::AbstractMatrix{<:Integer})
    PauliList(paulis::AbstractVector{<:UPauli})
    PauliList{T,Q}()
    PauliList{T,Q}(undef, n)

An ordered list of Pauli strings on `Q` qubits — an algebra, a generator set, a subalgebra.

Behaves as an `AbstractVector` of the raw bit strings, so `push!`, `append!`, `deleteat!`,
indexing and iteration all work; `v[i:j]` returns another `PauliList`. Use
[`tostring`](@ref) to read one back as text.

Elements are the bare `T` bit patterns rather than [`UPauli`](@ref) objects: the qubit count
is a property of the list, stored once, so a list of a million strings is a million
integers. Wrap an element as `UPauli(v[i], v.qubits)` when a standalone string is needed.
Only `Q` and the storage type `T` — chosen as in [`UPauli`](@ref), `UInt` by default — are
part of the type; a list built from a matrix of `1`–`4` indices reads one string per
*column*.

Two keywords control what the constructor does with the vector it is handed. `check=false`
skips verifying that every element fits in `2Q` bits, and `iscopy=false` takes ownership of
the vector instead of copying it — both worth it in a hot loop that has just built the
vector itself, and both unsafe otherwise.

# Examples
```jldoctest
julia> v = PauliList(["XX-", "-YZ"]);

julia> length(v), v.qubits
(2, 3)

julia> tostring(v)
2-element Vector{String}:
 "XX-"
 "-YZ"

julia> tostring(push!(copy(v), UPauli("ZZZ").string))
3-element Vector{String}:
 "XX-"
 "-YZ"
 "ZZZ"

julia> tostring(PauliList([2 1; 2 3; 1 4]))     # one string per column
2-element Vector{String}:
 "XX-"
 "-YZ"

julia> tostring(UPauli(v[2], v.qubits))         # an element as a standalone string
"-YZ"
```
"""
struct PauliList{T<:Unsigned,Q} <: AbstractVector{T}
    strings::Vector{T}
    qubits::Int
    function PauliList{T,Q}(vector; iscopy=true, check=true) where {T,Q}
        if check && any(i -> !iszero(i >> (2 * Q)), vector)
            throw(ArgumentError("Elements of the vector must be less than $(big(4)^Q)."))
        elseif iscopy
            new{T,Q}(copy(vector), Q)
        else
            new{T,Q}(vector, Q)
        end
    end
end
Base.show(io::IO, v::PauliList) = print(io, tostring(v))
Base.size(v::PauliList) = (length(v.strings),)
Base.length(v::PauliList) = length(v.strings)
Base.IndexStyle(::Type{<:PauliList}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::PauliList, i::Integer) = v.strings[i]
Base.copy(v::PauliList{T,Q}) where {T,Q} =
    PauliList{T,Q}(copy(v.strings), iscopy=false, check=false)
Base.push!(v::PauliList, item) = (push!(v.strings, item); v)
Base.append!(v::PauliList, items) = (append!(v.strings, items); v)
Base.empty!(v::PauliList) = (empty!(v.strings); v)
Base.sizehint!(v::PauliList, n::Integer) = (sizehint!(v.strings, n); v)
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
Base.pop!(v::PauliList) = pop!(v.strings)
Base.popat!(v::PauliList, i) = popat!(v.strings, i)
Base.popfirst!(v::PauliList) = popfirst!(v.strings)
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
) where {T} =
    PauliList{T,length(v[1])}(map(x -> UPauli{T}(x).string, v), iscopy=false, check=false)

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
