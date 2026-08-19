@doc raw"""
    PauliSentence(strings, coefficients, Q)
    PauliSentence(v::PauliList, coefficients)
    PauliSentence(m::AbstractMatrix)

An operator: Pauli strings with coefficients, behaving as an `AbstractDict` from string to
coefficient.

Coefficients multiply the *bare* embedding of their string (`X` ↦ ``σ_1``, `Z` ↦ ``σ_3``,
`Y` position ↦ the real antisymmetric ``σ_2``), so building one from a
[`PauliList`](@ref) folds in the ``i^{\#Y}`` phase that makes each term Hermitian — a real
coefficient vector gives a Hermitian operator. [`tomatrix`](@ref) inverts the
representation, and constructing from a matrix decomposes it into Paulis.

# Examples
```jldoctest
julia> H = PauliSentence(PauliList(["ZZ", "X-"]), [1.0, -0.5]);

julia> length(H), H.qubits
(2, 2)

julia> H[UPauli("ZZ").string]
1.0 + 0.0im

julia> tomatrix(H)
4×4 Matrix{ComplexF64}:
  1.0+0.0im   0.0+0.0im  -0.5+0.0im   0.0+0.0im
  0.0+0.0im  -1.0+0.0im   0.0+0.0im  -0.5+0.0im
 -0.5+0.0im   0.0+0.0im  -1.0+0.0im   0.0+0.0im
  0.0+0.0im  -0.5+0.0im   0.0+0.0im   1.0+0.0im
```
"""
struct PauliSentence{T<:Unsigned,N<:Number,Q} <: AbstractDict{T,N}
    sentence::Dict{T,N}
    qubits::Int
    function PauliSentence{T,N,Q}(
        sentence::AbstractDict;
        iscopy::Bool=true,
        check::Bool=true,
    ) where {T,N,Q}
        check && _check_keys(keys(sentence), Q)
        # `sentence isa Dict{T,N}` is a compile-time test: when the dictionary already has
        # the target element types it is stored (or copied) as is, otherwise the conversion
        # itself produces the fresh dictionary and no extra copy is made.
        if sentence isa Dict{T,N}
            return iscopy ? new{T,N,Q}(copy(sentence), Q) : new{T,N,Q}(sentence, Q)
        else
            return new{T,N,Q}(Dict{T,N}(sentence), Q)
        end
    end
end

_check_keys(ks, Q::Integer) = if all(k -> iszero(k >> (2 * Q)), ks)
    nothing
else
    throw(ArgumentError("String must not exceed $(big(4)^Q - 1)."))
end

Base.show(io::IO, s::PauliSentence) = print(io, tostring(s))
Base.iterate(s::PauliSentence, i=1) = iterate(s.sentence, i)
Base.length(s::PauliSentence) = length(s.sentence)
Base.isempty(s::PauliSentence) = isempty(s.sentence)
Base.get(s::PauliSentence, key, default) = get(s.sentence, key, default)
Base.getindex(s::PauliSentence, key) = getindex(s.sentence, key)
Base.setindex!(s::PauliSentence, value, key) = setindex!(s.sentence, value, key)
Base.haskey(s::PauliSentence, key) = haskey(s.sentence, key)
Base.keys(s::PauliSentence) = keys(s.sentence)
Base.values(s::PauliSentence) = values(s.sentence)
Base.sizehint!(s::PauliSentence, n::Integer) = (sizehint!(s.sentence, n); s)
Base.empty(::PauliSentence{T,N,Q}) where {T,N,Q} =
    PauliSentence{T,N,Q}(Dict{T,N}(), iscopy=false, check=false)
Base.copy(s::PauliSentence{T,N,Q}) where {T,N,Q} =
    PauliSentence{T,N,Q}(copy(s.sentence), iscopy=false, check=false)
Base.delete!(s::PauliSentence, key) = (delete!(s.sentence, key); s)
Base.pop!(s::PauliSentence, key, default) = pop!(s.sentence, key, default)
Base.pop!(s::PauliSentence, key) = pop!(s.sentence, key)
Base.filter!(f, s::PauliSentence) = (filter!(f, s.sentence); s)
Base.filter(f, s::PauliSentence{T,N,Q}) where {T,N,Q} =
    PauliSentence{T,N,Q}(filter(f, s.sentence), iscopy=false, check=false)

PauliSentence{T,N}(
    s::AbstractDict{<:Unsigned,<:Number},
    Q::Integer;
    iscopy=true,
    check=true,
) where {T,N} = PauliSentence{T,N,Q}(s, iscopy=iscopy, check=check)
PauliSentence(
    s::AbstractDict{T,N},
    Q::Integer;
    iscopy=true,
    check=true,
) where {T<:Unsigned,N<:Number} = PauliSentence{T,N,Q}(s, iscopy=iscopy, check=check)
function PauliSentence{T,N,Q}(
    paulis::AbstractVector{<:Unsigned},
    coeffs::AbstractVector{<:Number},
) where {T,N,Q}
    length(paulis) == length(coeffs) ||
        throw(DimensionMismatch("Length of paulis and coeffs must be the same."))
    _check_keys(paulis, Q)
    return PauliSentence{T,N,Q}(Dict{T,N}(Pair.(paulis, coeffs)), iscopy=false, check=false)
end
PauliSentence{T,N}(
    paulis::AbstractVector{<:Unsigned},
    coeffs::AbstractVector{<:Number},
    Q::Integer,
) where {T,N} = PauliSentence{T,N,Q}(paulis, coeffs)
PauliSentence(
    paulis::AbstractVector{T},
    coeffs::AbstractVector{N},
    Q::Integer,
) where {T<:Unsigned,N<:Number} = PauliSentence{T,N,Q}(paulis, coeffs)
PauliSentence{T,N}(paulis::PauliList, coeffs::AbstractVector{<:Number}) where {T,N} =
    PauliSentence{T,N}(UPauli.(paulis, paulis.qubits), coeffs)
PauliSentence(paulis::PauliList{T,Q}, coeffs::AbstractVector{N}) where {T,Q,N<:Number} =
    PauliSentence{T,promote_type(C8, N)}(paulis, coeffs)
PauliSentence{T,N}(
    paulis::AbstractVector{<:UPauli{<:Unsigned,Q}},
    coeffs::AbstractVector{<:Number},
) where {T,N,Q} =
    PauliSentence{T,N,Q}(map(x -> x.string, paulis), im .^ (county.(paulis)) .* coeffs)
PauliSentence(
    paulis::AbstractVector{UPauli{T,Q}},
    coeffs::AbstractVector{N},
) where {T,Q,N<:Number} =
    if any(p -> isodd(county(p)), paulis)
        PauliSentence{T,promote_type(C8, N)}(paulis, coeffs)
    else
        PauliSentence{T,N}(paulis, coeffs)
    end
function PauliSentence{T,N}(
    paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}},
    coeffs::AbstractVector{<:Number},
) where {T,N}
    ps = Pauli{T}.(paulis)
    # The convention phase is `i` on every string with an odd number of `Y`s, so it has to
    # widen the coefficients rather than be written back into a copy of them: a real vector
    # in gives a complex vector out.
    c = getfield.(ps, :sign) .* coeffs
    return PauliSentence{T,N,length(paulis[1])}(map(x -> x.string, ps), c)
end
PauliSentence(
    paulis::AbstractVector{<:Union{AbstractString,AbstractVector{<:Integer}}},
    coeffs::AbstractVector{N},
) where {N<:Number} =
    if any(p -> isodd(county(p)), Pauli.(paulis))
        PauliSentence{UInt,promote_type(C8, N)}(paulis, coeffs)
    else
        PauliSentence{UInt,N}(paulis, coeffs)
    end
PauliSentence{T,N}(
    paulis::AbstractMatrix{<:Integer},
    coeffs::AbstractVector{<:Number},
) where {T,N} = PauliSentence{T,N}(eachcol(paulis), coeffs)
# Through the vector method rather than straight to `PauliSentence{UInt,N}`: the coefficient
# type has to widen when the convention phase of a string is imaginary, and that is where
# the widening lives.
PauliSentence(paulis::AbstractMatrix{<:Integer}, coeffs::AbstractVector{<:Number}) =
    PauliSentence(eachcol(paulis), coeffs)

PauliSentence{T,N}(
    sentence::AbstractDict{
        <:Union{UPauli,AbstractString,AbstractVector{<:Integer}},
        <:Number,
    },
) where {T,N} = PauliSentence{T,N}(collect(keys(sentence)), collect(values(sentence)))
PauliSentence(
    sentence::AbstractDict{
        <:Union{UPauli,AbstractString,AbstractVector{<:Integer}},
        <:Number,
    },
) = PauliSentence(collect(keys(sentence)), collect(values(sentence)))
function PauliSentence{T,N}(m::AbstractMatrix{<:Number}) where {T,N}
    (x, y) = size(m)
    Q = log2(x)
    (isequal(x, y) & isinteger(Q)) || throw(
        ArgumentError("Matrix must be square and of size 2^Q x 2^Q for some integer Q."),
    )
    sentence = PauliSentence(Dict{T,N}(), Int(Q), iscopy=false)
    for i in 0:(x^2-1)
        c = conj(tr(tomatrix(T(i), Int(Q)) * m') / x)
        abs(c) > eps(Float64) && (sentence[T(i)] = c)
    end
    return sentence
end
PauliSentence(m::AbstractMatrix{<:Number}) = PauliSentence{UInt,ComplexF64}(m)
PauliSentence{T,N}(s::PauliSentence{<:Unsigned,<:Number,Q}) where {T,N,Q} =
    PauliSentence{T,N,Q}(s.sentence, check=false)
PauliSentence(s::PauliSentence) = copy(s)
