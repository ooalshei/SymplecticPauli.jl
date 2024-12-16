struct PauliSentenceStyle <: Broadcast.AbstractArrayStyle{1} end
Base.BroadcastStyle(::Type{<:PauliSentence}) = PauliSentenceStyle()
Base.BroadcastStyle(s::PauliSentenceStyle, ::Broadcast.DefaultArrayStyle) = s
# Base.broadcastable(s::PauliSentence) = replace(s, missing => zero(eltype(s)))
Base.similar(s::Broadcast.Broadcasted{PauliSentenceStyle}, ::Type{T}) where {T} = PauliSentence(Vector{T}(undef, length(s)))

function Base.:+(A::PauliSentence, Bs::PauliSentence...)
    for B in Bs
        promote_shape(A, B) # check size compatibility
    end
    replace(Base.broadcast_preserving_zero_d(Base.:+, replace(A, missing => zero(eltype(A))), replace.(Bs, missing => Int8(0))...), 0 => missing)
end
function Base.:-(A::PauliSentence, B::PauliSentence)
    promote_shape(A, B) # check size compatibility
    replace(Base.broadcast_preserving_zero_d(Base.:-, replace(A, missing => zero(eltype(A))), replace(B, missing => zero(eltype(B)))), 0 => missing)
end