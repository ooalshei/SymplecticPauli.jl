struct ZeroTo{T<:Integer} <: AbstractUnitRange{T}
    stop::T
    ZeroTo{T}(stop) where {T<:Integer} = new(max(zero(T) + 1, stop))
end
ZeroTo(stop::T) where {T<:Integer} = ZeroTo{T}(stop)
Base.length(r::ZeroTo) = Integer(r.stop - zero(r.stop))
Base.first(::ZeroTo{T}) where {T} = zero(T)
Base.last(r::ZeroTo{T}) where {T} = convert(T, r.stop - 1)
Base.unsafe_getindex(::ZeroTo{T}, i::Integer) where {T} = convert(T, i - 1)
Base.show(io::IO, r::ZeroTo) = print(io, "ZeroTo(", r.stop, ")")
Base.:(==)(r::ZeroTo, s::ZeroTo) = last(r) == last(s)
Base.intersect(r::ZeroTo, s::ZeroTo) = ZeroTo(min(r.stop, s.stop))
Base.union(r::ZeroTo, s::ZeroTo) = ZeroTo(max(r.stop, s.stop))
Base.issubset(r::ZeroTo, s::ZeroTo) = r.stop <= s.stop
Base.to_shape(r::ZeroTo) = Int(r.stop)
# Base.similar(A::AbstractArray, T::Type, shape::Tuple{ZeroTo,Vararg{ZeroTo}}) = similar(A, T, Base.to_shape(shape))
# Base.similar(::Type{T}, shape::Tuple{ZeroTo,Vararg{ZeroTo}}) where {T<:AbstractArray} = similar(T, Base.to_shape(shape))
