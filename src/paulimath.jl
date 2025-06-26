# const SIGNS = C8[1, 1im, -1, -1im]

function Base.:*(scalar::Number, p::AbstractPauli)
    typeof(p) <: UPauli && return Pauli(p, im^county(p) * C8(scalar))
    return Pauli(p.string, p.sign * C8(scalar))
end

unsigned_prod(p::UPauli{<:Unsigned,Q}, q::UPauli{<:Unsigned,Q}) where {Q} = UPauli(p.string ⊻ q.string, Q)
Base.xor(p::UPauli, q::UPauli) = unsigned_prod(p, q)

function com(p::AbstractPauli{T,Q}, q::AbstractPauli{U,Q})::UPauli{<:Unsigned,Q} where {T,U,Q}
    S = promote_type(T, U)
    # mask::S = 2^Q - 1
    # bmask = ~amask
    # a1 = p.string & mask
    b1 = p.string >> Q
    # a2 = q.string & mask
    b2 = q.string >> Q
    (isodd(count_ones(p.string & b2)) ⊻ isodd(count_ones(b1 & q.string))) && (return p ⊻ q)
    return UPauli(S(0), Q)
end
function com(p::Unsigned, q::Unsigned, Q::Integer)
    # mask::UInt = 2^Q - 1
    # bmask = ~amask
    # a1 = p & mask
    b1 = p >> Q
    # a2 = q & mask
    b2 = q >> Q
    (isodd(count_ones(p & b2)) ⊻ isodd(count_ones(b1 & q))) && (return p ⊻ q)
    return UInt(0)
end
Base.:^(p::AbstractPauli, q::AbstractPauli) = com(p, q)

function Base.:*(p::AbstractPauli{<:Unsigned,Q}, q::AbstractPauli{<:Unsigned,Q})::Pauli{<:Unsigned,Q} where {Q}
    # S = promote_type(T, U)
    # mask::S = 2^Q - 1
    # bmask = ~amask
    # a1 = p.string & mask
    b1 = p.string >> Q
    # a2 = q.string & mask
    # b2 = q.string >> Q

    # overlap = count_ones((a1 | b1) & (a2 | b2))
    # signstring = ((b1 & a1) & (~b2 & a2)) | ((b1 & ~a1) & (b2 & a2)) | ((~b1 & a1) & (b2 & ~a2))
    # equalUPaulis = (a1 & a2 & ~b1 & ~b2) | (a1 & b1 & a2 & b2) | (~a1 & ~a2 & b1 & b2)
    # ind = (overlap - count_ones(equalUPaulis)) % 4 + 1
    # sign = SIGNS[ind] * Pauli(p).sign * Pauli(q).sign
    # isodd(count_ones(signstring)) && (sign *= C8(-1))
    result = p.string ⊻ q.string
    sign::C8 = Pauli(p).sign * Pauli(q).sign * (-1)^count_ones(b1 & q)
    return Pauli(result, sign, Q)
end
function _symplectic_prod(p::T, q::T, Q::Integer) where {T<:Unsigned}
    # mask::T = 2^Q - 1
    # bmask = ~amask
    # a1 = p & amask
    b1 = p >> Q
    # a2 = q & mask
    # b2 = (q & bmask) >> Q

    # overlap = count_ones((a1 | b1) & (a2 | b2))
    # signstring = ((b1 & a1) & (~b2 & a2)) | ((b1 & ~a1) & (b2 & a2)) | ((~b1 & a1) & (b2 & ~a2))
    # equalUPaulis = (a1 & a2 & ~b1 & ~b2) | (a1 & b1 & a2 & b2) | (~a1 & ~a2 & b1 & b2)
    # ind = (overlap - count_ones(equalUPaulis)) % 4 + 1
    # sign = SIGNS[ind]
    # isodd(count_ones(signstring)) && (sign *= C8(-1))
    # a = result & amask
    # b = (result & bmask) >> Q
    return p ⊻ q => C8(-1)^count_ones(b1 & q)
end

Base.cis(p::UPauli{T,Q}, θ::Real) where {T,Q} = PauliSentence{T,ComplexF64,Q}([T(0), p.string], [cos(θ), im^(county(p) + 1) * sin(θ)])

Base.:+(s::PauliSentence{<:Unsigned,<:Number,Q}, r::PauliSentence{<:Unsigned,<:Number,Q}...) where {Q} = PauliSentence(mergewith(+, s, r...), Q)
Base.:-(s::PauliSentence{<:Unsigned,<:Number,Q}, r::PauliSentence{<:Unsigned,<:Number,Q}) where {Q} = PauliSentence(mergewith(-, s, r), Q)

function Base.:*(c::N, s::PauliSentence{T,Ns,Q}) where {T,N,Ns,Q}
    result = PauliSentence{T,promote_type(N, Ns)}(s)
    for key in keys(result)
        result[key] *= c
    end
    return result
end
Base.:*(s::PauliSentence, c::Number) = c * s

function Base.:*(s::PauliSentence{Ts,<:Number,Q}, r::PauliSentence{Tr,<:Number,Q}) where {Ts,Tr,Q}
    result = PauliSentence(Dict{promote_type(Ts, Tr),ComplexF64}(), Q)
    for (key1, value1) in s
        for (key2, value2) in r
            string = _symplectic_prod(key1, key2, Q)
            haskey(result, string.first) ? result[string.first] += string.second * value1 * value2 : result[key] = string.second * value1 * value2
        end
    end
    return result
end

Base.:^(s::PauliSentence{T,N,Q}, x::Integer) where {T,N,Q} = x >= 0 ? prod(s for _ in 1:x; init=PauliSentence{T,N,Q}(T[0], [1])) : throw(ArgumentError("Exponent must be a non-negative integer."))

function ad(s::PauliSentence{Ts,<:Number,Q}, generator::UPauli{T,Q}, cosine::Real, sine::Real; atol::Real=0) where {Ts,T,Q}
    result = PauliSentence{promote_type(T, Ts),ComplexF64,Q}(s)
    (iszero(sine) | iszero(generator.string)) && return result
    modsine = (im)^(county(generator) + 1) * sine
    for (key, value) in s
        string = _symplectic_prod(generator.string, key, Q)
        if string.second != _symplectic_prod(key, generator.string, Q).second
            result[key] += (cosine - 1) * value
            haskey(result, string.first) ? result[string.first] += modsine * string.second * value : result[string.first] = modsine * string.second * value
        end
    end
    return filter!(p -> (abs(p.second) > atol), result)
end
ad(s::PauliSentence, generator::UPauli, angle::Real; atol::Real=0) = ad(s, generator, cos(2 * angle), sin(2 * angle), atol=atol)

function ad!(s::PauliSentence{<:Unsigned,ComplexF64,Q}, generator::UPauli{<:Unsigned,Q}, cosine::Real, sine::Real; atol::Real=0) where {Q}
    (iszero(sine) | iszero(generator.string)) && return s
    modsine = (im)^(county(generator) + 1) * sine
    keylist = keys(s)
    valuelist = values(s)
    for (key, value) in zip(keylist, valuelist)
        string = _symplectic_prod(generator.string, key, Q)
        if string.second != _symplectic_prod(key, generator.string, Q).second
            s[key] += (cosine - 1) * value
            haskey(s, string.first) ? s[string.first] += modsine * string.second * value : s[string.first] = modsine * string.second * value
        end
    end
    return filter!(p -> (abs(p.second) > atol), s)
end
ad!(s::PauliSentence, generator::UPauli, angle::Real; atol::Real=0) = ad!(s, generator, cos(2 * angle), sin(2 * angle), atol=atol)

function ad(s::PauliSentence, generators::AbstractVector{<:UPauli}, cosines::AbstractVector{<:Real}, sines::AbstractVector{<:Real}; atol::Real=0)
    length(generators) == length(cosines) == length(sines) || throw(DimensionMismatch("Generators and angles need to be equal size ($(length(generators)), $(length(cosines)), $(length(sines)))"))
    result = copy(s)
    length(generators) == 0 && return result
    for (generator, cosine, sine) in zip(reverse(generators), reverse(cosines), reverse(sines))
        result = ad(result, generator, cosine, sine, atol=atol)
    end
    return result
end
ad(s::PauliSentence, generators::AbstractVector{<:UPauli}, angles::AbstractVector{<:Real}; atol::Real=0) = ad(s, generators, cos.(2 .* angles), sin.(2 .* angles), atol=atol)

function ad!(s::PauliSentence, generators::AbstractVector{<:UPauli}, cosines::AbstractVector{<:Real}, sines::AbstractVector{<:Real}; atol::Real=0)
    length(generators) == length(cosines) == length(sines) || throw(DimensionMismatch("Generators and angles need to be equal size ($(length(generators)), $(length(cosines)), $(length(sines)))"))
    length(generators) == 0 && return s
    for (generator, cosine, sine) in zip(reverse(generators), reverse(cosines), reverse(sines))
        ad!(s, generator, cosine, sine, atol=atol)
    end
    return s
end
ad!(s::PauliSentence, generators::AbstractVector{<:UPauli}, angles::AbstractVector{<:Real}; atol::Real=0) = ad!(s, generators, cos.(2 .* angles), sin.(2 .* angles), atol=atol)

trace(s::PauliSentence) = haskey(s, 0) ? s[0] / 2^s.qubits : zero(eltype(s))