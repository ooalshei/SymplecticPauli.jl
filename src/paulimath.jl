# const SIGNS = C8[1, 1im, -1, -1im]

function Base.:*(scalar::Number, p::AbstractPauli)
    typeof(p) <: UPauli && return Pauli(p, C8(scalar))
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
    return UPauli{S,Q}(0)
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
    a = result & amask
    b = (result & bmask) >> Q
    sign::C8 = (-1im)^count_ones(a & b) * (1im)^count_ones(a1 & b1) * SignedPauli(p).sign * (1im)^count_ones(a2 & b2) * SignedPauli(q).sign * (-1)^count_ones(b1 & a2)
    return SignedPauli(result, sign, Q)
end
function _symplectic_prod(p::T, q::T, Q::Integer) where {T<:Unsigned}
    amask::T = 2^Q - 1
    bmask = ~amask
    a1 = p & amask
    b1 = (p & bmask) >> Q
    a2 = q & amask
    b2 = (q & bmask) >> Q

    # overlap = count_ones((a1 | b1) & (a2 | b2))
    # signstring = ((b1 & a1) & (~b2 & a2)) | ((b1 & ~a1) & (b2 & a2)) | ((~b1 & a1) & (b2 & ~a2))
    # equalUPaulis = (a1 & a2 & ~b1 & ~b2) | (a1 & b1 & a2 & b2) | (~a1 & ~a2 & b1 & b2)
    # ind = (overlap - count_ones(equalUPaulis)) % 4 + 1
    # sign = SIGNS[ind]
    # isodd(count_ones(signstring)) && (sign *= C8(-1))
    result = p ⊻ q
    a = result & amask
    b = (result & bmask) >> Q
    sign::C8 = (-1im)^count_ones(a & b) * (1im)^count_ones(a1 & b1) * (1im)^count_ones(a2 & b2) * (-1)^count_ones(b1 & a2)
    return p ⊻ q => sign
end

function Base.:*(A::PauliSentence, B::PauliSentence)
    length(A) == length(B) || throw(DimensionMismatch("PauliSentences must have the same length."))
    strings = skipmissing(B)
    result = PauliSentence(zeros(ComplexF64, length(A)))
    Q = Int(log2(length(A)) ÷ 2)
    for (key, value) in pairs(skipmissing(A))
        products = _symplectic_prod.(UInt(key), UInt.(keys(strings)), Q)
        result[first.(products)] += value .* last.(products) .* values(strings)
    end
    return replace(result, zero(ComplexF64) => missing)
end

function ad(s::PauliSentence, generator::Pauli{T,Q}, cosine::Real, sine::Real; atol::Real=0) where {T,Q}
    length(s) == 4^Q || throw(DimensionMismatch("PauliSentence must have length $(4^Q)"))
    result = copy(s)
    (iszero(sine) | iszero(generator.string)) && return result
    sentence = skipmissing(s)
    products = _symplectic_prod.(generator.string, T.(eachindex(sentence)), Q)
    noncomind = findall(p -> !isreal(p.second), products)
    noncomkeys = @view collect(eachindex(sentence))[noncomind]
    noncomprods = @view products[noncomind]
    vals = @view collect(sentence)[noncomind]
    result[noncomkeys] = cosine * vals
    result[first.(noncomprods)] = replace!(@view(result[first.(noncomprods)]), missing => 0) .+ sine .* -imag(last.(noncomprods)) .* vals
    return replace(x -> (ismissing(x) || abs(x) <= atol ? missing : x), result)
end
function ad(s::PauliSentence, generators::AbstractVector{<:Pauli}, cosines::AbstractVector{<:Real}, sines::AbstractVector{<:Real}; atol::Real=0)
    length(generators) == length(cosines) == length(sines) || throw(DimensionMismatch("Generators and angles need to be equal size ($(length(generators)), $(length(cosines)), $(length(sines)))"))
    result = copy(s)
    length(generators) == 0 && return result
    for (generator, cosine, sine) in zip(reverse(generators), reverse(cosines), reverse(sines))
        result = ad(result, generator, cosine, sine, atol=atol)
    end
    return result
end
ad(s::PauliSentence, generator::Pauli, angle::Real; atol::Real=0) = ad(s, generator, cos(2 * angle), sin(2 * angle), atol=atol)
ad(s::PauliSentence, generators::AbstractVector{<:Pauli}, angles::AbstractVector{<:Real}; atol::Real=0) = ad(s, generators, cos.(2 .* angles), sin.(2 .* angles), atol=atol)
