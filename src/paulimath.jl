const SIGNS = C8[1, 1im, -1, -1im]

function Base.:*(scalar::Number, p::AbstractPauli)
    typeof(p) <: Pauli && return SignedPauli(p, C8(scalar))
    return SignedPauli(p.string, p.sign * C8(scalar))
end

unsigned_prod(p::Pauli{<:Unsigned,Q}, q::Pauli{<:Unsigned,Q}) where {Q} = Pauli(p.string ⊻ q.string, Q)
Base.xor(p::Pauli, q::Pauli) = unsigned_prod(p, q)

function com(p::AbstractPauli{T,Q}, q::AbstractPauli{U,Q})::Pauli{<:Unsigned,Q} where {T,U,Q}
    S = promote_type(T, U)
    amask::S = 2^Q - 1
    bmask = ~amask
    a1 = p.string & amask
    b1 = (p.string & bmask) >> Q
    a2 = q.string & amask
    b2 = (q.string & bmask) >> Q
    (isodd(count_ones(a1 & b2)) ⊻ isodd(count_ones(b1 & a2))) && (return p ⊻ q)
    return Pauli{S,Q}(0)
end
function com(p::Unsigned, q::Unsigned, Q::Integer)
    amask::UInt = 2^Q - 1
    bmask = ~amask
    a1 = p & amask
    b1 = (p & bmask) >> Q
    a2 = q & amask
    b2 = (q & bmask) >> Q
    (isodd(count_ones(a1 & b2)) ⊻ isodd(count_ones(b1 & a2))) && (return p ⊻ q)
    return UInt(0)
end
Base.:^(p::AbstractPauli, q::AbstractPauli) = com(p, q)

function Base.:*(p::AbstractPauli{T,Q}, q::AbstractPauli{U,Q})::SignedPauli{<:Unsigned,Q} where {T,U,Q}
    S = promote_type(T, U)
    amask::S = 2^Q - 1
    bmask = ~amask
    a1 = p.string & amask
    b1 = (p.string & bmask) >> Q
    a2 = q.string & amask
    b2 = (q.string & bmask) >> Q

    overlap = count_ones((a1 | b1) & (a2 | b2))
    signstring = ((b1 & a1) & (~b2 & a2)) | ((b1 & ~a1) & (b2 & a2)) | ((~b1 & a1) & (b2 & ~a2))
    equalpaulis = (a1 & a2 & ~b1 & ~b2) | (a1 & b1 & a2 & b2) | (~a1 & ~a2 & b1 & b2)

    ind = (overlap - count_ones(equalpaulis)) % 4 + 1
    sign = SIGNS[ind] * SignedPauli(p).sign * SignedPauli(q).sign
    isodd(count_ones(signstring)) && (sign *= C8(-1))
    return SignedPauli(p.string ⊻ q.string, sign)
end
function _symplectic_prod(p::Unsigned, q::Unsigned, Q::Integer)
    amask::UInt = 2^Q - 1
    bmask = ~amask
    a1 = p & amask
    b1 = (p & bmask) >> Q
    a2 = q & amask
    b2 = (q & bmask) >> Q

    overlap = count_ones((a1 | b1) & (a2 | b2))
    signstring = ((b1 & a1) & (~b2 & a2)) | ((b1 & ~a1) & (b2 & a2)) | ((~b1 & a1) & (b2 & ~a2))
    equalpaulis = (a1 & a2 & ~b1 & ~b2) | (a1 & b1 & a2 & b2) | (~a1 & ~a2 & b1 & b2)

    ind = (overlap - count_ones(equalpaulis)) % 4 + 1
    sign = SIGNS[ind]
    isodd(count_ones(signstring)) && (sign *= C8(-1))
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

function ad(s::PauliSentence, generator::Pauli{<:Unsigned,Q}, angle::Real; atol::Real=0) where {Q}
    length(s) == 4^Q || throw(DimensionMismatch("PauliSentence must have length $(4^Q)"))
    result = copy(s)
    (iszero(angle) | iszero(generator.string)) && return result
    sentence = skipmissing(s)
    products = _symplectic_prod.(generator.string, UInt.(eachindex(sentence)), Q)
    noncomind = findall(p -> !isreal(p.second), products)
    noncomkeys = @view collect(eachindex(sentence))[noncomind]
    noncomprods = @view products[noncomind]
    vals = @view collect(sentence)[noncomind]
    result[noncomkeys] = cos(2 * angle) * vals
    result[first.(noncomprods)] = replace!(@view(result[first.(noncomprods)]), missing => 0.0) .+ sin(2 * angle) .* -imag(last.(noncomprods)) .* vals
    return replace(x -> (ismissing(x) || abs(x) <= atol ? missing : x), result)
end
function ad(s::PauliSentence, generator::Pauli{<:Unsigned,Q}, cosine::Real, sine::Real; atol::Real=0) where {Q}
    length(s) == 4^Q || throw(DimensionMismatch("PauliSentence must have length $(4^Q)"))
    result = copy(s)
    (iszero(sine) | iszero(generator.string)) && return result
    sentence = skipmissing(s)
    products = _symplectic_prod.(generator.string, UInt.(eachindex(sentence)), Q)
    noncomind = findall(p -> !isreal(p.second), products)
    noncomkeys = @view collect(eachindex(sentence))[noncomind]
    noncomprods = @view products[noncomind]
    vals = @view collect(sentence)[noncomind]
    result[noncomkeys] = cosine * vals
    result[first.(noncomprods)] = replace!(@view(result[first.(noncomprods)]), missing => 0.0) .+ sine .* -imag(last.(noncomprods)) .* vals
    return replace(x -> (ismissing(x) || abs(x) <= atol ? missing : x), result)
end
function ad(s::PauliSentence, generators::AbstractVector{<:Pauli}, angles::AbstractVector{<:Real}; atol::Real=0)
    length(generators) == length(angles) || throw(DimensionMismatch("Generators and angles need to be equal size ($(length(generators)) and $(length(angles)))"))
    result = copy(s)
    length(angles) == 0 && return result
    for (generator, angle) in zip(reverse(generators), reverse(angles))
        result = ad(result, generator, angle, atol=atol)
    end
    return result
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