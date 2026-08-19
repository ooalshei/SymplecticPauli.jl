# Powers of `i` and of `-1` appear in essentially every product, commutator and rotation.
# `im^k` / `(-1)^k` go through `Base.power_by_squaring` on `Complex{Int8}`, which costs a
# handful of complex multiplications; the phase only depends on `k mod 4` (resp. `k mod 2`),
# so a table lookup and a parity test give the identical value for a fraction of the work.
const _IPOWERS = (C8(1), C8(im), C8(-1), C8(-im))
@inline _ipow(k::Integer) = @inbounds _IPOWERS[(k&3)+1]
@inline _mpow(k::Integer) = ifelse(isodd(k), C8(-1), C8(1))

Base.:*(scalar::Number, p::UPauli) = Pauli(p, _ipow(county(p)) * scalar)
Base.:*(scalar::Number, p::Pauli) = Pauli(p.string, p.sign * scalar, p.qubits)
Base.:*(p::AbstractPauli, scalar::Number) = scalar * p

unsigned_prod(p::UPauli{<:Unsigned,Q}, q::UPauli{<:Unsigned,Q}) where {Q} =
    UPauli(p.string ⊻ q.string, Q)
Base.xor(p::UPauli, q::UPauli) = unsigned_prod(p, q)

@doc raw"""
    com(p, q)
    com(p::Unsigned, q::Unsigned, Q)

Product string of two Paulis paired with whether they commute: `p ⊻ q => commutes`.

Both come out of the same handful of bit operations, so code that needs the product only for
the anticommuting pairs — commutators, rotations, algebra closures — gets the test for free.

# Examples
```jldoctest
julia> com(UPauli("XX"), UPauli("ZZ")).second      # two anticommutations cancel
true

julia> com(UPauli("X-"), UPauli("Z-")).second
false

julia> tostring(com(UPauli("X-"), UPauli("Z-")).first)
"Y-"
```
"""
function com(
    p::AbstractPauli{<:Unsigned,Q},
    q::AbstractPauli{<:Unsigned,Q},
)::Pair{UPauli{<:Unsigned,Q},Bool} where {Q}
    b1 = p.string >> Q
    b2 = q.string >> Q
    return p ⊻ q => iseven(count_ones((p.string & b2) ⊻ (b1 & q.string)))
end
function com(p::Unsigned, q::Unsigned, Q::Integer)
    b1 = p >> Q
    b2 = q >> Q
    return p ⊻ q => iseven(count_ones((p & b2) ⊻ (b1 & q)))
end
Base.:^(p::AbstractPauli, q::AbstractPauli) = com(p, q)

function Base.:*(p::UPauli{<:Unsigned,Q}, q::UPauli{<:Unsigned,Q}) where {Q}
    b1 = p.string >> Q
    result = p.string ⊻ q.string
    sign::C8 = _ipow(county(p) + county(q)) * _mpow(count_ones(b1 & q.string))
    return Pauli(result, sign, Q)
end
function Base.:*(p::Pauli{<:Unsigned,Q}, q::UPauli{<:Unsigned,Q}) where {Q}
    b1 = p.string >> Q
    result = p.string ⊻ q.string
    sign::C8 = _ipow(county(q)) * p.sign * _mpow(count_ones(b1 & q.string))
    return Pauli(result, sign, Q)
end
function Base.:*(p::UPauli{<:Unsigned,Q}, q::Pauli{<:Unsigned,Q}) where {Q}
    b1 = p.string >> Q
    result = p.string ⊻ q.string
    sign::C8 = _ipow(county(p)) * q.sign * _mpow(count_ones(b1 & q.string))
    return Pauli(result, sign, Q)
end

function Base.:*(p::Pauli{<:Unsigned,Q}, q::Pauli{<:Unsigned,Q}) where {Q}
    b1 = p.string >> Q
    result = p.string ⊻ q.string
    sign::C8 = p.sign * q.sign * _mpow(count_ones(b1 & q.string))
    return Pauli(result, sign, Q)
end
@inline function _symplectic_prod(p::Unsigned, q::Unsigned, Q::Integer)
    return p ⊻ q => _prodsign(p, q, Q)
end
"""
    _prodsign(p, q, Q)

Real ``±1`` phase picked up by the product of the two bare symplectic strings `p` and `q`,
i.e. the `second` field of [`_symplectic_prod`](@ref) without forming the product string.
"""
@inline _prodsign(p::Unsigned, q::Unsigned, Q::Integer) = _mpow(count_ones((p >> Q) & q))

function Base.cis(p::UPauli{T,Q}, θ::Real) where {T,Q}
    # The identity is its own rotation axis: both terms land on the same key, and writing
    # them as two entries would silently keep only the second.
    iszero(p.string) && return PauliSentence{T,ComplexF64,Q}([T(0)], [cis(θ)])
    return PauliSentence{T,ComplexF64,Q}(
        [T(0), p.string],
        [cos(θ), _ipow(county(p) + 1) * sin(θ)],
    )
end

function Base.:+(s::PauliSentence{Ts,Ns,Q}, r::PauliSentence{Tr,Nr,Q}) where {Ts,Ns,Tr,Nr,Q}
    T = promote_type(Ts, Tr)
    N = promote_type(Ns, Nr)
    result = Dict{T,N}()
    sizehint!(result, length(s) + length(r))
    for (key, value) in s.sentence
        result[key] = value
    end
    for (key, value) in r.sentence
        result[key] = get(result, key, zero(N)) + value
    end
    return PauliSentence{T,N,Q}(result, iscopy=false, check=false)
end
Base.:+(
    s::PauliSentence{<:Unsigned,<:Number,Q},
    r::PauliSentence{<:Unsigned,<:Number,Q},
    t::PauliSentence{<:Unsigned,<:Number,Q},
    u::PauliSentence{<:Unsigned,<:Number,Q}...,
) where {Q} = +(+(s, r), t, u...)
function Base.:-(s::PauliSentence{Ts,Ns,Q}, r::PauliSentence{Tr,Nr,Q}) where {Ts,Ns,Tr,Nr,Q}
    T = promote_type(Ts, Tr)
    N = promote_type(Ns, Nr)
    result = Dict{T,N}()
    sizehint!(result, length(s) + length(r))
    for (key, value) in s.sentence
        result[key] = value
    end
    for (key, value) in r.sentence
        result[key] = get(result, key, zero(N)) - value
    end
    return PauliSentence{T,N,Q}(result, iscopy=false, check=false)
end

function Base.:*(c::N, s::PauliSentence{T,Ns,Q}) where {T,N,Ns,Q}
    V = promote_type(N, Ns)
    result = Dict{T,V}()
    sizehint!(result, length(s))
    for (key, value) in s.sentence
        result[key] = c * value
    end
    return PauliSentence{T,V,Q}(result, iscopy=false, check=false)
end
Base.:*(s::PauliSentence, c::Number) = c * s

function Base.:*(
    s::PauliSentence{Ts,<:Number,Q},
    r::PauliSentence{Tr,<:Number,Q},
) where {Ts,Tr,Q}
    T = promote_type(Ts, Tr)
    result = Dict{T,ComplexF64}()
    sizehint!(result, length(s) * length(r))
    for (key1, value1) in s.sentence
        for (key2, value2) in r.sentence
            key = T(key1 ⊻ key2)
            result[key] =
                get(result, key, zero(ComplexF64)) +
                _prodsign(key1, key2, Q) * value1 * value2
        end
    end
    return PauliSentence{T,ComplexF64,Q}(result, iscopy=false, check=false)
end

Base.:^(s::PauliSentence{T,N,Q}, x::Integer) where {T,N,Q} =
    if x >= 0
        prod(s for _ in 1:x; init=PauliSentence{T,N,Q}(T[0], [1]))
    else
        throw(ArgumentError("Exponent must be a non-negative integer."))
    end

Base.cis(s::PauliSentence, θ::Real) =
    PauliSentence{keytype(s),ComplexF64}(cis(θ * tomatrix(s)))

# Dropping the coefficients that a rotation has sent (back) to zero is the common case, and
# `abs` on a complex number is a `hypot` call; comparing against an explicit tolerance is
# only needed when one was asked for.
_prune!(s::PauliSentence, atol::Real) = if iszero(atol)
    filter!(p -> !iszero(p.second), s)
else
    filter!(p -> abs(p.second) > atol, s)
end

# Whether `_prune!` would drop this coefficient. Checking it as each coefficient is written
# costs a comparison already in registers and saves sweeping the whole sentence afterwards,
# which a rotation rarely has anything to show for.
@inline _prunable(v::Number, atol::Real) = iszero(atol) ? iszero(v) : !(abs(v) > atol)

@doc raw"""
    ad(s::PauliSentence, generator::UPauli, angle; atol=0)
    ad(s::PauliSentence, generators::PauliList, angles; atol=0)

Conjugate `s` by the rotation generated by `generator`:

```math
s \;\longmapsto\; U s U^\dagger, \qquad U = e^{i\theta A},
```

with `A` the Hermitian Pauli of the generator. Given a `PauliList`, the rotations are applied
from the **last** generator to the first, so undoing a sequence means reversing both the
generators and the angles.

Coefficients whose magnitude drops to `atol` or below are dropped, which is what keeps
long rotation sequences from filling up with numerical dust. [`ad!`](@ref) does the same in
place.

# Examples
```jldoctest
julia> s = PauliSentence(PauliList(["Z-"]), [1.0]);

julia> tostring(ad(s, UPauli("Y-"), pi / 4, atol=1e-12))   # a quarter turn: Z ↦ -X
Dict{String, ComplexF64} with 1 entry:
  "X-" => -1.0+0.0im

julia> gens = PauliList(["Y-", "-Y"]); angles = [0.3, -0.7];

julia> back = ad(ad(s, gens, angles), reverse(gens), -reverse(angles));

julia> tostring(back)
Dict{String, ComplexF64} with 1 entry:
  "Z-" => 1.0+0.0im
```
"""
function ad(
    s::PauliSentence{Ts,<:Number,Q},
    generator::UPauli{T,Q},
    cosine::Real,
    sine::Real;
    atol::Real=0,
) where {Ts,T,Q}
    result = PauliSentence{promote_type(T, Ts),ComplexF64,Q}(s)
    (iszero(sine) | iszero(generator.string)) && return result
    modsine = _ipow(county(generator) + 1) * sine
    g = generator.string
    dict = result.sentence
    prunable = false
    for (key, value) in s.sentence
        # `com` returns the product string together with the commutation bit, so the
        # anticommutation test and the target key come out of a single pass over the bits.
        product = com(g, key, Q)
        product.second && continue
        newkey = keytype(dict)(product.first)
        kept = dict[key] + (cosine - 1) * value
        dict[key] = kept
        moved = get(dict, newkey, zero(ComplexF64)) + modsine * _prodsign(g, key, Q) * value
        dict[newkey] = moved
        prunable |= _prunable(kept, atol) | _prunable(moved, atol)
    end
    return prunable ? _prune!(result, atol) : result
end
ad(s::PauliSentence, generator::UPauli, angle::Real; atol::Real=0) =
    ad(s, generator, cos(2 * angle), sin(2 * angle), atol=atol)

"""
    ad!(s::PauliSentence, generator, angle; atol=0)
    ad!(s::PauliSentence, generators::PauliList, angles; atol=0)

In-place [`ad`](@ref): rotate `s` and return it, without allocating a new sentence.

`s` must already hold `ComplexF64` coefficients, since a rotation generally produces complex
ones.
"""
function ad!(
    s::PauliSentence{<:Unsigned,ComplexF64,Q},
    generator::UPauli{<:Unsigned,Q},
    cosine::Real,
    sine::Real;
    atol::Real=0,
) where {Q}
    (iszero(sine) | iszero(generator.string)) && return s
    modsine = _ipow(county(generator) + 1) * sine
    g = generator.string
    dict = s.sentence
    # The rotation mixes each term with its partner `g ⊻ key`, so the update has to read the
    # coefficients as they were before this call: one snapshot of the pairs, taken up front.
    snapshot = collect(dict)
    prunable = false
    for (key, value) in snapshot
        product = com(g, key, Q)
        product.second && continue
        newkey = keytype(dict)(product.first)
        kept = dict[key] + (cosine - 1) * value
        dict[key] = kept
        moved = get(dict, newkey, zero(ComplexF64)) + modsine * _prodsign(g, key, Q) * value
        dict[newkey] = moved
        prunable |= _prunable(kept, atol) | _prunable(moved, atol)
    end
    return prunable ? _prune!(s, atol) : s
end
ad!(s::PauliSentence, generator::UPauli, angle::Real; atol::Real=0) =
    ad!(s, generator, cos(2 * angle), sin(2 * angle), atol=atol)

function ad(
    s::PauliSentence{Ts,<:Number,Q},
    generators::PauliList{Tg},
    cosines::AbstractVector{<:Real},
    sines::AbstractVector{<:Real};
    atol::Real=0,
) where {Ts,Tg,Q}
    length(generators) == length(cosines) == length(sines) || throw(
        DimensionMismatch(
            "Generators and angles must have equal lengths " *
            "($(length(generators)), $(length(cosines)), $(length(sines))).",
        ),
    )
    # One conversion up front, then rotate in place: chaining the out-of-place `ad` would
    # copy the whole sentence once per generator.
    result = PauliSentence{promote_type(Ts, Tg),ComplexF64,Q}(s)
    for i in reverse(eachindex(generators))
        ad!(result, UPauli(generators[i], Q), cosines[i], sines[i], atol=atol)
    end
    return result
end
ad(s::PauliSentence, generators::PauliList, angles::AbstractVector{<:Real}; atol::Real=0) =
    ad(s, generators, cos.(2 .* angles), sin.(2 .* angles), atol=atol)

function ad!(
    s::PauliSentence{<:Unsigned,<:Number,Q},
    generators::PauliList{T},
    cosines::AbstractVector{<:Real},
    sines::AbstractVector{<:Real};
    atol::Real=0,
) where {T,Q}
    length(generators) == length(cosines) == length(sines) || throw(
        DimensionMismatch(
            "Generators and angles must have equal lengths " *
            "($(length(generators)), $(length(cosines)), $(length(sines))).",
        ),
    )
    for i in reverse(eachindex(generators))
        ad!(s, UPauli{T,Q}(generators[i]), cosines[i], sines[i], atol=atol)
    end
    return s
end
ad!(s::PauliSentence, generators::PauliList, angles::AbstractVector{<:Real}; atol::Real=0) =
    ad!(s, generators, cos.(2 .* angles), sin.(2 .* angles), atol=atol)

trace(s::PauliSentence) = haskey(s, 0) ? 2^s.qubits * s[0] : zero(valtype(s))
