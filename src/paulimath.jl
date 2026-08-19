# Powers of `i` and of `-1` appear in essentially every product, commutator and rotation.
# `im^k` / `(-1)^k` go through `Base.power_by_squaring` on `Complex{Int8}`, which costs a
# handful of complex multiplications; the phase only depends on `k mod 4` (resp. `k mod 2`),
# so a table lookup and a parity test give the identical value for a fraction of the work.
const _IPOWERS = (C8(1), C8(im), C8(-1), C8(-im))

"""
    _ipow(k)

``i^k`` as a `Complex{Int8}`, by table lookup on `k mod 4`.
"""
@inline _ipow(k::Integer) = @inbounds _IPOWERS[(k&3)+1]

"""
    _mpow(k)

``(-1)^k`` as a `Complex{Int8}`, by a parity test on `k`.
"""
@inline _mpow(k::Integer) = ifelse(isodd(k), C8(-1), C8(1))

@doc raw"""
    scalar * p

Attach a scalar to a Pauli string, giving a [`Pauli`](@ref).

The scalar multiplies the *Hermitian* operator, so `c * p` for a [`UPauli`](@ref) folds in
the ``i^{\#Y}`` convention phase on top of `c`; for a [`Pauli`](@ref) it multiplies the
phase the string already carries. Since a `Pauli` only holds a fourth root of unity,
anything else throws — use a one-term [`PauliSentence`](@ref) for a general coefficient.

# Examples
```jldoctest
julia> -1 * UPauli("Y-")
Pauli{UInt64}((-)Y-)

julia> UPauli("Y-") * im
Pauli{UInt64}((i)Y-)

julia> 0.5 * UPauli("Z-")
ERROR: ArgumentError: Sign must be 1, -1, im, or -im
```
"""
Base.:*(scalar::Number, p::UPauli) = Pauli(p, _ipow(county(p)) * scalar)
Base.:*(scalar::Number, p::Pauli) = Pauli(p.string, p.sign * scalar, p.qubits)
Base.:*(p::AbstractPauli, scalar::Number) = scalar * p

"""
    unsigned_prod(p, q)
    p ⊻ q

Product of two Pauli strings with the phase thrown away: the bit pattern alone.

This is the group operation of the symplectic representation — one `xor` — and is what
[`com`](@ref) and the rotations use when the phase is accounted for separately. `p * q`
gives the same string with its phase attached.

# Examples
```jldoctest
julia> tostring(UPauli("XX") ⊻ UPauli("ZZ"))
"YY"

julia> UPauli("XX") * UPauli("ZZ")     # the same string, with its phase
Pauli{UInt64}((-)YY)
```
"""
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

@doc raw"""
    p * q

Product of two Pauli strings, as a [`Pauli`](@ref) carrying the phase it picks up.

Every combination of [`UPauli`](@ref) and [`Pauli`](@ref) is accepted and the result is
always a `Pauli`, since even two bare strings generally multiply to a phased one. Both
factors must be strings of the same length and the same storage type.

The bits are `p ⊻ q`; the phase is ``i^{\#Y_p + \#Y_q}`` from the Hermiticity convention
times ``(-1)^{z_p \cdot x_q}`` from reordering the single-qubit factors.

# Examples
```jldoctest
julia> UPauli("X") * UPauli("Y")
Pauli{UInt64}((i)Z)

julia> UPauli("Y") * UPauli("X")       # anticommuting: the other sign
Pauli{UInt64}((-i)Z)

julia> p = UPauli("XY-Z"); tomatrix(p * p) == tomatrix(UPauli("----"))
true
```
"""
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
"""
    _symplectic_prod(p, q, Q)

Product of two bare symplectic strings as `string => sign`: the `⊻` of the bits paired with
the real ``±1`` the reordering costs.

The phase convention of [`PauliSentence`](@ref) is chosen exactly so that this real sign is
the whole phase of a product of bare strings, which is what lets the sentence arithmetic
stay in the coefficient type it started in.
"""
@inline function _symplectic_prod(p::Unsigned, q::Unsigned, Q::Integer)
    return p ⊻ q => _prodsign(p, q, Q)
end
"""
    _prodsign(p, q, Q)

Real ``±1`` phase picked up by the product of the two bare symplectic strings `p` and `q`,
i.e. the `second` field of [`_symplectic_prod`](@ref) without forming the product string.
"""
@inline _prodsign(p::Unsigned, q::Unsigned, Q::Integer) = _mpow(count_ones((p >> Q) & q))

@doc raw"""
    cis(p::UPauli, θ)
    cis(s::PauliSentence, θ)

The rotation ``e^{i\theta P}`` as a [`PauliSentence`](@ref).

A Pauli squares to the identity, so the exponential closes on two terms,

```math
e^{i\theta P} = \cos\theta \, \mathbb{1} + i \sin\theta \, P ,
```

with ``P`` the Hermitian Pauli of the string — exact, and built in constant time. For a
general sentence there is no such closed form: `cis(s, θ)` goes through the dense matrix
exponential and back, which costs ``2^Q``.

Conjugating an operator by such a rotation is [`ad`](@ref), which does not form the rotation
at all.

# Examples
```jldoctest
julia> u = cis(UPauli("Z-"), pi / 3);

julia> tomatrix(u) ≈ exp(im * (pi / 3) * tomatrix(Pauli("Z-")))
true

julia> tostring(cis(UPauli("Z-"), pi / 2))["Z-"]    # a quarter turn: only the Pauli is left
0.0 + 1.0im
```
"""
function Base.cis(p::UPauli{T,Q}, θ::Real) where {T,Q}
    # The identity is its own rotation axis: both terms land on the same key, and writing
    # them as two entries would silently keep only the second.
    iszero(p.string) && return PauliSentence{T,ComplexF64,Q}([T(0)], [cis(θ)])
    return PauliSentence{T,ComplexF64,Q}(
        [T(0), p.string],
        [cos(θ), _ipow(county(p) + 1) * sin(θ)],
    )
end

"""
    s + r
    s - r

Add or subtract two operators on the same number of qubits.

Terms are matched by their string and their coefficients combined; a term that appears in
only one of the two is carried over unchanged. Key and coefficient types promote, so adding
a `Float64` sentence to a `ComplexF64` one gives a `ComplexF64` one. Coefficients that
cancel are kept as explicit zeros — [`ad`](@ref) prunes, arithmetic does not — so `filter!`
them out if the sparsity matters.

# Examples
```jldoctest
julia> a = PauliSentence(["ZZ", "X-"], [1.0, 0.5]);

julia> b = PauliSentence(["X-"], [0.25]);

julia> tostring(a - b)["X-"]
0.25

julia> length(a - a)          # the cancelled terms are still there, as zeros
2

julia> length(filter!(p -> !iszero(p.second), a - a))
0
```
"""
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

"""
    c * s
    s * c

Scale every coefficient of an operator by the number `c`.

Unlike scaling a bare [`Pauli`](@ref), any number will do: the coefficient type promotes to
hold it.

# Examples
```jldoctest
julia> s = PauliSentence(["ZZ"], [1.0]);

julia> tostring(2im * s)["ZZ"]
0.0 + 2.0im
```
"""
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

@doc raw"""
    s * r

Operator product, term by term.

Each pair of strings contributes one term, `key1 ⊻ key2` with the ``\pm 1`` its reordering
costs, and terms landing on the same string are summed — so the cost is the product of the
two lengths, and the result is generally shorter than that. Coefficients come back as
`ComplexF64`.

The product is not commutative; `s * r - r * s` is the commutator, and vanishes term by term
exactly when every pair of strings [`com`](@ref)mutes.

# Examples
```jldoctest
julia> x, z = PauliSentence(["X"], [1.0]), PauliSentence(["Z"], [1.0]);

julia> tostring(x * z)["Y"]        # XZ = -iY
0.0 - 1.0im

julia> tostring(x * z - z * x)["Y"]   # the commutator [X, Z] = -2iY
0.0 - 2.0im

julia> length(x * x), tostring(x * x)["-"]
(1, 1.0 + 0.0im)
```
"""
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

"""
    s^n

Repeated operator product, for a non-negative integer `n`; `s^0` is the identity.

# Examples
```jldoctest
julia> s = PauliSentence(["Z-", "-X"], [1.0, 2.0]);

julia> tomatrix(s^3) ≈ tomatrix(s)^3
true

julia> tostring(s^0)
Dict{String, Float64} with 1 entry:
  "--" => 1.0
```
"""
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
"""
    _prune!(s, atol)

Drop the terms of `s` whose coefficient has magnitude `atol` or less, in place.

`atol = 0` drops the exact zeros only, which is the common case after a rotation and avoids
a `hypot` call per term.
"""
_prune!(s::PauliSentence, atol::Real) =
    if iszero(atol)
        filter!(p -> !iszero(p.second), s)
    else
        filter!(p -> abs(p.second) > atol, s)
    end

# Whether `_prune!` would drop this coefficient. Checking it as each coefficient is written
# costs a comparison already in registers and saves sweeping the whole sentence afterwards,
# which a rotation rarely has anything to show for.
"""
    _prunable(v, atol)

Whether [`_prune!`](@ref) would drop a coefficient of value `v`.
"""
@inline _prunable(v::Number, atol::Real) = iszero(atol) ? iszero(v) : !(abs(v) > atol)

@doc raw"""
    ad(s::PauliSentence, generator::UPauli, angle; atol=0)
    ad(s::PauliSentence, generators::PauliList, angles; atol=0)
    ad(s::PauliSentence, generator, cosine, sine; atol=0)

Conjugate `s` by the rotation generated by `generator`:

```math
s \;\longmapsto\; U s U^\dagger, \qquad U = e^{i\theta A},
```

with `A` the Hermitian Pauli of the generator. Given a `PauliList`, the rotations are
applied from the **last** generator to the first, so undoing a sequence means reversing both
the generators and the angles.

Only the terms that *anticommute* with the generator move, and each one mixes with exactly
one partner, `generator ⊻ key`, so the cost is one pass over the sentence and no matrix
anywhere. The rotation is exact: nothing is truncated except by `atol`.

Coefficients whose magnitude drops to `atol` or below are dropped, which is what keeps
long rotation sequences from filling up with numerical dust. [`ad!`](@ref) does the same in
place. The four-argument form takes ``\cos 2\theta`` and ``\sin 2\theta`` directly, for a
sweep that reuses one pair of angles across many sentences.

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
    ad!(s::PauliSentence, generator, cosine, sine; atol=0)

In-place [`ad`](@ref): rotate `s` and return it, without allocating a new sentence.

`s` must already hold `ComplexF64` coefficients, since a rotation generally produces complex
ones. Rotating by a `PauliList` of generators applies them last to first, as `ad` does, and
does not copy in between — which is what makes it the right call in an optimization loop.

# Examples
```jldoctest
julia> s = PauliSentence{UInt,ComplexF64}(["Z-"], [1.0]);

julia> ad!(s, UPauli("Y-"), pi / 4, atol = 1e-12);

julia> tostring(s)
Dict{String, ComplexF64} with 1 entry:
  "X-" => -1.0+0.0im
```
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

@doc raw"""
    trace(s::PauliSentence)

Trace of the operator, ``\mathrm{Tr}\,\mathcal{S}``.

Every Pauli string but the identity is traceless, so the whole trace sits in the identity
coefficient: `trace(s) == 2^Q * s[0]`, and zero when `s` has no identity term. It agrees
with `tr(tomatrix(s))` by construction, without building the matrix.

# Examples
```jldoctest
julia> s = PauliSentence(["--", "ZZ"], [0.5, 1.0]);

julia> trace(s)
2.0

julia> trace(PauliSentence(["ZZ"], [1.0]))
0.0
```
"""
trace(s::PauliSentence) = haskey(s, 0) ? 2^s.qubits * s[0] : zero(valtype(s))
