```@meta
CurrentModule = SymplecticPauli
DocTestSetup = :(using SymplecticPauli)
```

# Pauli strings

## Building one

A single Pauli string is a [`UPauli`](@ref) — a bit pattern with no phase — or a
[`Pauli`](@ref), which carries one of $1, -1, i, -i$ alongside it. Either can be written as
text, as indices, or as the raw bits:

```@example strings
using SymplecticPauli

UPauli("XY-Z")           # 'X', 'Y', 'Z' and '-' (or 'I'), leftmost qubit first
```

```@example strings
UPauli([2, 3, 1, 4])     # 1-4 for I, X, Y, Z
```

```@example strings
UPauli(UInt(0xa3), 4)    # the bits, plus the number of qubits they describe
```

The three are the same string. Reading it back is [`tostring`](@ref):

```@example strings
tostring(UPauli("XY-Z"))
```

The leftmost character is qubit 1, and it is the *lowest* bit of each half — so `"X-"` is
`0b0001` and `"-X"` is `0b0010`, while `"Z-"` is `0b0100`:

```@example strings
[tostring(UPauli(UInt(b), 2)) => string(b, base = 2, pad = 4) for b in 0x0:0xf]
```

## Storage

The integer the bits live in is a type parameter, and so is the number of qubits: a
`UPauli{UInt16,4}` is a four-qubit string held in a `UInt16`. Since a string needs two bits
per qubit, `UInt` covers 32 qubits and `UInt128` covers 64. Ask for more than the type can
hold and you get an error rather than silent truncation:

```@example strings
UPauli{UInt16}("XY-Z")
```

```jldoctest
julia> UPauli{UInt8}("XXXXX")
ERROR: ArgumentError: String length cannot exceed 4. Consider using a larger unsigned integer type.
```

Both types are `isbitstype`, so a string costs no allocation and an array of them is a flat
block of memory.

## The phase

`UPauli` deliberately has no phase: it is a bit pattern, and two of them multiply by `⊻`.
The phase appears the moment you ask what operator the pattern *is*. The convention is the
one that makes each string Hermitian, $i^{\#Y}$:

```@example strings
Pauli(UPauli("Y")).sign, Pauli(UPauli("YY")).sign, Pauli(UPauli("XZ")).sign
```

so that [`tomatrix`](@ref) of any string squares to the identity and equals its own adjoint:

```@example strings
m = tomatrix(Pauli("Y"))
```

A `Pauli` tracks its phase exactly through products, and prints it in front of the string.
Note that the printed phase is *relative* to the Hermitian convention, so a plain `"Y"` shows
as `(+)`:

```@example strings
Pauli("Y"), Pauli("Y", -1), UPauli("X") * UPauli("Y")
```

Multiplying by a scalar gives a `Pauli` too, which means the scalar has to be a fourth root
of unity — for anything else, use a one-term [`PauliSentence`](@ref):

```@example strings
-1 * UPauli("Y-")
```

## Products and commutation

Every combination of `UPauli` and `Pauli` multiplies, and the result is a `Pauli`:

```@example strings
UPauli("XX") * UPauli("ZZ")
```

Dropping the phase — the bit pattern alone — is `⊻`:

```@example strings
tostring(UPauli("XX") ⊻ UPauli("ZZ"))
```

[`com`](@ref) answers the question that actually comes up in an algebra: *do these two
commute*, and *what is their product*? Both fall out of the same handful of bit operations,
so asking for one gives the other for free:

```@example strings
com(UPauli("XX"), UPauli("ZZ")), com(UPauli("X-"), UPauli("Z-"))
```

A commutator is then two lines with no matrices in sight: $[P, Q] = 0$ when they commute,
and $2PQ$ when they do not.

## Counting

[`countx`](@ref), [`county`](@ref), [`countz`](@ref) and [`counti`](@ref) are popcounts on
the two halves. A `Y` sets both bits of its qubit, so it is counted by `countx` and `countz`
as well as by `county`:

```@example strings
p = UPauli("XY-Z")
countx(p), county(p), countz(p), counti(p)
```

`counti` is `Q` minus the other three, which subtracts a `Y` three times over. It is a weight
ordering — larger means closer to the identity — rather than a literal count of identity
factors, and for a string with no `Y`s the two coincide.

## Lists of strings

A [`PauliList`](@ref) is an ordered list of strings on the same number of qubits: a
generating set, an algebra, a subalgebra. It is an `AbstractVector` of the *raw* bit
patterns, with the qubit count stored once for the whole list:

```@example strings
v = PauliList(["XX-", "-YZ", "Z-Z"])
```

```@example strings
length(v), v.qubits, eltype(v)
```

The whole `Base` vector vocabulary works, and slices come back as `PauliList`s:

```@example strings
push!(v, UPauli("YYY").string)
tostring(v[2:end])
```

Elements are bare integers, so wrap one to treat it as a string on its own:

```@example strings
UPauli(v[1], v.qubits)
```

That is the trade the type makes: a million strings is a million integers, not a million
objects, and the qubit count is not repeated a million times.
