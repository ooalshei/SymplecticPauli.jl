```@meta
CurrentModule = SymplecticPauli
```

# Operators

A [`PauliSentence`](@ref) is an operator: Pauli strings with coefficients,

```math
\mathcal{S} = \sum_j c_j P_j ,
```

stored as an `AbstractDict` from the bit pattern of $P_j$ to $c_j$.

```@example operators
using SymplecticPauli

H = PauliSentence(["ZZ-", "-ZZ", "X--", "-X-", "--X"], [1.0, 1.0, 0.5, 0.5, 0.5])
```

The keys are strings on three qubits and the values their couplings — this is the
transverse-field Ising Hamiltonian at $J = 1$, $g = 0.5$. The same list of strings can come
from a [`PauliList`](@ref), from `UPauli`s, from index vectors or from a matrix of indices,
one string per column:

```@example operators
PauliSentence(PauliList(["ZZ", "X-"]), [1.0, -0.5])
```

## The phase convention

This is the one thing worth reading carefully. A coefficient multiplies the **bare** embedding
of its string — `X` ↦ $σ_1$, `Z` ↦ $σ_3$, and a `Y` position ↦ the *real* antisymmetric
[`σ₂real`](@ref) rather than $σ_2$. That is what makes a product of two strings an `xor` and a
$\pm 1$, with no $i$ anywhere.

Constructing a sentence from anything that *names* its strings folds in the $i^{\#Y}$ phase
that turns each bare string into the Hermitian Pauli operator, so a real coefficient vector
gives a Hermitian operator:

```@example operators
s = PauliSentence(["Y-"], [1.0])
```

```@example operators
tomatrix(s) ≈ tomatrix(s)'
```

The stored coefficient is $i$, and [`tostring`](@ref) divides that phase back out, so what
you read is what you entered:

```@example operators
tostring(s)
```

Only the raw-bits constructor takes coefficients exactly as given, without a convention:

```@example operators
PauliSentence(UInt[3], [1.0], 1)[UInt(3)]     # the bare Y, coefficient untouched
```

One consequence to keep in mind: a string with an odd number of `Y`s forces the coefficient
type to widen, since the phase is imaginary. `PauliSentence(["Y-"], [1.0])` stores
`ComplexF64`, while `PauliSentence(["YY"], [1.0])` stays `Float64`.

## As a dictionary

Everything `Base` offers on a dictionary is available:

```@example operators
haskey(H, UPauli("ZZ-").string), H[UPauli("X--").string], get(H, UInt(0), 0.0)
```

```@example operators
length(filter(p -> p.second > 0.9, H))
```

```@example operators
sum(H.qubits - counti(k, H.qubits) for k in keys(H))   # total weight of the Hamiltonian
```

`copy`, `empty`, `delete!`, `pop!`, `filter!`, `sizehint!` and iteration over `key => value`
pairs all behave as expected; the ones that return a container return a `PauliSentence`, not
a `Dict`.

## Arithmetic

Sentences add, subtract, scale, multiply and take non-negative integer powers. Types promote
rather than being pinned to whatever came first:

```@example operators
A = PauliSentence(["X-"], [1.0])
B = PauliSentence(["Z-"], [2.0])

tostring(A + B)
```

```@example operators
tostring(A * B)          # XZ = -iY
```

```@example operators
tostring(A * B - B * A)  # the commutator [X, Z] = -2iY
```

The product costs one term per pair of strings, and terms landing on the same string are
summed — so a product of two sparse operators stays sparse if their strings collide, and
grows if they do not. Cancelling terms are kept as explicit zeros; arithmetic does not prune,
so `filter!` them out when the sparsity matters.

## Trace

Every Pauli string but the identity is traceless, so the trace of a sentence is the identity
coefficient scaled by the dimension — no matrix needed:

```@example operators
T = PauliSentence(["--", "ZZ"], [0.5, 1.0])
trace(T), trace(H)
```

## Dense matrices

[`tomatrix`](@ref) builds the $2^Q \times 2^Q$ matrix, and the [`PauliSentence`](@ref)
constructor takes one apart again:

```@example operators
using LinearAlgebra

m = tomatrix(PauliSentence(["ZZ", "X-"], [1.0, -0.5]))
```

```@example operators
tostring(PauliSentence(m))
```

Both directions are exponential in the number of qubits. They are there to check small cases
against — an identity that holds for every string on three qubits is usually an identity —
not to compute with.
