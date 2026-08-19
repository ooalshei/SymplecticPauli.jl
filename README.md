# SymplecticPauli.jl

[![Build Status](https://github.com/ooalshei/SymplecticPauli.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ooalshei/SymplecticPauli.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ooalshei/SymplecticPauli.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ooalshei/SymplecticPauli.jl)
[![Coverage](https://coveralls.io/repos/github/ooalshei/SymplecticPauli.jl/badge.svg?branch=main)](https://coveralls.io/github/ooalshei/SymplecticPauli.jl?branch=main)

Pauli operators as **symplectic bit strings** — products, commutators and rotations at the
cost of a few machine instructions.

A Pauli string on `Q` qubits is one unsigned integer of `2Q` bits. The low `Q` bits say which
qubits carry an `X` factor, the high `Q` bits which carry a `Z`, and a qubit with both bits
set carries a `Y`:

```
P = P₁ ⊗ ⋯ ⊗ P_Q   ↔   (z | x),    (zⱼ, xⱼ) = (0,0), (0,1), (1,0), (1,1)  for  I, X, Z, Y
```

Everything else follows. The product of two strings is an exclusive or,

```
PQ ∝ (z_P ⊻ z_Q | x_P ⊻ x_Q)
```

its phase is `(-1)^(z_P · x_Q)` times a power of `i` counting the `Y`s, and the two strings
commute exactly when the symplectic form `z_P · x_Q + x_P · z_Q` is even. A product is one
`xor`; a commutation test is two `and`s, an `xor` and a popcount. No matrix is ever built —
`tomatrix` exists for checking small cases, and that is all it is for.

## Installation

The package is not registered:

```julia
using Pkg
Pkg.add(url="https://github.com/ooalshei/SymplecticPauli.jl")
```

## Example

```julia
using SymplecticPauli

p, q = UPauli("XY-Z"), UPauli("Z--Z")

p * q             # Pauli{UInt64}((-i)YY--)   the product, and the phase it picks up
com(p, q).second  # false                     …and whether they commute
```

An operator is a `PauliSentence` — strings with coefficients, behaving as a dictionary from
one to the other:

```julia
# the transverse-field Ising Hamiltonian on three qubits, J = 1, g = 0.5
H = PauliSentence(["ZZ-", "-ZZ", "X--", "-X-", "--X"], [1.0, 1.0, 0.5, 0.5, 0.5])

length(H), H.qubits    # (5, 3)
trace(H)               # 0.0
```

Conjugating it by `exp(iθA)` is `ad`, which touches only the terms that anticommute with the
generator — one pass over the sentence, no matrix exponential:

```julia
rotated = ad(H, UPauli("Y--"), pi/4, atol=1e-12)

tostring(rotated)
# Dict{String, ComplexF64} with 5 entries:
#   "-ZZ" => 1.0+0.0im
#   "--X" => 0.5+0.0im
#   "Z--" => 0.5+0.0im      X-- rotated onto Z--
#   "-X-" => 0.5+0.0im
#   "XZ-" => -1.0+0.0im     ZZ- rotated onto XZ-
```

Given a `PauliList` of generators and a vector of angles, `ad` applies them from the last
generator to the first, so undoing a sequence means reversing both:

```julia
gens, angles = PauliList(["Y--", "-YZ", "ZZ-"]), [0.4, -0.9, 1.3]

back = ad(ad(H, gens, angles), reverse(gens), -reverse(angles))
maximum(abs(back[k] - get(H, k, 0.0im)) for k in keys(back))   # 5.6e-17
```

`ad!` does the same in place, for a loop that rotates the same operator over and over.

## What is in the box

| | |
|:--|:--|
| `UPauli` | one string, bits only |
| `Pauli` | one string with a phase `±1, ±i` |
| `PauliList` | an ordered list of strings — a generating set, an algebra |
| `PauliSentence` | strings with coefficients — an operator |
| `com` | commutation test, returning the product string with it |
| `ad`, `ad!` | conjugation by `∏ⱼ exp(iθⱼAⱼ)` |
| `trace` | trace of an operator, without the matrix |
| `countx`, `county`, `countz`, `counti` | bit counts on a string |
| `tostring`, `tomatrix` | text and dense matrix, for reading and checking |

`*`, `+`, `-`, `^` and `cis` come from `Base` and are extended to all of these. `PauliList`
is an `AbstractVector` of the raw strings and `PauliSentence` an `AbstractDict` from string
to coefficient, so the whole `Base` vocabulary — `push!`, `filter!`, `get`, iteration —
applies to both.

The integer type is a parameter: `UInt` covers 32 qubits, `UInt128` covers 64, and
`UPauli{UInt8}` is available for a four-qubit string that has to fit in a byte.

## One convention to know

A `PauliSentence` coefficient multiplies the **bare** embedding of its string — `X` ↦ `σ₁`,
`Z` ↦ `σ₃`, and a `Y` position ↦ the *real* antisymmetric `-iσ₂`. That is what makes a
product of two strings an `xor` and a `±1`, with no `i` anywhere.

Constructing a sentence from anything that names its strings folds in the `i^#Y` phase that
turns each bare string into the Hermitian Pauli operator, so a real coefficient vector gives
a Hermitian operator, and `tostring` divides that phase back out — what you read is what you
entered. Only the raw-bits constructor takes coefficients exactly as given.

## Documentation

Full documentation, a manual and the API reference are built with
[Documenter.jl](https://documenter.juliadocs.org):

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

then open `docs/build/index.html`. Every example in the documentation is executed during that
build, and every docstring example is a doctest verified against the code.

## Tests

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite runs one `@safetestset` per source file, and checks the algebra — products,
commutators, phases, rotations — against the dense matrices `tomatrix` produces, exhaustively
over every string on up to three qubits.

## Related

[RedCarD.jl](https://github.com/ooalshei/RedCarD.jl) builds Cartan decompositions of Pauli
algebras and the fixed-depth quantum circuits they produce, on top of this representation.
