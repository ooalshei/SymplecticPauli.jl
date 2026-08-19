```@meta
CurrentModule = SymplecticPauli
```

# SymplecticPauli.jl

Pauli operators as **symplectic bit strings** — products, commutators and rotations at the
cost of a few machine instructions.

A Pauli string on $Q$ qubits is one unsigned integer of $2Q$ bits. The low $Q$ bits say which
qubits carry an $X$ factor, the high $Q$ bits which carry a $Z$, and a qubit with both bits
set carries a $Y$:

```math
P = P_1 \otimes \cdots \otimes P_Q \;\longleftrightarrow\; (z \,|\, x), \qquad
(z_j, x_j) = (0,0),\,(0,1),\,(1,0),\,(1,1) \;\text{ for } P_j = I, X, Z, Y .
```

Everything else follows from that. The product of two strings is an exclusive or,

```math
PQ \;\propto\; (z_P \oplus z_Q \,|\, x_P \oplus x_Q),
```

its phase is $(-1)^{z_P \cdot x_Q}$ times a power of $i$ counting the $Y$s, and the two
strings commute exactly when the symplectic form

```math
\langle P, Q \rangle = z_P \cdot x_Q + x_P \cdot z_Q
```

is even. A product is one `xor`; a commutation test is two `and`s, an `xor` and a popcount.
No matrix is ever built — [`tomatrix`](@ref) exists for checking small cases, and that is all
it is for.

## Installation

The package is not registered:

```julia
using Pkg
Pkg.add(url="https://github.com/ooalshei/SymplecticPauli.jl")
```

## Quick start

```jldoctest quickstart
julia> using SymplecticPauli

julia> p, q = UPauli("XY-Z"), UPauli("Z--Z")
(Pauli{UInt64}(XY-Z), Pauli{UInt64}(Z--Z))

julia> p * q                          # a product, and the phase it picks up
Pauli{UInt64}((-i)YY--)

julia> com(p, q).second               # the same string, plus: do they commute?
false
```

An operator is a [`PauliSentence`](@ref) — strings with coefficients:

```jldoctest quickstart
julia> H = PauliSentence(["ZZ-", "-ZZ", "X--", "-X-", "--X"], [1.0, 1.0, 0.5, 0.5, 0.5]);

julia> length(H), H.qubits
(5, 3)

julia> trace(H)
0.0
```

and rotating it by $e^{i\theta A}$ is [`ad`](@ref), which touches only the terms that
anticommute with the generator:

```jldoctest quickstart
julia> rotated = ad(H, UPauli("Y--"), pi / 4, atol = 1e-12);

julia> tostring(rotated)["Z--"]       # X-- has been turned into Z--
0.5 + 0.0im
```

## What is in the box

| | |
|:--|:--|
| [`UPauli`](@ref) | one string, bits only |
| [`Pauli`](@ref) | one string with a phase $\pm 1, \pm i$ |
| [`PauliList`](@ref) | an ordered list of strings — a generating set, an algebra |
| [`PauliSentence`](@ref) | strings with coefficients — an operator |
| [`com`](@ref) | commutation test, returning the product string with it |
| [`ad`](@ref), [`ad!`](@ref) | conjugation by $\prod_j e^{i\theta_j A_j}$ |
| [`trace`](@ref) | trace of an operator, without the matrix |
| [`countx`](@ref), [`county`](@ref), [`countz`](@ref), [`counti`](@ref) | bit counts on a string |
| [`tostring`](@ref), [`tomatrix`](@ref) | text and dense matrix, for reading and checking |

`*`, `+`, `-`, `^` and `cis` come from `Base` and are extended to all of these.

## Where to go next

The manual works through the three layers in order — [Pauli strings](@ref), the
[Operators](@ref) built from them, and the [Rotations](@ref) that act on those — and the
[API reference](@ref) lists everything with its signatures.

This package is the representation layer of
[RedCarD.jl](https://github.com/ooalshei/RedCarD.jl), which builds Cartan decompositions of
Pauli algebras and the fixed-depth circuits they produce; the operations here are the ones
that inner loop needs.
