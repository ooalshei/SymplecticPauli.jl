```@meta
CurrentModule = SymplecticPauli
```

# Rotations

The operation this representation exists for: conjugating an operator by

```math
U = e^{i\theta A}, \qquad A \text{ a Hermitian Pauli string},
```

which is [`ad`](@ref).

## Why it is cheap

$A^2 = 1$, so the exponential closes on two terms,

```math
e^{i\theta A} = \cos\theta \, \mathbb{1} + i\sin\theta \, A ,
```

and conjugating a single string $P$ by it splits in two cases. If $[P, A] = 0$ the two
factors cancel and $P$ is untouched. If they anticommute,

```math
U P U^\dagger = \cos(2\theta)\, P + i \sin(2\theta) \, A P ,
```

and $AP$ is one `xor` away. So a rotation is a single pass over the sentence: skip the terms
that commute, and for each of the rest move some amplitude onto exactly one partner string.
No matrix exponential, no truncation.

```@example rotations
using SymplecticPauli

H = PauliSentence(["ZZ-", "-ZZ", "X--", "-X-", "--X"], [1.0, 1.0, 0.5, 0.5, 0.5])
tostring(ad(H, UPauli("Y--"), pi / 4, atol = 1e-12))
```

A quarter turn about $Y$ on the first qubit has taken `X--` to `Z--` and `ZZ-` to `-XZ-`,
and left the rest alone.

## Checking it

[`cis`](@ref) builds the rotation itself — as a two-term sentence, exactly — which is how the
identity above can be checked against dense matrices for small cases:

```@example rotations
using LinearAlgebra

A = UPauli("Y--")
θ = 0.37
U = tomatrix(cis(A, θ))

U * tomatrix(H) * U' ≈ tomatrix(ad(H, A, θ))
```

## Sequences and pruning

Given a [`PauliList`](@ref) of generators and a vector of angles, `ad` applies them from the
**last** generator to the first — the order in which $K = \prod_j e^{i\theta_j A_j}$ acts on
what it conjugates. Undoing a sequence therefore means reversing both:

```@example rotations
gens = PauliList(["Y--", "-YZ", "ZZ-"])
angles = [0.4, -0.9, 1.3]

rotated = ad(H, gens, angles)
back = ad(rotated, reverse(gens), -reverse(angles))

maximum(abs(back[k] - get(H, k, 0.0im)) for k in keys(back))
```

A rotation generally lands amplitude on new strings, so a long sequence fills a sentence up.
The `atol` keyword drops any coefficient whose magnitude falls to it or below, as the
rotation writes them — which is what keeps a sequence from accumulating numerical dust:

```@example rotations
length(ad(H, gens, angles)), length(ad(H, gens, angles, atol = 1e-10))
```

With `atol = 0`, only exact zeros are dropped.

## In place

[`ad!`](@ref) does the same without allocating a new sentence, which is the form to use
inside an optimization loop that rotates the same operator over and over. It needs a sentence
that already holds `ComplexF64` coefficients, since a rotation generally produces complex
ones:

```@example rotations
s = PauliSentence{UInt,ComplexF64}(["Z--"], [1.0])
ad!(s, gens, angles, atol = 1e-12)
length(s)
```

Rotating by a list of generators in place applies them last-to-first as well, and does not
copy in between.

## Precomputed angles

Every `ad`/`ad!` method has a four-argument form taking $\cos 2\theta$ and $\sin 2\theta$
directly, for a sweep that reuses one pair of angles across many sentences:

```@example rotations
tostring(ad(H, A, cos(2θ), sin(2θ))) == tostring(ad(H, A, θ))
```
