@doc raw"""
    SymplecticPauli

Pauli operators as symplectic bit strings, and the algebra that comes with them.

A Pauli string on ``Q`` qubits is stored in a single unsigned integer of ``2Q`` bits: the
low ``Q`` bits mark the qubits carrying an ``X`` factor, the high ``Q`` bits the qubits
carrying a ``Z`` factor, and a qubit with both bits set carries a ``Y``. Written for the
string ``P = P_1 \otimes \cdots \otimes P_Q``,

```math
P \;\longleftrightarrow\; (z \,|\, x), \qquad
x = \sum_{j} x_j 2^{\,j-1}, \quad z = \sum_{j} z_j 2^{\,j-1},
```

so that qubit ``j`` holds ``I, X, Z, Y`` according to ``(z_j, x_j) = (0,0), (0,1), (1,0),
(1,1)``. Multiplication is then exclusive or on the bits,

```math
P Q \;\propto\; (z_P \oplus z_Q \,|\, x_P \oplus x_Q),
```

and the phase in front is fixed by counting bits: ``\pm 1`` from ``z_P \cdot x_Q`` and a
power of ``i`` from the number of ``Y`` factors. Two strings commute exactly when their
symplectic form ``z_P \cdot x_Q + x_P \cdot z_Q`` is even. Every operation in this package
is built from those word-sized primitives, so products, commutation tests and rotations cost
a few machine instructions each rather than a matrix multiplication.

# Types

| Type | What it is |
|:--|:--|
| [`UPauli`](@ref) | one string, bits only, no phase |
| [`Pauli`](@ref) | one string with a phase ``\pm 1, \pm i`` |
| [`PauliList`](@ref) | an ordered list of strings — a generating set, an algebra |
| [`PauliSentence`](@ref) | strings with coefficients — an operator |

`PauliList` is an `AbstractVector` of the raw strings and `PauliSentence` an `AbstractDict`
from string to coefficient, so the whole `Base` vocabulary (`push!`, `filter!`, `get`,
iteration, …) applies to both.

# Operations

[`com`](@ref) tests commutation and returns the product string with it; `*` multiplies
strings and sentences, `+`, `-` and scalar `*` combine sentences, `cis` builds
``e^{i\theta P}``, and [`ad`](@ref)/[`ad!`](@ref) conjugate a sentence by such a rotation.
[`tostring`](@ref) and [`tomatrix`](@ref) convert back to something readable — text, or a
dense matrix for checking small cases.

# Examples
```jldoctest
julia> p, q = UPauli("XY-Z"), UPauli("Z--Z");

julia> tostring(p * q)      # a product, with the phase it picks up
"(-i)YY--"

julia> com(p, q).second     # do they commute?
false

julia> H = PauliSentence(PauliList(["ZZ-", "-ZZ", "X--"]), [1.0, 1.0, 0.5]);

julia> rotated = ad(H, UPauli("Y--"), pi / 4, atol = 1e-12);

julia> tostring(rotated)["XZ-"]     # a quarter turn about Y has taken ZZ- to -XZ-
-1.0 + 0.0im
```
"""
module SymplecticPauli

export I,
    σ₁,
    σ₂,
    σ₃,
    σ₂real,
    ⊗,
    AbstractPauli,
    UPauli,
    Pauli,
    PauliList,
    PauliSentence,
    tostring,
    com,
    ad,
    ad!,
    trace,
    countx,
    county,
    countz,
    counti,
    tomatrix

using LinearAlgebra: I, norm, tr

"""
    C8

`Complex{Int8}`, the type the ``\\pm 1, \\pm i`` phase of a [`Pauli`](@ref) is stored in.

Exact and one word wide, so a phase never introduces rounding and promotes to whatever
coefficient type a [`PauliSentence`](@ref) carries.
"""
const C8 = Complex{Int8}

"""
    σ₁

The Pauli matrix ``σ_1 = \\begin{pmatrix} 0 & 1 \\\\ 1 & 0 \\end{pmatrix}``, as `Int8`.

The single-qubit factor [`tomatrix`](@ref) uses for an `X`.
"""
const σ₁ = Int8[0 1; 1 0]

"""
    σ₂

The Pauli matrix ``σ_2 = \\begin{pmatrix} 0 & -i \\\\ i & 0 \\end{pmatrix}``, as
`Complex{Int8}`.

The Hermitian ``Y``. [`tomatrix`](@ref) builds its bare strings from the real
[`σ₂real`](@ref) instead and restores the ``i^{\\#Y}`` phase separately, which keeps the
factors integer-valued.
"""
const σ₂ = C8[0 -im; im 0]

"""
    σ₃

The Pauli matrix ``σ_3 = \\begin{pmatrix} 1 & 0 \\\\ 0 & -1 \\end{pmatrix}``, as `Int8`.

The single-qubit factor [`tomatrix`](@ref) uses for a `Z`.
"""
const σ₃ = Int8[1 0; 0 -1]

"""
    σ₂real

The real antisymmetric ``-i σ_2 = \\begin{pmatrix} 0 & -1 \\\\ 1 & 0 \\end{pmatrix}``, as
`Int8`.

The single-qubit factor of the *bare* embedding of a `Y`: `X` ↦ [`σ₁`](@ref), `Z` ↦
[`σ₃`](@ref), `Y` ↦ `σ₂real`. That embedding is exactly the product ``σ_1 σ_3`` the bit
pattern of a `Y` describes, and it is what a [`PauliSentence`](@ref) coefficient multiplies;
the missing ``i`` per `Y` is the phase [`Pauli`](@ref) carries.
"""
const σ₂real = Int8[0 -1; 1 0]

"""
    ⊗

`kron`, the tensor product, written infix: `A ⊗ B`.

Used to assemble the dense matrices of [`tomatrix`](@ref) one qubit at a time.

# Examples
```jldoctest
julia> σ₁ ⊗ σ₃
4×4 Matrix{Int8}:
 0   0  1   0
 0   0  0  -1
 1   0  0   0
 0  -1  0   0
```
"""
const ⊗ = kron

include("pauli.jl")
include("paulilist.jl")
include("paulisentence.jl")
include("paulimath.jl")
include("utils.jl")

end
