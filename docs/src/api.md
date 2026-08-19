```@meta
CurrentModule = SymplecticPauli
```

# API reference

```@docs
SymplecticPauli
```

```@contents
Pages = ["api.md"]
Depth = 3
```

## Types

```@docs
AbstractPauli
UPauli
Pauli
PauliList
PauliSentence
```

## Products and commutation

```@docs
com
unsigned_prod
Base.:*(::UPauli{<:Unsigned,Q}, ::UPauli{<:Unsigned,Q}) where {Q}
Base.:*(::Number, ::UPauli)
```

## Operator arithmetic

```@docs
Base.:+(::PauliSentence{Ts,Ns,Q}, ::PauliSentence{Tr,Nr,Q}) where {Ts,Ns,Tr,Nr,Q}
Base.:*(::Number, ::PauliSentence)
Base.:*(::PauliSentence{Ts,<:Number,Q}, ::PauliSentence{Tr,<:Number,Q}) where {Ts,Tr,Q}
Base.:^(::PauliSentence, ::Integer)
trace
```

## Rotations and conjugation

```@docs
Base.cis(::UPauli, ::Real)
ad
ad!
```

## Counting

```@docs
countx
county
countz
counti
```

## Conversions

```@docs
tostring
tomatrix
```

## Single-qubit matrices

```@docs
σ₁
σ₂
σ₃
σ₂real
⊗
```

`using SymplecticPauli` also brings in `LinearAlgebra`'s `I`, which
[`tomatrix`](@ref) uses for the identity factors.

## Internals

Not exported, and not part of the public interface — documented because the packages built on
this one reach for them.

```@docs
SymplecticPauli.C8
SymplecticPauli._check_string_length
SymplecticPauli._check_type
SymplecticPauli._check_keys
SymplecticPauli._ipow
SymplecticPauli._mpow
SymplecticPauli._symplectic_prod
SymplecticPauli._prodsign
SymplecticPauli._prune!
SymplecticPauli._prunable
```

## Index

```@index
```
