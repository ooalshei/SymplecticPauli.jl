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

const C8 = Complex{Int8}
const σ₁ = Int8[0 1; 1 0]
const σ₂ = C8[0 -im; im 0]
const σ₃ = Int8[1 0; 0 -1]
const σ₂real = Int8[0 -1; 1 0]
const ⊗ = kron

include("pauli.jl")
include("paulilist.jl")
include("paulisentence.jl")
include("paulimath.jl")
include("utils.jl")

end
