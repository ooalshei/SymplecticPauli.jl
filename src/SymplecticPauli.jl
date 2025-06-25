module SymplecticPauli

export I, σx, σy, σz, ⊗, AbstractPauli, UPauli, Pauli, PauliSentence, toint, tostring, com, ad, ad!, countx, county, countz, tomatrix

using LinearAlgebra: I

const C8 = Complex{Int8}
const σx = [0 1; 1 0]
const σy = [0 -im; im 0]
const σz = [1 0; 0 -1]
const σy_real = -im * σy
const ⊗ = kron

include("pauli.jl")
include("paulisentence.jl")
include("paulimath.jl")
include("utils.jl")

end
