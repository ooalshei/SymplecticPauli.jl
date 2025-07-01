module SymplecticPauli

export I, σx, σy, σz, σy_real, ⊗, AbstractPauli, UPauli, Pauli, PauliList, PauliSentence, toint, tostring, com, ad, ad!, countx, county, countz, counti, tomatrix

using LinearAlgebra: I, norm, tr

const C8 = Complex{Int8}
const σx = Int8[0 1; 1 0]
const σy = C8[0 -im; im 0]
const σz = Int8[1 0; 0 -1]
const σy_real = Int8[0 -1; 1 0]
const ⊗ = kron

include("pauli.jl")
include("paulilist.jl")
include("paulisentence.jl")
include("paulimath.jl")
include("utils.jl")

end
