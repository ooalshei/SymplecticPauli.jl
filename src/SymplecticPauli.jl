module SymplecticPauli

const C8 = Complex{Int8}

export AbstractPauli, UPauli, Pauli, PauliSentence, toint, tostring, com, ad, countx, county, countz

include("pauli.jl")
include("paulisentence.jl")
include("paulimath.jl")
include("utils.jl")

end
