module SymplecticPauli

const C8 = Complex{Int8}

export AbstractPauli, Pauli, SignedPauli, PauliSentence, tostring, com, ad

include("zerorange.jl")
include("pauli.jl")
include("paulisentence.jl")
include("broadcast.jl")
include("paulimath.jl")
include("utils.jl")

end


