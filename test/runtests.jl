using SafeTestsets

# One safetestset per source file: each runs in a module of its own, so a name defined by
# one group of tests cannot leak into another and change what a later test means.

@safetestset "Pauli strings (pauli.jl)" begin
    include("pauli_tests.jl")
end

@safetestset "PauliList (paulilist.jl)" begin
    include("paulilist_tests.jl")
end

@safetestset "PauliSentence (paulisentence.jl)" begin
    include("paulisentence_tests.jl")
end

@safetestset "Pauli algebra (paulimath.jl)" begin
    include("paulimath_tests.jl")
end

@safetestset "Counts and conversions (utils.jl)" begin
    include("utils_tests.jl")
end
