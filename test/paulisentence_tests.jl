using SymplecticPauli
using LinearAlgebra: I, tr
using Test

# A `PauliSentence` coefficient multiplies the *bare* embedding of its string (X ↦ σ₁,
# Z ↦ σ₃, Y position ↦ σ₂real), so `tomatrix` is what every claim below is checked against.

dense(s::PauliSentence) = tomatrix(s)

@testset "construction from raw strings" begin
    s = PauliSentence(UInt[1, 6], [1.0, 2.0], 2)
    @test s isa PauliSentence{UInt,Float64,2}
    @test s.qubits == 2
    @test length(s) == 2
    @test s[UInt(1)] == 1.0
    @test PauliSentence{UInt,ComplexF64}(UInt[1], [1.0], 2) isa
          PauliSentence{UInt,ComplexF64,2}
    @test PauliSentence{UInt,ComplexF64,2}(UInt[1], [1.0]) isa
          PauliSentence{UInt,ComplexF64,2}
    # Raw strings are taken as given: no convention phase is folded in.
    @test PauliSentence(UInt[3], [1.0], 1)[UInt(3)] == 1.0
    @test_throws DimensionMismatch PauliSentence(UInt[1, 2], [1.0], 1)
    @test_throws ArgumentError PauliSentence{UInt,ComplexF64,1}(UInt[4], [1.0])
end

@testset "construction from named strings" begin
    # Anything that names its strings folds in the i^#Y that makes each term Hermitian, so a
    # real coefficient vector gives a Hermitian operator.
    for ctor in (
        c -> PauliSentence(PauliList(["Y-", "ZZ"]), c),
        c -> PauliSentence(["Y-", "ZZ"], c),
        c -> PauliSentence([[3, 1], [4, 4]], c),
        c -> PauliSentence([3 4; 1 4], c),
        c -> PauliSentence([UPauli("Y-"), UPauli("ZZ")], c),
    )
        s = ctor([1.0, -0.5])
        @test dense(s) ≈ adjoint(dense(s))
        @test dense(s) ≈ tomatrix(Pauli("Y-")) - 0.5 * tomatrix(Pauli("ZZ"))
        @test tostring(s)["Y-"] ≈ 1.0
    end
    # An odd number of Ys makes the stored coefficient complex, whatever came in.
    @test PauliSentence(["Y-"], [1.0])[UPauli("Y-").string] == im
    @test valtype(PauliSentence(["Y-"], [1.0])) <: Complex
    # An even number leaves a real vector real.
    @test valtype(PauliSentence(["YY"], [2.0])) == Float64
    @test PauliSentence(["YY"], [2.0])[UPauli("YY").string] == -2.0
    @test PauliSentence{UInt16,ComplexF64}(["XY"], [1.0]) isa
          PauliSentence{UInt16,ComplexF64,2}
    @test PauliSentence(PauliList(["Y-"]), [1.0])[UPauli("Y-").string] == im
end

@testset "construction from dictionaries and sentences" begin
    @test PauliSentence(Dict(UInt(1) => 2.0), 1)[UInt(1)] == 2.0
    @test PauliSentence(Dict("Y-" => 1.0))[UPauli("Y-").string] == im
    @test PauliSentence(Dict(UPauli("Y-") => 1.0))[UPauli("Y-").string] == im
    @test PauliSentence(Dict([3, 1] => 1.0))[UPauli("Y-").string] == im
    d = Dict(UInt(1) => 2.0)
    s = PauliSentence(d, 1)
    d[UInt(1)] = 0.0
    @test s[UInt(1)] == 2.0                        # copied by default
    shared = Dict(UInt(1) => 2.0)
    t = PauliSentence(shared, 1, iscopy=false)
    shared[UInt(1)] = 0.0
    @test t[UInt(1)] == 0.0                        # ownership taken
    u = PauliSentence(["ZZ"], [1.0])
    @test PauliSentence{UInt16,ComplexF64}(u) isa PauliSentence{UInt16,ComplexF64,2}
    @test PauliSentence(u) == u
    @test PauliSentence(u).sentence !== u.sentence
    @test_throws ArgumentError PauliSentence(Dict(UInt(4) => 1.0), 1)
end

@testset "construction from a matrix" begin
    for m in (
        [1.0 0.0; 0.0 -1.0],                       # Z
        [0.0 -im; im 0.0],                         # Y
        ComplexF64[1 2; 3 4],                      # not Hermitian: still decomposable
    )
        s = PauliSentence(m)
        @test s isa PauliSentence{UInt,ComplexF64,1}
        @test dense(s) ≈ m
    end
    s = PauliSentence(PauliList(["ZZ", "X-", "-Y"]), [1.0, 2.0, 3.0])
    @test dense(PauliSentence(dense(s))) ≈ dense(s)
    @test PauliSentence{UInt16,ComplexF64}(dense(s)) isa PauliSentence{UInt16,ComplexF64,2}
    @test length(PauliSentence(zeros(ComplexF64, 4, 4))) == 0
    @test_throws ArgumentError PauliSentence(ComplexF64[1 2 3; 4 5 6])
    @test_throws ArgumentError PauliSentence(ComplexF64[1 2 3; 4 5 6; 7 8 9])
end

@testset "dictionary interface" begin
    Q = 2
    s = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 6], ComplexF64[1.0, 2.0])
    @test length(s) == 2
    @test !isempty(s)
    @test haskey(s, UInt(1)) && !haskey(s, UInt(3))
    @test s[UInt(6)] == 2.0
    @test get(s, UInt(3), 0.0im) == 0.0im
    @test sort(collect(keys(s))) == UInt[1, 6]
    @test sum(values(s)) == 3.0
    @test keytype(s) == UInt && valtype(s) == ComplexF64
    @test sort([k for (k, _) in s]) == UInt[1, 6]
    @test sizehint!(s, 10) === s

    c = copy(s)
    c[UInt(3)] = 5.0
    @test length(s) == 2 && length(c) == 3          # the copy is independent
    @test isempty(empty(s)) && !isempty(s)
    @test empty(s) isa typeof(s)
    @test length(filter(p -> real(p.second) > 1.5, s)) == 1
    @test typeof(filter(p -> true, s)) === typeof(s)
    @test length(filter!(p -> real(p.second) > 1.5, c)) == 2
    @test delete!(c, UInt(3)) === c
    @test !haskey(c, UInt(3))
    @test pop!(c, UInt(6)) == 2.0
    @test pop!(c, UInt(6), -1.0) == -1.0
    @test_throws KeyError pop!(c, UInt(6))
end

@testset "tostring and show" begin
    s = PauliSentence(PauliList(["ZZ", "-Y"]), [0.5, 1.0])
    @test tostring(s) == Dict("ZZ" => 0.5 + 0.0im, "-Y" => 1.0 + 0.0im)
    @test sprint(show, s) == sprint(show, tostring(s))
    # The convention phase goes in on construction and comes back out here.
    @test tostring(PauliSentence(["YYY"], [1.0]))["YYY"] ≈ 1.0
end

@testset "trace" begin
    Q = 2
    s = PauliSentence(PauliList(["--", "ZZ"]), [0.5, 1.0])
    @test trace(s) ≈ 2^Q * s[UInt(0)]
    @test trace(s) ≈ tr(dense(s))
    @test trace(PauliSentence(PauliList(["ZZ"]), [1.0])) == 0
    @test trace(PauliSentence{UInt,ComplexF64,2}(UInt[], ComplexF64[])) === zero(ComplexF64)
    @test trace(PauliSentence(PauliList(["ZZ"]), [1.0])) isa Number
end
