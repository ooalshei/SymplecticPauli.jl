using SymplecticPauli
using Test

@testset "construction" begin
    v = PauliList(["XX", "ZZ", "-Y"])
    @test v isa PauliList{UInt,2}
    @test length(v) == 3
    @test v.qubits == 2
    @test eltype(v) == UInt
    @test v[1] == UPauli("XX").string
    @test tostring(v) == ["XX", "ZZ", "-Y"]
    # The four ways of naming the same list.
    @test v == PauliList([UPauli("XX"), UPauli("ZZ"), UPauli("-Y")])
    @test v == PauliList([[2, 2], [4, 4], [1, 3]])
    @test v == PauliList([2 4 1; 2 4 3])              # one string per column
    @test v == PauliList(UInt[3, 12, 10], 2)
    @test PauliList{UInt16}(["XX", "ZZ"]) isa PauliList{UInt16,2}
    @test PauliList{UInt16}(["XX", "ZZ"])[1] === UInt16(3)
    @test PauliList{UInt8}([[2, 2], [4, 4]])[2] === UInt8(12)
    @test isempty(PauliList{UInt,2}())
    @test length(PauliList{UInt,2}(undef, 3)) == 3
    @test length(PauliList{UInt}(undef, 3, 2)) == 3
    @test PauliList(undef, 3, 2) isa PauliList{UInt,2}
    @test_throws ArgumentError PauliList(UInt[16], 1)  # element wider than 2Q bits
end

@testset "copying and ownership" begin
    strings = UInt[3, 12]
    v = PauliList(strings, 2)
    strings[1] = 0
    @test v[1] == 3                                    # copied by default
    shared = UInt[3, 12]
    w = PauliList(shared, 2, iscopy=false)
    shared[1] = 0
    @test w[1] == 0                                    # ownership taken
    @test PauliList(UInt[16], 1, check=false)[1] == 16 # the bound check can be skipped
    c = copy(v)
    @test c == v
    @test c.strings !== v.strings
    push!(c, UInt(0))
    @test length(v) == 2 && length(c) == 3
end

@testset "vector interface" begin
    v = PauliList(["XX", "ZZ", "-Y"])
    @test collect(v) == [UPauli("XX").string, UPauli("ZZ").string, UPauli("-Y").string]
    @test size(v) == (3,)
    @test firstindex(v) == 1 && lastindex(v) == 3
    @test v[end] == UPauli("-Y").string
    @test v[2:3] == PauliList(["ZZ", "-Y"])
    @test v[2:3] isa PauliList
    @test v[[1, 3]] == PauliList(["XX", "-Y"])
    @test [tostring(UPauli(s, v.qubits)) for s in v] == ["XX", "ZZ", "-Y"]

    w = copy(v)
    @test push!(w, UPauli("YY").string) === w
    @test w[end] == UPauli("YY").string
    @test append!(w, PauliList(["X-"])) === w
    @test length(w) == 5
    w[1] = UPauli("--").string
    @test w[1] == 0
    @test pop!(w) == UPauli("X-").string
    @test popfirst!(w) == 0
    @test popat!(w, 1) == UPauli("ZZ").string
    @test length(deleteat!(copy(w), 1)) == length(w) - 1
    @test length(resize!(copy(w), 1)) == 1
    @test isempty(empty!(copy(w)))
    @test sizehint!(copy(w), 100) isa PauliList
    @test length(unique(PauliList(["XX", "XX", "ZZ"]))) == 2
    @test unique(PauliList(["XX", "XX"])) isa PauliList
end

@testset "similar" begin
    v = PauliList(["XX", "ZZ", "-Y"])
    @test similar(v) isa PauliList{UInt,2}
    @test length(similar(v)) == 3
    @test length(similar(v, 5)) == 5
    @test length(similar(v, Base.OneTo(4))) == 4
    @test similar(v, UInt16) isa PauliList{UInt16,2}
    @test similar(v, UInt16, 2) isa PauliList{UInt16,2}
    @test length(similar(v, UInt16, 2)) == 2
end

@testset "show and tostring" begin
    v = PauliList(["XX", "-Y"])
    @test tostring(v) == ["XX", "-Y"]
    @test sprint(show, v) == sprint(show, ["XX", "-Y"])
    @test tostring(PauliList{UInt,2}()) == String[]
end
