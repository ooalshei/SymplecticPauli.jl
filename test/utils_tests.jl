using SymplecticPauli
using SymplecticPauli: _check_string_length, _check_type
using LinearAlgebra: I
using Test

"All `2Q`-bit strings on `Q` qubits."
allstrings(Q) = UInt.(0:(4^Q-1))

# A `Y` sets both the X and the Z bit of its qubit, so it is counted by `countx` and by
# `countz` as well as by `county`. `counti` is the odd one out: it folds the two halves
# together, so a `Y` is one non-identity qubit like any other factor and the count is the
# number of identity factors the printed string shows.
@testset "counting" begin
    p = UPauli("XYZI-")
    @test countx(p) == 2
    @test county(p) == 1
    @test countz(p) == 2
    @test counti(p) == 2                              # the I and the - of "XYZI-"
    @test countx(p.string, p.qubits) == countx(p)
    @test county(p.string, p.qubits) == county(p)
    @test countz(p.string, p.qubits) == countz(p)
    @test counti(p.string, p.qubits) == counti(p)
    @test countx(Pauli("XYZI-")) == 2                 # the phase makes no difference
    for Q in 1:3, s in allstrings(Q)
        chars = tostring(UPauli(s, Q))
        @test countx(s, Q) == count(∈("XY"), chars)
        @test county(s, Q) == count(==('Y'), chars)
        @test countz(s, Q) == count(∈("ZY"), chars)
        @test counti(s, Q) == count(==('-'), chars)
        # …which is `Q` minus the weight, and the three halves-counts overcount a `Y` once.
        @test counti(s, Q) ==
              Q - countx(s, Q) - countz(s, Q) + county(s, Q)
    end
    @test_throws ArgumentError countx(UInt(4), 1)
    @test_throws ArgumentError county(UInt(4), 1)
    @test_throws ArgumentError countz(UInt(4), 1)
    @test_throws ArgumentError counti(UInt(4), 1)
    # The masks have to be built in the string's own type: an `Int` one overflows at Q = 64.
    wide = UPauli(typemax(UInt128), 64)
    @test countx(wide) == 64
    @test county(wide) == 64
    @test countz(wide) == 64
    @test countx(UPauli(UInt128(1) << 64, 64)) == 0   # a lone Z on the last qubit
    @test countz(UPauli(UInt128(1) << 64, 64)) == 1
    @test counti(wide) == 0                           # every qubit carries a Y
    @test counti(UPauli(UInt128(0), 64)) == 64
end

@testset "tostring" begin
    @test tostring(UPauli("XY-Z")) == "XY-Z"
    @test tostring(UPauli("XYIZ")) == "XY-Z"          # 'I' comes back as '-'
    @test tostring(UPauli(UInt(0), 3)) == "---"
    @test tostring(Pauli("Z-")) == "(+)Z-"
    @test tostring(Pauli("Z-", -1)) == "(-)Z-"
    @test tostring(Pauli("Z-", im)) == "(i)Z-"
    @test tostring(Pauli("Z-", -im)) == "(-i)Z-"
    @test tostring(Pauli("Y-")) == "(+)Y-"            # the convention phase reads as +
    @test tostring(PauliList(["ZZ", "-X"])) == ["ZZ", "-X"]
    @test tostring(PauliSentence(PauliList(["ZZ"]), [0.5])) == Dict("ZZ" => 0.5 + 0.0im)
    for Q in 1:3, s in allstrings(Q)
        @test UPauli(tostring(UPauli(s, Q))).string == s
    end
end

@testset "tomatrix" begin
    @test tomatrix(UInt(0), 2) == I(4)
    @test tomatrix(Pauli("X")) == [0 1; 1 0]
    @test tomatrix(Pauli("Y")) == [0 -im; im 0]
    @test tomatrix(Pauli("Z")) == [1 0; 0 -1]
    @test tomatrix(UPauli("Y")) == [0 -im; im 0]
    @test tomatrix(UInt(3), 1) == [0 -1; 1 0]         # the bare embedding of a Y is real
    @test tomatrix(Pauli("X", -1)) == -[0 1; 1 0]
    @test tomatrix(UPauli("XZ")) == kron([0 1; 1 0], [1 0; 0 -1])
    @test size(tomatrix(UPauli(UInt(0), 4))) == (16, 16)
    @test_throws ArgumentError tomatrix(UInt(4), 1)
    for Q in 1:3, s in allstrings(Q)
        m = tomatrix(UPauli(s, Q))
        @test m ≈ adjoint(m)                          # Hermitian
        @test m * m ≈ I(2^Q)                          # squares to the identity
        # …and the bare form is the phaseless one
        @test tomatrix(s, Q) ≈ (-im)^county(s, Q) * m
    end
    s = PauliSentence(PauliList(["ZZ", "X-"]), [1.0, -0.5])
    @test tomatrix(s) ≈ tomatrix(Pauli("ZZ")) - 0.5 * tomatrix(Pauli("X-"))
    @test tomatrix(PauliSentence{UInt,ComplexF64,2}(UInt[], ComplexF64[])) ==
          zeros(ComplexF64, 4, 4)
end
