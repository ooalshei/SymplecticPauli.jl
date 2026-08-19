using SymplecticPauli
using SymplecticPauli: _check_type, _check_string_length
using Test

# `tomatrix` is the definition of the representation, so the identities below are checked
# against it: a `Pauli`'s `sign` is what turns the bare symplectic string into the Hermitian
# Pauli operator.

"All `2Q`-bit strings on `Q` qubits."
allstrings(Q) = UInt.(0:(4^Q-1))

@testset "construction from text and indices" begin
    @test UPauli("XZ").string == UPauli([2, 4]).string
    @test UPauli("XYZ-").string == UPauli([2, 3, 4, 1]).string
    @test UPauli("--").string == 0
    @test UPauli("II").string == UPauli("--").string     # 'I' and '-' both mean identity
    @test UPauli("XY-Z").string == 0xa3
    @test tostring(UPauli([2, 3, 1, 4])) == "XY-Z"
    # The leftmost character is qubit 1: the low bit of the X half.
    @test UPauli("X-").string == 1
    @test UPauli("-X").string == 2
    @test UPauli("Z-").string == 4
    @test UPauli("Y-").string == 5
end

@testset "storage type" begin
    @test UPauli("XY-Z") isa UPauli{UInt,4}
    @test UPauli("XY-Z") isa AbstractPauli{UInt,4}
    @test Pauli("XY-Z") isa Pauli{UInt,4}
    @test isbitstype(UPauli{UInt,3})
    @test isbitstype(Pauli{UInt,3})
    # A narrower type has to give the same bits, not silently promote to `Int` on the way.
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
        @test UPauli{T}("XY-Z").string === T(0xa3)
        @test UPauli{T}([2, 3, 1, 4]).string === T(0xa3)
        @test UPauli{T}("XY-Z") isa UPauli{T,4}
        @test Pauli{T}("XY-Z").string === T(0xa3)
    end
    # A string needs 2 bits per qubit, so `UInt8` tops out at four qubits.
    @test _check_type(UInt8, 4) === nothing
    @test_throws ArgumentError _check_type(UInt8, 5)
    @test_throws ArgumentError UPauli{UInt8}("XXXXX")
    @test_throws ArgumentError UPauli{UInt8}([2, 2, 2, 2, 2])
    # Wide types have to work all the way out to their last bit.
    @test UPauli(typemax(UInt32), 16).string == typemax(UInt32)
    @test UPauli(typemax(UInt128), 64).qubits == 64
end

@testset "conversion between the two types" begin
    p = UPauli("XY-Z")
    @test UPauli(Pauli(p)) === p                  # the phase is dropped, the bits are kept
    @test Pauli(p).string == p.string
    @test Pauli(Pauli(p)) === Pauli(p)
    @test UPauli{UInt16}(p).string === UInt16(p.string)
    @test Pauli{UInt16}(p) isa Pauli{UInt16,4}
    @test Pauli{UInt16}(Pauli(p)).sign == Pauli(p).sign
    @test copy(p) === p                           # immutable: a copy is the thing itself
end

@testset "the phase" begin
    @test Pauli(UPauli("Y")).sign == im
    @test Pauli(UPauli("YY")).sign == -1
    @test Pauli(UPauli("YYY")).sign == -im
    @test Pauli(UPauli("XZ")).sign == 1
    @test Pauli("Y").sign == im
    @test Pauli(UPauli("Y"), -1).sign == -1      # an explicit sign beats the convention
    @test Pauli("Y", -1).sign == -1
    @test Pauli([3], -1).sign == -1
    @test Pauli(UInt(1), 2) isa Pauli{UInt,2}
    @test Pauli(UInt(1), 2).sign == 1
    for Q in 1:3, s in allstrings(Q)
        # The convention phase is exactly the one that makes the operator Hermitian.
        m = tomatrix(UPauli(s, Q))
        @test m ≈ adjoint(m)
        @test Pauli(UPauli(s, Q)).sign == im^county(UPauli(s, Q))
    end
end

@testset "validation" begin
    @test_throws ArgumentError UPauli("XW")
    @test_throws ArgumentError UPauli("xy")
    @test_throws ArgumentError UPauli([1, 5])
    @test_throws ArgumentError UPauli([0, 1])
    @test_throws ArgumentError UPauli(UInt(4^2), 1)   # string wider than 2Q bits
    @test_throws ArgumentError Pauli(UInt(4^2), 1, 1)
    @test_throws ArgumentError Pauli(UInt(1), 2, 1)   # sign not a fourth root of 1
    @test_throws ArgumentError Pauli(UInt(1), 0, 1)
    @test _check_string_length(UInt(3), 1) === nothing
    # Full-width strings must not be rejected by an overflowing bound check.
    @test _check_string_length(typemax(UInt), 32) === nothing
    @test _check_string_length(typemax(UInt128), 64) === nothing
end

# The default storage type is `UInt`, which is 32 bits wide on a 32-bit platform, so the
# expected text interpolates it rather than naming `UInt64`.
@testset "show" begin
    @test sprint(show, UPauli("XY-Z")) == "Pauli{$UInt}(XY-Z)"
    @test sprint(show, UPauli{UInt16}("XY-Z")) == "Pauli{UInt16}(XY-Z)"
    @test sprint(show, Pauli("Z")) == "Pauli{$UInt}((+)Z)"
    @test sprint(show, Pauli("Z", -1)) == "Pauli{$UInt}((-)Z)"
    @test sprint(show, Pauli("Z", im)) == "Pauli{$UInt}((i)Z)"
    @test sprint(show, Pauli("Z", -im)) == "Pauli{$UInt}((-i)Z)"
end

@testset "equality and hashing" begin
    @test UPauli("XY") == UPauli("XY")
    @test UPauli("XY") != UPauli("YX")
    @test Pauli("XY", 1) != Pauli("XY", -1)
    @test length(Set([UPauli("XY"), UPauli("XY"), UPauli("YX")])) == 2
end
