using SymplecticPauli
using SymplecticPauli: _ipow, _mpow, _prodsign, _symplectic_prod, _check_type
using LinearAlgebra: I, tr
using Test

# Every algebraic identity below is checked against the dense matrices produced by
# `tomatrix`, which is the definition of the representation: a `PauliSentence` coefficient
# multiplies the *bare* embedding of its string (X ↦ σ₁, Z ↦ σ₃, Y position ↦ σ₂real), and a
# `Pauli`'s `sign` turns that bare string into the Hermitian Pauli operator.

"All `2Q`-bit strings on `Q` qubits."
allstrings(Q) = UInt.(0:(4^Q-1))

dense(s::PauliSentence) = tomatrix(s)

@testset "SymplecticPauli.jl" begin

    @testset "phase helpers" begin
        for k in 0:16
            @test _ipow(k) == im^k
            @test _mpow(k) == (-1)^k
        end
        for Q in 1:3, p in allstrings(Q), q in allstrings(Q)
            @test _prodsign(p, q, Q) == _symplectic_prod(p, q, Q).second
            @test _symplectic_prod(p, q, Q).first == p ⊻ q
        end
    end

    # A `Y` sets both the X and the Z bit of its qubit, so it is counted by `countx` and by
    # `countz` as well as by `county`; `counti` is `Q` minus the three, which weighs a `Y`
    # three times over. `subalgfind` relies on that ordering, so the convention is pinned
    # here rather than reinterpreted.
    @testset "counting" begin
        p = UPauli("XYZI-")
        @test countx(p) == 2
        @test county(p) == 1
        @test countz(p) == 2
        @test counti(p) == 5 - 2 - 1 - 2
        for Q in 1:3, s in allstrings(Q)
            chars = tostring(UPauli(s, Q))
            @test countx(s, Q) == count(∈("XY"), chars)
            @test county(s, Q) == count(==('Y'), chars)
            @test countz(s, Q) == count(∈("ZY"), chars)
            @test counti(s, Q) == Q - countx(s, Q) - county(s, Q) - countz(s, Q)
        end
    end

    @testset "construction and validation" begin
        @test UPauli("XZ").string == UPauli([2, 4]).string
        @test UPauli("--").string == 0
        @test Pauli(UPauli("Y")).sign == im
        @test_throws ArgumentError UPauli("XW")
        @test_throws ArgumentError UPauli([1, 5])
        @test_throws ArgumentError UPauli(UInt(4^2), 1)  # string wider than 2Q bits
        @test_throws ArgumentError Pauli(UInt(1), 2, 1)  # sign not a fourth root of 1
        # A string needs 2 bits per qubit, so `UInt8` tops out at four qubits.
        @test _check_type(UInt8, 4) === nothing
        @test_throws ArgumentError _check_type(UInt8, 5)
        # Full-width strings must not be rejected by an overflowing bound check.
        @test UPauli(typemax(UInt32), 16).string == typemax(UInt32)
        @test isbitstype(UPauli{UInt,3})
        @test isbitstype(Pauli{UInt,3})
    end

    @testset "tostring / tomatrix" begin
        @test tostring(UPauli("XY-Z")) == "XY-Z"
        @test tomatrix(UInt(0), 2) == I(4)
        @test tomatrix(Pauli("X")) == [0 1; 1 0]
        @test tomatrix(Pauli("Y")) == [0 -im; im 0]
        @test tomatrix(Pauli("Z")) == [1 0; 0 -1]
        for Q in 1:3, s in allstrings(Q)
            m = tomatrix(UPauli(s, Q))
            @test m ≈ adjoint(m)                     # Hermitian
            @test m * m ≈ I(2^Q)                     # squares to the identity
        end
    end

    @testset "Pauli products" begin
        for Q in 1:2, p in allstrings(Q), q in allstrings(Q)
            up, uq = UPauli(p, Q), UPauli(q, Q)
            expected = tomatrix(up) * tomatrix(uq)
            @test tomatrix(up * uq) ≈ expected
            @test tomatrix(Pauli(up) * uq) ≈ expected
            @test tomatrix(up * Pauli(uq)) ≈ expected
            @test tomatrix(Pauli(up) * Pauli(uq)) ≈ expected
            # `com` reports whether the two commute, alongside their product string.
            @test com(up, uq).second ==
                  isapprox(tomatrix(up) * tomatrix(uq), tomatrix(uq) * tomatrix(up))
            @test com(p, q, Q).second == com(up, uq).second
            @test com(p, q, Q).first == (up ⊻ uq).string
        end
    end

    @testset "sentence arithmetic" begin
        Q = 3
        a = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 9, 26], ComplexF64[0.3, -0.7im, 0.2])
        b = PauliSentence{UInt,ComplexF64,Q}(UInt[9, 5], ComplexF64[1.5, 0.25])
        @test dense(a + b) ≈ dense(a) + dense(b)
        @test dense(a - b) ≈ dense(a) - dense(b)
        @test dense(a + b + a) ≈ 2 * dense(a) + dense(b)
        @test dense(2.5 * a) ≈ 2.5 * dense(a)
        @test dense(a * 2.5) ≈ 2.5 * dense(a)
        @test dense((1 + 2im) * a) ≈ (1 + 2im) * dense(a)
        @test dense(a * b) ≈ dense(a) * dense(b)
        @test dense(a^3) ≈ dense(a)^3
        @test dense(a^0) ≈ I(2^Q)
        negative = -1
        @test_throws ArgumentError a^negative
        # Coefficients live in the bare basis, where the identity carries the whole trace.
        @test tr(dense(a * a)) ≈ 2^Q * (a*a)[UInt(0)]
    end

    @testset "sentence container" begin
        Q = 2
        s = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 6], ComplexF64[1.0, 2.0])
        @test length(s) == 2
        @test haskey(s, UInt(1)) && !haskey(s, UInt(3))
        @test s[UInt(6)] == 2.0
        @test get(s, UInt(3), 0.0im) == 0.0im
        @test sort(collect(keys(s))) == UInt[1, 6]
        @test sum(values(s)) == 3.0
        c = copy(s)
        c[UInt(3)] = 5.0
        @test length(s) == 2 && length(c) == 3     # the copy is independent
        @test isempty(empty(s)) && !isempty(s)
        @test length(filter(p -> real(p.second) > 1.5, s)) == 1
        @test length(filter!(p -> real(p.second) > 1.5, c)) == 2
        @test typeof(filter(p -> true, s)) === typeof(s)
        delete!(c, UInt(3))
        @test !haskey(c, UInt(3))
        @test pop!(c, UInt(6)) == 2.0
        @test pop!(c, UInt(6), -1.0) == -1.0
        @test_throws ArgumentError PauliSentence{UInt,ComplexF64,1}(UInt[4], [1.0])
        # Round trip through a dense matrix.
        @test dense(PauliSentence{UInt,ComplexF64}(dense(s))) ≈ dense(s)
    end

    @testset "PauliList" begin
        v = PauliList(["XX", "ZZ", "-Y"])
        @test length(v) == 3
        @test v[1] == UPauli("XX").string
        @test tostring(v) == ["XX", "ZZ", "-Y"]
        @test v[2:3] == PauliList(["ZZ", "-Y"])
        @test copy(v) == v && copy(v).strings !== v.strings
        push!(v, UPauli("YY").string)
        @test length(v) == 4 && v[end] == UPauli("YY").string
        append!(v, PauliList(["X-"]))
        @test length(v) == 5
        @test length(deleteat!(copy(v), 1)) == 4
        @test length(unique(PauliList(["XX", "XX"]))) == 1
        @test v.qubits == 2
        @test eltype(v) == UInt
        @test_throws ArgumentError PauliList(UInt[16], 1)
        @test isempty(PauliList{UInt,2}())
        @test length(PauliList{UInt,2}(undef, 3)) == 3
    end

    @testset "rotations" begin
        Q = 3
        s = PauliSentence{UInt,ComplexF64,Q}(
            UInt[1, 9, 26, 63],
            ComplexF64[0.3, -0.7im, 0.2, 1.1],
        )
        for gs in UInt[5, 9, 63], θ in (0.0, 0.37, pi / 4, -1.2)
            g = UPauli(gs, Q)
            # exp(iθA)·s·exp(iθA)† with A the Hermitian Pauli of the generator.
            u = dense(cis(g, θ))
            @test u * adjoint(u) ≈ I(2^Q)
            expected = u * dense(s) * adjoint(u)
            @test dense(ad(s, g, θ)) ≈ expected
            @test dense(ad!(copy(s), g, θ)) ≈ expected
            # …and the two-argument form agrees with the (cosine, sine) form.
            @test dense(ad(s, g, cos(2θ), sin(2θ))) ≈ expected
        end
        # The identity generates a global phase, which conjugation cannot see.
        @test dense(ad(s, UPauli(UInt(0), Q), 0.8)) ≈ dense(s)

        gens = PauliList(["X--", "-YZ", "ZZ-"])
        angles = [0.4, -0.9, 1.3]
        expected = dense(s)
        # `ad` applies the generators from the last to the first.
        for i in reverse(eachindex(gens))
            u = dense(cis(UPauli(gens[i], Q), angles[i]))
            expected = u * expected * adjoint(u)
        end
        @test dense(ad(s, gens, angles)) ≈ expected
        @test dense(ad!(copy(s), gens, angles)) ≈ expected
        @test dense(ad(s, gens, cos.(2angles), sin.(2angles))) ≈ expected
        @test_throws DimensionMismatch ad(s, gens, [0.1, 0.2])
        @test dense(ad(s, PauliList{UInt,Q}(), Float64[])) ≈ dense(s)

        # Undoing a sequence of rotations means walking it backwards: the generators do not
        # commute, so the order has to be reversed along with the angles.
        @test dense(ad(ad(s, gens, angles), reverse(gens), -reverse(angles))) ≈ dense(s)
        @test all(v -> abs(v) > 0.5, values(ad(s, gens, angles, atol=0.5)))
    end
end
