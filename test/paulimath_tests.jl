using SymplecticPauli
using SymplecticPauli:
    _ipow, _mpow, _prodsign, _symplectic_prod, _prune!, _prunable, unsigned_prod
using LinearAlgebra: I, tr
using Test

# Every algebraic identity below is checked against the dense matrices produced by
# `tomatrix`, which is the definition of the representation.

"All `2Q`-bit strings on `Q` qubits."
allstrings(Q) = UInt.(0:(4^Q-1))

dense(s::PauliSentence) = tomatrix(s)

@testset "phase helpers" begin
    for k in 0:16
        @test _ipow(k) == im^k
        @test _mpow(k) == (-1)^k
    end
    for Q in 1:3, p in allstrings(Q), q in allstrings(Q)
        @test _prodsign(p, q, Q) == _symplectic_prod(p, q, Q).second
        @test _symplectic_prod(p, q, Q).first == p ⊻ q
        # The sign is the whole phase of a product of *bare* strings.
        @test tomatrix(p, Q) * tomatrix(q, Q) ≈ _prodsign(p, q, Q) * tomatrix(p ⊻ q, Q)
    end
end

@testset "scalar multiplication" begin
    @test (-1 * UPauli("Y-")).sign == -im
    @test (im * UPauli("Y-")).sign == -1
    @test UPauli("Y-") * im == im * UPauli("Y-")
    @test (im * Pauli("Z-", -1)).sign == -im
    @test Pauli("Z-", -1) * im == im * Pauli("Z-", -1)
    for c in (1, -1, im, -im)
        @test tomatrix(c * UPauli("XY")) ≈ c * tomatrix(UPauli("XY"))
        @test tomatrix(c * Pauli("XY", -1)) ≈ c * tomatrix(Pauli("XY", -1))
    end
    # A `Pauli` only holds a fourth root of unity.
    @test_throws ArgumentError 0.5 * UPauli("Z-")
    @test_throws ArgumentError 2 * Pauli("Z-")
end

@testset "products and commutation" begin
    for Q in 1:2, p in allstrings(Q), q in allstrings(Q)
        up, uq = UPauli(p, Q), UPauli(q, Q)
        expected = tomatrix(up) * tomatrix(uq)
        @test tomatrix(up * uq) ≈ expected
        @test tomatrix(Pauli(up) * uq) ≈ expected
        @test tomatrix(up * Pauli(uq)) ≈ expected
        @test tomatrix(Pauli(up) * Pauli(uq)) ≈ expected
        @test (up * uq).string == (up ⊻ uq).string
        @test unsigned_prod(up, uq) === up ⊻ uq
        # `com` reports whether the two commute, alongside their product string.
        @test com(up, uq).second ==
              isapprox(tomatrix(up) * tomatrix(uq), tomatrix(uq) * tomatrix(up))
        @test com(p, q, Q).second == com(up, uq).second
        @test com(p, q, Q).first == (up ⊻ uq).string
        @test (up^uq) == com(up, uq)
        @test com(up, uq).first === up ⊻ uq
    end
    # Squaring gives the identity, phase and all.
    for Q in 1:3, s in allstrings(Q)
        p = UPauli(s, Q)
        @test (p * p).string == 0
        @test (p * p).sign == 1
    end
end

@testset "cis of a single string" begin
    Q = 2
    for s in allstrings(Q), θ in (0.0, 0.37, pi / 4, -1.2)
        p = UPauli(s, Q)
        u = dense(cis(p, θ))
        @test u ≈ exp(im * θ * tomatrix(p))
        @test u * adjoint(u) ≈ I(2^Q)
    end
    @test length(cis(UPauli("Z-"), 0.0)) == 2
    @test tostring(cis(UPauli("Z-"), pi / 2))["Z-"] ≈ im
end

@testset "sentence arithmetic" begin
    Q = 3
    a = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 9, 26], ComplexF64[0.3, -0.7im, 0.2])
    b = PauliSentence{UInt,ComplexF64,Q}(UInt[9, 5], ComplexF64[1.5, 0.25])
    @test dense(a + b) ≈ dense(a) + dense(b)
    @test dense(a - b) ≈ dense(a) - dense(b)
    @test dense(a + b + a) ≈ 2 * dense(a) + dense(b)
    @test dense(a + b + a + b) ≈ 2 * dense(a) + 2 * dense(b)
    @test dense(2.5 * a) ≈ 2.5 * dense(a)
    @test dense(a * 2.5) ≈ 2.5 * dense(a)
    @test dense((1 + 2im) * a) ≈ (1 + 2im) * dense(a)
    @test dense(a * b) ≈ dense(a) * dense(b)
    @test dense(b * a) ≈ dense(b) * dense(a)
    @test dense(a * b - b * a) ≈ dense(a) * dense(b) - dense(b) * dense(a)
    @test dense(a^3) ≈ dense(a)^3
    @test dense(a^0) ≈ I(2^Q)
    negative = -1
    @test_throws ArgumentError a^negative
    # Coefficients live in the bare basis, where the identity carries the whole trace.
    @test tr(dense(a * a)) ≈ 2^Q * (a*a)[UInt(0)]
    @test trace(a * a) ≈ tr(dense(a * a))

    # Types promote rather than being pinned to whatever came first.
    r = PauliSentence{UInt16,Float64,Q}(UInt16[1], [1.0])
    @test valtype(r + a) == ComplexF64
    @test keytype(r + a) == UInt
    @test valtype(2im * r) == ComplexF64
    @test dense(r + a) ≈ dense(r) + dense(a)
    @test dense(r - a) ≈ dense(r) - dense(a)
    # Cancelling terms are kept as explicit zeros; arithmetic does not prune.
    @test length(a - a) == length(a)
    @test all(iszero, values(a - a))
end

@testset "cis of a sentence" begin
    Q = 2
    s = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 6], ComplexF64[0.3, 0.8])
    @test dense(cis(s, 0.7)) ≈ cis(0.7 * dense(s))
end

@testset "pruning" begin
    Q = 2
    s = PauliSentence{UInt,ComplexF64,Q}(UInt[1, 6, 3], ComplexF64[0.0, 1e-9, 1.0])
    @test _prunable(0.0, 0) && !_prunable(1e-9, 0)
    @test _prunable(1e-9, 1e-6) && !_prunable(1.0, 1e-6)
    @test length(_prune!(copy(s), 0)) == 2
    @test length(_prune!(copy(s), 1e-6)) == 1
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
        @test dense(ad!(copy(s), g, cos(2θ), sin(2θ))) ≈ expected
    end
    # `ad` leaves its argument alone; `ad!` does not.
    before = copy(s)
    ad(s, UPauli(UInt(5), Q), 0.4)
    @test s == before
    @test ad!(copy(s), UPauli(UInt(5), Q), 0.4) != before
    # The identity generates a global phase, which conjugation cannot see.
    @test dense(ad(s, UPauli(UInt(0), Q), 0.8)) ≈ dense(s)
    # A term is left alone by a rotation about itself.
    single = PauliSentence{UInt,ComplexF64,Q}(UInt[9], ComplexF64[1.0])
    @test dense(ad(single, UPauli(UInt(9), Q), 0.8)) ≈ dense(single)

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
    @test dense(ad!(copy(s), gens, cos.(2angles), sin.(2angles))) ≈ expected
    @test_throws DimensionMismatch ad(s, gens, [0.1, 0.2])
    @test_throws DimensionMismatch ad!(copy(s), gens, [0.1, 0.2])
    @test dense(ad(s, PauliList{UInt,Q}(), Float64[])) ≈ dense(s)

    # Undoing a sequence of rotations means walking it backwards: the generators do not
    # commute, so the order has to be reversed along with the angles.
    @test dense(ad(ad(s, gens, angles), reverse(gens), -reverse(angles))) ≈ dense(s)
    @test all(v -> abs(v) > 0.5, values(ad(s, gens, angles, atol=0.5)))
    @test all(v -> abs(v) > 0.5, values(ad!(copy(s), gens, angles, atol=0.5)))
    # A real-coefficient sentence is widened rather than rejected.
    real_s = PauliSentence{UInt,Float64,Q}(UInt[1], [1.0])
    @test dense(ad(real_s, UPauli(UInt(9), Q), 0.3)) ≈
          dense(cis(UPauli(UInt(9), Q), 0.3)) *
          dense(real_s) *
          adjoint(dense(cis(UPauli(UInt(9), Q), 0.3)))
end
