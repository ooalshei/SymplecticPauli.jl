_check_string_length(string::Unsigned, Q::Integer) = string > 4^Q - 1 ? throw(ArgumentError("String must not exceed $(4^Q - 1).")) : nothing
function _check_type(::Type{T}, s::Integer) where {T<:Unsigned}
    Base.hastypemax(T) && typemax(T) < 4^s - 1 && throw(ArgumentError("String length cannot exceed $(count_ones(typemax(T)) ÷ 2)). Consider using a larger unsigned integer type."))
    nothing
end

countx(string::Unsigned, Q::Integer) = (_check_string_length(string, Q); count_ones(string & (2^Q - 1)))
countx(p::AbstractPauli) = countx(p.string, p.qubits)

county(string::Unsigned, Q::Integer) = (_check_string_length(string, Q); count_ones(string & (string >> Q)))
county(p::AbstractPauli) = county(p.string, p.qubits)

countz(string::Unsigned, Q::Integer) = (_check_string_length(string, Q); count_ones(string >> Q))
countz(p::AbstractPauli) = countz(p.string, p.qubits)

function tomatrix(string::Unsigned, Q::Integer)
    _check_string_length(string, Q)
    result = I(1)
    string = digits(string, base=2, pad=2*Q)
    for i in 1:Q
        if string[i] == string[i + Q] == 1
            result = σy_real ⊗ result  # σy
        elseif string[i] == 1
            result = σx ⊗ result  # σx
        elseif string[i + Q] == 1
            result = σz ⊗ result  # σz
        else
            result = I(2) ⊗ result  # I
        end
    end
    return result
end
tomatrix(p::AbstractPauli) = Pauli(p).sign * tomatrix(p.string, p.qubits)

function tomatrix(p::PauliSentence)
    result = zeros(ComplexF64, 2^p.qubits, 2^p.qubits)
    for (key, value) in p
        result .+= value .* tomatrix(key, p.qubits)
    end
    return result
end