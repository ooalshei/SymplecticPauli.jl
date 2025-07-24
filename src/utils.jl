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

counti(string::Unsigned, Q::Integer) = Q - countx(string, Q) - county(string, Q) - countz(string, Q)
counti(p::AbstractPauli) = counti(p.string, p.qubits)

function tostring(p::UPauli)::String
    result = ""
    string = digits(p.string, base=2, pad=2*p.qubits)
    for i in 1:p.qubits
        if string[i] == string[i + p.qubits] == 1
            result *= "Y"
        elseif string[i] == 1
            result *= "X"
        elseif string[i + p.qubits] == 1
            result *= "Z"
        else
            result *= "-"
        end
    end
    return result
end
function tostring(p::Pauli)::String
    sign = p.sign * C8(-im)^county(p)
    sign == 1 && return "(+)" * tostring(UPauli(p))
    sign == -1 && return "(-)" * tostring(UPauli(p))
    sign == im && return "(i)" * tostring(UPauli(p))
    sign == -im && return "(-i)" * tostring(UPauli(p))
end
function tostring(s::PauliSentence)
    result = Dict{String,valtype(s)}()
    for (key, value) in pairs(s)
        result[tostring(UPauli(UInt(key), s.qubits))] = (-im)^county(key, s.qubits) * value
    end
    return result
end
tostring(v::PauliList) = tostring.(UPauli.(v.strings, v.qubits))

function tomatrix(string::Unsigned, Q::Integer)
    _check_string_length(string, Q)
    result = I(1)
    string = digits(string, base=2, pad=2*Q)
    for i in 1:Q
        if string[i] == string[i + Q] == 1
            result = result ⊗ σ₂real  # σy
        elseif string[i] == 1
            result = result ⊗ σ₁  # σx
        elseif string[i + Q] == 1
            result = result ⊗ σ₃  # σz
        else
            result = result ⊗ I(2)  # I
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