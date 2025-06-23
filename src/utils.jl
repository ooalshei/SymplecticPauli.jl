countx(string::Unsigned, Q::Integer) = count_ones(string & (2^Q - 1))
countx(p::AbstractPauli) = countx(p.string, p.qubits)

county(string::Unsigned, Q::Integer) = count_ones(string & (string >> Q))
county(p::AbstractPauli) = county(p.string, p.qubits)

countz(string::Unsigned, Q::Integer) = count_ones(string >> Q)
countz(p::AbstractPauli) = countz(p.string, p.qubits)
