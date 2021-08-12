
"""
    NotImplementedError{M}(m)

`Exception` thrown when a method from a NetworkEpidemics interface
is not implemented by a given type. Borrowed from the `LightGraphs` source code.
"""
struct NotImplementedError{M} <: Exception
    m::M
    NotImplementedError(m::M) where {M} = new{M}(m)
end

Base.showerror(io::IO, ie::NotImplementedError) = print(io, "method $(ie.m) not implemented.")

_NI(m) = throw(NotImplementedError(m))