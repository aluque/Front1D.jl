function pretty_print(io::IO, x, level=0)
    @show nfields(x)
    for i in 1:nfields(x)
        k = fieldname(typeof(x), i)
        v = getfield(x, k)
        
        if nfields(v) > 0
            print(io, join(fill(" ", level * 8)))
            printstyled(io, String(k), color=:light_yellow, bold=true)            
            printstyled(io, " => \n", color=:light_black)            
            pretty_print(io, v, level + 1)
        else
            print(io, join(fill(" ", level * 8)))
            printstyled(io, String(k), color=:light_yellow, bold=true)
            printstyled(io, " => ", color=:light_black)
            printstyled(io, repr(v) * "\n", color=:blue) 
        end
    end
    nothing
end

pretty_print(x) = pretty_print(stdout, x, kw...)
