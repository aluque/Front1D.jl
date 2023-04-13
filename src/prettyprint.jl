using REPL

function pretty_print(io::IO, x, level=0)
    
    for i in 1:nfields(x)
        k = fieldname(typeof(x), i)
        v = getfield(x, k)
        
        if nfields(v) > 0
            print(io, join(fill(" ", level * 4)))
            printstyled(io, String(k), color=:light_yellow, bold=true)            
            printstyled(io, " => \n", color=:light_black)            
            pretty_print(io, v, level + 1)
        else
            doc = sprint(show, REPL.fielddoc(typeof(x), k))
            print(io, join(fill(" ", level * 4)))
            printstyled(io, "# " * doc, color=:light_black)

            print(io, join(fill(" ", level * 4)))
            printstyled(io, String(k), color=:light_yellow, bold=true)
            printstyled(io, " => ", color=:light_black)
            printstyled(io, repr(v) * "\n\n", color=:blue) 
        end
    end
    nothing
end

pretty_print(x, kw...) = pretty_print(stdout, x, kw...)
