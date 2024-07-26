


project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))

const adverbs= unique(lowercase.(open(pkgdir(JLDispense,"src/RandomProtocolNames/adverbs.txt")) do f readlines(f) end ))
const verbs=unique(lowercase.(open(pkgdir(JLDispense,"src/RandomProtocolNames/verbs_present_participle.txt")) do f readlines(f) end ))
const max_n = 0.1*length(adverbs)* length(verbs)

function random_protocol_name()
    return string(rand(adverbs),"_",rand(verbs))
end 

function random_protocol_name(n::Integer) 
    n < max_n ? nothing : error("Cannot generate $n random protocol names at once. Entry must be less than $max_n")
    out=String[]
    while length(out) < n  
        x=random_protocol_name()
        if !in(x,out)
            push!(out,x)
        end 
    end


    return out 
end 






