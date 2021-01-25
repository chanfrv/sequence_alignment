#!/usr/bin/julia

module CmdLine

using ArgParse


# enums
@enum SequenceEnum Any ADN ARN PROT
@enum CommandEnum score align

# types
GammaType = Tuple{Int,Int}
SequenceType = Tuple{SequenceEnum,AbstractString}


# Gamma parsing
function ArgParse.parse_item(::Type{GammaType}, x::AbstractString)
    return GammaType(map(x -> parse(Int, x), split(x, ",")))
end

# Command parsing
function ArgParse.parse_item(::Type{CommandEnum}, x::AbstractString)
    if lowercase(x) == "score"
        return score
    elseif lowercase(x) == "align"
        return align
    else
        throw(error("Not a valid command!"))
    end
end


function ismatch(regex, str)
    return match(regex, str).match == str
end


# Sequence parsing
function ArgParse.parse_item(::Type{SequenceType}, x::AbstractString)
    if ismatch(r"[AGC]+", x)
        return (Any, x) # will depends on the other
    elseif ismatch(r"[ATGC]+", x)
        return (ADN, x) # ADN
    elseif ismatch(r"[AUGC]+", x)
        return (ARN, x) # ARN
    elseif ismatch(r"[ARNDCQEGHILKMFPSTWYVBZX]+", x)
        return (PROT, x)
    else
        throw(error("Not a valid sequence!"))
    end
end

# Command line parsing
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--gamma"
            help = "Paramètres de la fonction d'évaluation γ(x) = ex + o"
            arg_type = GammaType
            default = (-1,0)
            range_tester = x -> all(x -> x <= 0, x)
        "cmd"
            help = "Commande (score ou align)"
            arg_type = CommandEnum
            required = true
        "x"
            help = "Séquence 1"
            arg_type = SequenceType
            required = true
        "y"
            help = "Séquence 2"
            arg_type = SequenceType
            required = true
    end

    return parse_args(s, as_symbols=true)
end

end
