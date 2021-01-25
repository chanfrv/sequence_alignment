#!/usr/bin/julia

include("cmdline.jl")
include("sigma.jl")


# Returns Sij
function s(i, j, Σ, Σ_idx, γ, x, y)
    (e, o) = γ

    if i == 0 && j == 0
        return (0.0, nothing)
    elseif i == 0
        return (j * e, (0, j-1))
    elseif j == 0
        return (i * e, (i-1, 0))
    else
        (sub, _) = s(i-1, j-1, Σ, Σ_idx, γ, x, y)
        (ins, _) = s(i-1, j, Σ, Σ_idx, γ, x, y)
        (del, _) = s(i, j-1, Σ, Σ_idx, γ, x, y)

        xi = Σ_idx[x[i]]
        yj = Σ_idx[y[j]]

        return max((Σ[xi, yj] + sub, (i-1, j-1)),
                   (e + ins, (i-1, j)),
                   (e + del, (i, j-1)))
    end
end

function backtrack(S, x, y)
    (i, j) = size(S)
    (i, j) = (i-1, j-1) # we go from 0 to size-1

    (_, kl) = S[end, end]

    if kl == nothing
        return ("", "")

    else
        if kl == (i-1, j-1)
            (x′, y′) = backtrack(S[1:end-1, 1:end-1], x[1:end-1], y[1:end-1])
            return (x′ * x[end], y′ * y[end])

        elseif kl == (i-1, j)
            (x′, y′) = backtrack(S[1:end-1, :], x[1:end-1], y[:])
            return (x′ * x[end], y′ * '-')

        elseif kl == (i, j-1)
            (x′, y′) = backtrack(S[:, 1:end-1], x[:], y[1:end-1])
            return (x′ * '-', y′ * y[end])

        else
            println("(k,l) = ", kl, " (i,j) = ", (i, j))
            throw(error("Wrong backtrack matrix"))
        end
    end
end

 
function main()
    # parse args
    parsed_args = CmdLine.parse_commandline()

    γ = parsed_args[:gamma]
    cmd = parsed_args[:cmd]
    (xtype, x) = parsed_args[:x]
    (ytype, y) = parsed_args[:y]

    # deal with ambiguities
    if xtype == CmdLine.Any && ytype == CmdLine.Any
        (xtype, ytype) = (CmdLine.ADN, CmdLine.ADN)
    elseif xtype == CmdLine.Any
        xtype = ytype
    elseif ytype == CmdLine.Any
        ytype = xtype
    end

    # check same type of sequence
    if xtype != ytype
        throw(error("x and y are different sequence types!"))
    end

    # get Σ
    (Σ, Σ_idx) = if xtype == CmdLine.ADN
        (Σ_nuc, Σ_adn_idx)
    elseif xtype == CmdLine.ARN
        (Σ_nuc, Σ_arn_idx)
    else
        (Σ_prot, Σ_prot_idx)
    end

    M = length(x)
    N = length(y)

    # get S
    S = [ s(i, j, Σ, Σ_idx, γ, x, y) for i=0:M, j=0:N ]

    # result
    if cmd == CmdLine.score
        (score, _) = S[end, end]
        println(score)
    elseif cmd == CmdLine.align
        (x′, y′) = backtrack(S, x, y)
        println(x′)
        println(y′)
    end
end

main()
