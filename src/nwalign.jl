#!/usr/bin/julia

include("cmdline.jl")
include("sigma.jl")


function get_σ(xtype, x, ytype, y)
    # deal with ambiguities
    if xtype == CmdLine.Any && ytype == CmdLine.Any
        (xtype, ytype) = (CmdLine.ADN, CmdLine.ADN)

    elseif xtype == CmdLine.Any || (xtype == CmdLine.ADN && ytype == CmdLine.PROT)
        xtype = ytype

    elseif ytype == CmdLine.Any || (xtype == CmdLine.PROT && ytype == CmdLine.ADN)
        ytype = xtype
    end

    # check same type of sequence
    if xtype != ytype
        throw(error("x and y are different sequence types!"))
    end

    # get Σ
    Σ = if xtype == CmdLine.ADN
            Σ_adn_dict
        elseif xtype == CmdLine.ARN
            Σ_arn_dict
        elseif xtype == CmdLine.PROT
            Σ_prot_dict
        else
            throw(error("Wrong sequence type"))
        end

    return (i,j) -> Σ[ (x[i], y[j]) ]

end


function Needleman_Wunsch_rec(i, j, S, B, σ, γ, M, N)

    if i > M || j > N
        return

    elseif i == 0 && j == 0
        S[i+1,j+1] = γ(0)
        B[i+1,j+1] = nothing

    elseif i == 0
        S[i+1,j+1] = γ(j)
        B[i+1,j+1] = (1,j)

    elseif j == 0
        S[i+1,j+1] = γ(i)
        B[i+1,j+1] = (i,1)

    else
        sub = (σ(i,j) + S[i-1+1,j-1+1], (i-1+1,j-1+1))
        ins = [ (γ(k) + S[i-k+1,j+1  ], (i-k+1,j+1  )) for k=1:i ]
        del = [ (γ(k) + S[i+1  ,j-k+1], (i+1  ,j-k+1)) for k=1:j ]

        (score, back) = max(sub, maximum(ins), maximum(del))

        S[i+1,j+1] = score
        B[i+1,j+1] = back
    end

    Needleman_Wunsch_rec(i  , j+1, S, B, σ, γ, M, N)
    Needleman_Wunsch_rec(i+1, j  , S, B, σ, γ, M, N)
    Needleman_Wunsch_rec(i+1, j+1, S, B, σ, γ, M, N)
end


# Returns Sij
function Needleman_Wunsch(σ, γ, M, N)
    S = Int[ 0 for i=0:M, j=0:N ]
    B = Union{Nothing, Tuple{Int,Int}}[ nothing for i=0:M, j=0:N ]

    Needleman_Wunsch_rec(0, 0, S, B, σ, γ, M, N)

    return (S, B)
end


function backtrack(B, x, y)
    (u, v) = size(B)
    back = B[end, end]

    if back == nothing
        return ("", "")
    end

    (bk, bl) = back

    if back == (u-1, v-1)
        (x′, y′) = backtrack(B[1:end-1, 1:end-1], x, y)
        return (x′ * x[u-1], y′ * y[v-1])

    elseif bl == v
        k = u - bk
        (x′, y′) = backtrack(B[1:u-k, :], x, y)
        return (x′ * x[u-k:u-1], y′ * ('-'^k))

    elseif bk == u
        k = v - bl
        (x′, y′) = backtrack(B[:, 1:v-k], x, y)
        return (x′ * ('-'^k), y′ * y[v-k:v-1])

    else
        throw(error("Wrong backtrack matrix"))
    end
end


function main()
    # parse args
    parsed_args = CmdLine.parse_commandline()

    (e, o) = parsed_args[:gamma]
    cmd = parsed_args[:cmd]
    (xtype, x) = parsed_args[:x]
    (ytype, y) = parsed_args[:y]

    γ = x -> e*x + o
    σ = get_σ(xtype, x, ytype, y)

    (M, N) = (length(x), length(y))

    # get S
    (S, B) = Needleman_Wunsch(σ, γ, M, N)

    # result
    if cmd == CmdLine.score
        score = S[end, end]
        println(score)
    elseif cmd == CmdLine.align
        (x′, y′) = backtrack(B, x, y)
        println(x′)
        println(y′)
    end
end

main()
