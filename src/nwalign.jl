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


# Returns Sij
function Needleman_Wunsch(σ, γ, M, N)
    S = Int[ 0 for i=0:M, j=0:N ]
    B = Union{Nothing, Tuple{Int,Int}}[ nothing for i=0:M, j=0:N ]

    s₀ = 1

    # Fill the borders
    S[s₀,s₀] = γ(0)
    B[s₀,s₀] = nothing

    for i = 1:M
        sᵢ = i+1

        S[sᵢ,s₀] = γ(i)
        B[sᵢ,s₀] = (sᵢ-1,s₀)
    end

    for j = 1:N
        sⱼ = j+1

        S[s₀,sⱼ] = γ(j)
        B[s₀,sⱼ] = (s₀,sⱼ-1)
    end

    for i = 1:M, j = 1:N
        (sᵢ,sⱼ) = (i+1,j+1)
 
        sub = (σ(i,j) + S[sᵢ-1,sⱼ-1], (sᵢ-1,sⱼ-1))
        ins = [ (γ(k) + S[sᵢ-k,sⱼ  ], (sᵢ-k,sⱼ  )) for k=1:i ]
        del = [ (γ(k) + S[sᵢ  ,sⱼ-k], (sᵢ  ,sⱼ-k)) for k=1:j ]

        (score, back) = max(sub, maximum(ins), maximum(del))

        S[sᵢ,sⱼ] = score
        B[sᵢ,sⱼ] = back
    end

    display(S)
    println()

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

    (e, o)     = parsed_args[:gamma]
    cmd        = parsed_args[:cmd]
    (xtype, x) = parsed_args[:x]
    (ytype, y) = parsed_args[:y]

    γ = x -> e*x + o
    σ = get_σ(xtype, x, ytype, y)

    (M, N) = (length(x), length(y))

    # result
    if cmd == CmdLine.score
        # get S
        (S, _) = Needleman_Wunsch(σ, γ, M, N)
        # print score
        score = S[end, end]
        println(score)

    elseif cmd == CmdLine.align
        # get S
        (_, B) = Needleman_Wunsch(σ, γ, M, N)
        # backtrack
        (x′, y′) = backtrack(B, x, y)
        # print alignment
        println(x′)
        println(y′)
    end
end

main()
