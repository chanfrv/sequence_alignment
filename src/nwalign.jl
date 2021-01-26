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

# Fill a cell of S and B according to the Needleman-Wunsch formula
function Needleman_Wunsch_cell(i, j, S, B, σ, γ, M, N)
    sub = (σ(i,j) + S[i-1+1,j-1+1], (i-1+1,j-1+1))
    ins = [ (γ(k) + S[i-k+1,j+1  ], (i-k+1,j+1  )) for k=1:i ]
    del = [ (γ(k) + S[i+1  ,j-k+1], (i+1  ,j-k+1)) for k=1:j ]

    (score, back) = max(sub, maximum(ins), maximum(del))

    S[i+1,j+1] = score
    B[i+1,j+1] = back
end

# Start from a cell Sᵢⱼ and fill both rows i and cols j
function Needleman_Wunsch_rec(zi, zj, S, B, σ, γ, M, N)

    if zi > M || zj > N
        return
    end

    # for each row including (zi,zj)
    for i = zi:M
        Needleman_Wunsch_cell(i, zj, S, B, σ, γ, M, N)
    end

    # for each col excluding (zi,zj) to avoid useless call
    for j = zj+1:N
        Needleman_Wunsch_cell(zi, j, S, B, σ, γ, M, N)
    end

    Needleman_Wunsch_rec(zi+1, zj+1, S, B, σ, γ, M, N)
end


# Returns Sij
function Needleman_Wunsch(σ, γ, M, N)
    S = Int[ 0 for i=0:M, j=0:N ]
    B = Union{Nothing, Tuple{Int,Int}}[ nothing for i=0:M, j=0:N ]

    # Fill the borders
    S[1,1] = γ(0)
    B[1,1] = nothing

    for i = 1:M
        S[i+1,1] = γ(i)
        B[i+1,1] = (i,1)
    end

    for j = 1:N
        S[1,j+1] = γ(j)
        B[1,j+1] = (1,j)
    end

    # Start the recursive algorithm
    Needleman_Wunsch_rec(1, 1, S, B, σ, γ, M, N)

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
