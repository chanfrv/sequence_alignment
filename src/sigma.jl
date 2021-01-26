#!/usr/bin/julia


const Σ_adn_idx = Dict('A' => 1,
                       'T' => 2,
                       'G' => 3,
                       'C' => 4)

const Σ_arn_idx = Dict('A' => 1,
                       'U' => 2,
                       'G' => 3,
                       'C' => 4)

const Σ_nuc = [ (i == j ? 1 : -1) for i=1:4, j=1:4 ]

const Σ_adn_dict = Dict([ ((kx, ky), Σ_nuc[vx,vy]) for (kx, vx) in Σ_adn_idx, (ky, vy) in Σ_adn_idx ])
const Σ_arn_dict = Dict([ ((kx, ky), Σ_nuc[vx,vy]) for (kx, vx) in Σ_arn_idx, (ky, vy) in Σ_arn_idx ])


const Σ_prot_idx = Dict('A' => 1,
                        'R' => 2,
                        'N' => 3,
                        'D' => 4,
                        'C' => 5,
                        'Q' => 6,
                        'E' => 7,
                        'G' => 8,
                        'H' => 9,
                        'I' => 10,
                        'L' => 11,
                        'K' => 12,
                        'M' => 13,
                        'F' => 14,
                        'P' => 15,
                        'S' => 16,
                        'T' => 17,
                        'W' => 18,
                        'Y' => 19,
                        'V' => 20,
                        'B' => 21,
                        'Z' => 22,
                        'X' => 23)

#                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
const Σ_prot = [ 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 ;
                -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 ; 
                -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 ; 
                -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 ; 
                 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 ; 
                -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 ; 
                -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 ; 
                 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 ; 
                -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 ; 
                -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 ;
                -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 ;
                -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 ; 
                -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 ; 
                -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 ;
                -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 ;
                 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 ; 
                 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 ;
                -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 ;
                -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 ; 
                 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 ;
                -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 ; 
                -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 ; 
                 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 ]

const Σ_prot_dict = Dict([ ((kx, ky), Σ_prot[vx,vy]) for (kx, vx) in Σ_prot_idx, (ky, vy) in Σ_prot_idx ])

