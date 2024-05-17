"""
This file contains the functions to generate the ecosystems for the HOI models.

"""


"""
Generates the community matrix components: R for the internal growth rate, D for the death rate
A for is the matrix with the base pairwise interactions. B for the 3rd order HOI
"""
function generate_base_parameters(; N::Int=20, µₐ::Real=-0.0, σₐ::Real=0.0, Γ::Real=0.0,
    µᵣ::Real=1, σᵣ::Real=0.0, µ_d::Real=-1, σ_d::Real=0.0)

    n = Array{Float64}(undef, N)
    m = Array{Float64}(undef, N, N)
    R = Array{Float64}(undef, N)
    d = Array{Float64}(undef, N)
    A = Array{Float64}(undef, N, N)

    generate_base_parameters!(n, m, A, R, d, N=N, µₐ=µₐ,
        σₐ=σₐ, Γ=Γ, µᵣ=µᵣ, σᵣ=σᵣ, µ_d=µ_d, σ_d=σ_d)
    return n, m, A, R, d
end

"""
Generates the community matrix components: R for the internal growth rate, D for the death rate
A for is the matrix with the base pairwise interactions. B for the 3rd order HOI.
Does everything in place, to reduce garbage collectrion time
"""
function generate_base_parameters!(n::Vector{Float64}, m::Matrix{Float64}, A::Matrix{<:Real}, R, d; N::Int=20, µₐ::Real=-0.0, σₐ::Real=0.0, Γ::Real=0.0,
    µᵣ::Real=1, σᵣ::Real=0.0, µ_d::Real=-1, σ_d::Real=0.0)

    if µₐ != 1
        n .= abs(µᵣ / (µ_d + µₐ))
    else
        n .= abs(µᵣ)
    end

    m .= 1.0

    # Pairwise interaction matrix
    if Γ != 0 && σₐ != 0
        pairwise_interaction_matrix_normal!(A; μ=μₐ, σ=σₐ, Γ=Γ)
    else
        rand!(Normal(µₐ / (N - 1), σₐ / sqrt(N - 1)), A)
    end

    A .= round.(A, digits=2)

    rand!(Normal(µᵣ, σᵣ), R)
    rand!(Normal(µ_d, σ_d), d)

    for i in 1:N
        A[i, i] = 0
    end
    return n, m, A, R, d
end
"""
Generates the community matrix components: R for the internal growth rate, D for the death rate
A for is the matrix with the base pairwise interactions. B for the 3rd order HOI
"""
function generate_full_HOI_parameters(; N::Int=20, µₐ::Real=-0.0, σₐ::Real=0.0, Γ::Real=0.0,
    µᵦ::Real=0.0, σᵦ::Real=0.0, µᵣ::Real=10^5, σᵣ::Real=0.0, µ_d::Real=-1, σ_d::Real=0.0)

    n = Array{Float64}(undef, N)
    m = Array{Float64}(undef, N, N)
    R = Array{Float64}(undef, N)
    d = Array{Float64}(undef, N)
    A = Array{Float64}(undef, N, N)
    B = Array{Float64}(undef, N, N, N)

    generate_full_HOI_parameters!(n, m, A, B, R, d, N=N, µₐ=µₐ,
        σₐ=σₐ, Γ=Γ, µᵦ=µᵦ, σᵦ=σᵦ, µᵣ=µᵣ, σᵣ=σᵣ, µ_d=µ_d, σ_d=σ_d)
    return n, m, A, B, R, d
end

"""
Generates the community matrix components: R for the internal growth rate, D for the death rate
A for is the matrix with the base pairwise interactions. B for the 3rd order HOI.
Does everything in place, to reduce garbage collectrion time
"""
function generate_full_HOI_parameters!(n, m, A::Matrix{<:Real}, B, R, d; N::Int=20, µₐ::Real=-0.0, σₐ::Real=0.0, Γ::Real=0.0,
    µᵦ::Real=0.0, σᵦ::Real=0.0, µᵣ::Real=10^5, σᵣ::Real=0.0, µ_d::Real=-1, σ_d::Real=0.0)

    generate_base_parameters!(n, m, A, R, d; N=N, µₐ=µₐ, σₐ=σₐ, Γ=Γ, µᵣ=µᵣ, σᵣ=σᵣ, µ_d=µ_d, σ_d=σ_d)
    rand!(Normal(µᵦ / N, σᵦ / sqrt(N)), B)
    B .= (1 / µᵣ) .* B
    for i in 1:N
        B[i, i, :] .= 0
        B[i, :, i] .= 0
        B[:, i, i] .= 0
    end
    return n, m, A, B, R, d
end


function generate_intransitive_parameters(; type::Int64=1)
    N = 4
    n = Array{Float64}(undef, N)
    m = Array{Float64}(undef, N, N)
    R = Array{Float64}(undef, N)
    d = Array{Float64}(undef, N)
    A = Array{Float64}(undef, N, N)

    generate_intransitive_parameters!(n, m, A, R, d; type=type)
    return n, m, A, R, d
end

function generate_intransitive_parameters!(n::Vector{<:Real}, m, A::Matrix{Float64}, R, d; type::Int64=1)
    N = 4
    n .= rand(0.5:0.1:1.5, N)
    m .= 1.0
    R .= 1.0
    d .= -1.0
    intransitive_matrix!(A; type)
end


function generate_3_cycle(; intransitive=true, superiority=(1, 2, 3))
    N = 3
    n = [0.7, 0.6, 1.3]
    m = ones(N, N)
    R = ones(N)
    d = -ones(N)
    if intransitive
        A = [0 1.0 -1.0; -1.0 0 1.0; 1.0 -1.0 0]
    else
        A = zeros(N, N)
        for i in 1:length(superiority)
            for j in i:length(superiority)
                if i < j
                    A[superiority[i], superiority[j]] = 1
                    A[superiority[j], superiority[i]] = -1
                end
            end
        end
    end
    return n, m, A, R, d
end

function generate_nonuniform_3_cycle(αs; intransitive=true, superiority=(1, 2, 3))
    """
    generates n, m, A, R and d for a 3-cycle with non-uniform α values
    αs is a vector of length 3, containing the α values for respectively α_12, α_13, α_23
    """
    N = 3
    n = [0.7, 0.6, 1.3]
    m = ones(N, N)
    R = ones(N)
    d = -ones(N)
    if intransitive
        A = [0 1.0 -1.0; -1.0 0 1.0; 1.0 -1.0 0]
    else
        A = zeros(N, N)
        for i in 1:length(superiority)
            for j in i:length(superiority)
                if i < j
                    A[superiority[i], superiority[j]] = 1
                    A[superiority[j], superiority[i]] = -1
                end
            end
        end
    end
    A[1, 2] *= αs[1]
    A[2, 1] *= αs[1]
    A[1, 3] *= αs[2]
    A[3, 1] *= αs[2]
    A[2, 3] *= αs[3]
    A[3, 2] *= αs[3]
    return n, m, A, R, d
end