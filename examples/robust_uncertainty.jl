#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# robust_uncertainty.jl
#
# Computes the Value at Risk 
#############################################################################


using JuMP

const R = 1
const d = 1
const ẟ = 0.05
const ɛ = 0.05
const N = ceil((2+2log(2/ẟ))^2) + 1

Γ1(ẟ,N) = (R/sqrt(N))*(2+sqrt(2*log(1/ẟ)))
Γ2(ẟ,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/ẟ)))

μhat = randn(d)
Σhat = randn(d,d)

m = Model()

@defSDPVar(m, Σ[d])
@defVar(m, u[1:d])
@defVar(m, μ[1:d])

# addConstraint(m, norm2 <= Γ1(ẟ/2,N))
norm2 = Γ1(ẟ/2,N)*eye(d+1,d+1)
for i in 1:d
    norm2 += (μhat[i] - μ[i]) * sparse([i,d],[d,i],[1.0,1.0],d+1,d+1)
end

addConstraint(m, Γ1(ẟ/2,N)*eye(d+1,d+1) +  >= 0)
addConstraint(m, normF <= Γ2(ẟ/2,N))

A = [(1-ɛ)/ɛ (u-μ)'; (u-μ) Σ]
addConstraint(m, A >= 0)
