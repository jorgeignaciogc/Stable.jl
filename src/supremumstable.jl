using Distributions, SpecialFunctions

struct SupremumStable <: ContinuousUnivariateDistribution
  α::Float64
  β::Float64
  θ::Float64
  ρ::Float64
  SupremumStable(α::Real,β::Real) = ( β <- 1 || β > 1 || α <= 0 || α > 2 ) ?
    error("Parameters' requirements unmet:\n (α,β)∈(0,2]×[-1,1]") :
    α == 2 ? PositiveStable(2,-1) :
    β == -1 && α > 1 ? PositiveStable(α,β) :
    new(α,β,β*(α <= 1 ? 1 :
    (α-2)/α),(1 + β*(α <= 1 ? 1 : (α-2)/α))/2)
end

import Distributions.params
function params(d::SupremumStable)
  return (d.α,d.β,d.θ,d.ρ)
end

import Distributions.minimum
function minimum(d::SupremumStable)
  return (β == -1 && α <= 1) ? 0 : β == 1 && α == 1 ? 1. : 0.
end

import Distributions.maximum
function maximum(d::SupremumStable)
  return (β == -1 && α <= 1) ? 0 : β == 1 && α == 1 ? 1. : Inf
end

import Distributions.insupport
function insupport(d::SupremumStable,x::Real)
  return (β == -1 && α <= 1) ? x == 0 : x >= 0
end

import Distributions.mean
function mean(d::SupremumStable)
  return d.α*d.ρ*mean(PositiveStable(d.α,d.β))
end

struct PerfectSupremumStable <: Sampleable{Univariate,Continuous}
  α::Float64
  β::Float64
  θ::Float64
  ρ::Float64
  d::Float64
  η::Float64
  δ::Float64
  γ::Float64
  κ::Float64
  Δ::Int64
  mAst::Int64
  PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real,Δ::Int,mAst::Int) =
    (β < -1 || β > 1 || (β == -1 && α <= 1) || α <= 0 || α > 2) ?
    error("Parameters' requirements unmet: (α,β)∈(0,2]×[-1,1]-(0,1]×{-1}") : (0 >= γ || γ >= α) ?
    error("Parameters' requirements unmet: α>γ>0") : (0 >= δ || δ >= d || d >= 2/(α + β*(α <= 1 ? α : α-2))) ?
    error("Parameters' requirements unmet: 1/(αρ)>d>δ>0") : κ < 0 ?
    error("Parameters' requirements unmet: κ≥0") : Δ <= 0 ?
    error("Parameters' requirements unmet: Δ≥0") : mAst < 0 ?
    error("Parameters' requirements unmet: m*≥0") : new(α, β, β*(α <= 1 ? 1 : (α-2)/α),(1+β*(α <= 1 ? 1 : (α-2)/α))/2 , d,
      etaF(d*(α + β*(α <= 1 ? α : α-2))/2)*(α + β*(α <= 1 ? α : α-2))/2 , δ, γ, κ + Int(ceil(max(2/(α+β*(α <= 1 ? α : α-2)),
      log(2)*2/(3*etaF(d*(α+β*(α <= 1 ? α : α-2))/2)*(α+β*(α <= 1 ? α : α-2)) ) ))), Δ, mAst)
  PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real,Δ::Int) = PerfectSupremumStable(α,β,d,δ,γ,κ,Δ,12)
  PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real,κ::Real) = PerfectSupremumStable(α,β,d,δ,γ,κ,40)
  PerfectSupremumStable(α::Real,β::Real,d::Real,δ::Real,γ::Real) = PerfectSupremumStable(α,β,d,δ,γ,3)
  PerfectSupremumStable(α::Real,β::Real) = PerfectSupremumStable(α,β,(2/3)*2/(α+β*(α <= 1 ? α : α-2)),(1/3)*2/(α+β*(α <= 1 ? α : α-2)),.95*α)
end

import Distributions.sampler
function sampler(d::SupremumStable,dr::Real,δ::Real,γ::Real,κ::Real,r::Real,lag::Int64,mAst::Int64)
  return PerfectSupremumStable(d.α,d.β,dr,δ,γ,κ,r,lag,mAst)
end

function sampler(d::SupremumStable,dr::Real,δ::Real,γ::Real,κ::Real,r::Real,lag::Int64)
  return PerfectSupremumStable(d.α,d.β,dr,δ,γ,κ,r,lag)
end

function sampler(d::SupremumStable,dr::Real,δ::Real,γ::Real)
  return PerfectSupremumStable(d.α,d.β,dr,δ,γ)
end

function sampler(d::SupremumStable)
  return PerfectSupremumStable(d.α,d.β)
end

function params(d::PerfectSupremumStable)
  return (α,β,θ,ρ,d,η,δ,γ,κ,r,lag,mAst)
end

using LambertW
# @Input: normalised drift d'=αρd∈(0,1)
# @Returns: normalised η'=η/(αρ)>0
function etaF(d::Real)
  return -lambertw(-d*exp(-d),-1)/d-1
end

# @Input: drift d, η, 1/(αρ+η), and a bound M we want to know if we exceed
# @Returns: a sample of 1{R_0 <= M}
function reflectedProcess(d::Real,et::Real,et1::Real,M::Real)
  T = rand(Exponential())/et
  if T <= M
    return true
  end
  Sn = [0.]
  while Sn[end] < T
    push!(Sn,Sn[end]-rand(Exponential())*et1+d)
  end
  return maximum(Sn[1:(end-1)]) <= M
end

# @Input: drift d, η, 1/(αρ+η), a bound M that it should not exceed,
# and another bound M1 we want to know if we exceed
# @Returns: a sample of 1{R_0 <= M1} given R_0<M
function reflectedProcess(d::Real,et::Real,et1::Real,M::Real,M1::Real)
  m = exp(-M-d)
  Sn = [0.]
  while true
    T = -log(rand(Uniform(m,1)))/et # Exponential conditioned on <= M+d, since T>M+d always yields a rejection
    if T <= M1
      return true
    end
    while Sn[end] < T
      push!(Sn,Sn[end]-rand(Exponential())*et1+d)
    end
    max = maximum(Sn[1:(end-1)])
    if max <= M
      return max <= M1
    else
      Sn = [0.]
    end
  end
end

# @Input: drift d, η, 1/(αρ+η), κ and an upper bound x
# @Returns: A simulated run of the random walk up to the hitting time of the barrier
# -2κ conditioned on never exceeding x
function downRW(d::Real,rar::Real,et::Real,et1::Real,κ::Real,x::Real)
  Sn = [0.]
  INCn = [0.]
  if isfinite(x)
    while true
      while Sn[end] > -2 *κ
        push!(INCn,rand(Exponential())*rar)
        push!(Sn,Sn[end]-INCn[end]+d)
      end
      if maximum(Sn) < x && reflectedProcess(d,et,et1,x-Sn[end])
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.]
        INCn = [0.]
      end
    end
  else
    while Sn[end] > -2 *κ
      push!(INCn,rand(Exponential())*rar)
      push!(Sn,Sn[end]-INCn[end]+d)
    end
    return (Sn[2:end],INCn[2:end])
  end
end

# @Input: drift d, η, 1/(αρ+η), κ, an upper bound x that will be attained
# @Returns: A simulated run of the random walk up to the hitting time of the barrier
# κ conditioned on never exceeding x
function upRW(d::Real,et::Real,et1::Real,κ::Real,x::Real)
  Sn = [0.]
  INCn = [0.]
  if x == Inf
    while true
      while Sn[end] < κ
        push!(INCn,rand(Exponential())*et1)
        push!(Sn,Sn[end]-INCn[end]+d)
      end
      if rand(Uniform()) <= exp(-et*Sn[end])
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.]
        INCn = [0.]
      end
    end
  else
    while true
      while Sn[end] < κ
        push!(INCn,rand(Exponential())*et1)
        push!(Sn,Sn[end]-INCn[end]+d)
      end
      if rand(Uniform()) <= exp(-et*Sn[end]) && reflectedProcess(d,et,et1,x-Sn[end]) # The condition maximum(Sn) < x1 is trivially satisfied by construction of κ
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.]
        INCn = [0.]
      end
    end
  end
end

# Blanchet-Sigman algorithm
# @Input: drift d, η, 1/(αρ+η), κ and an upper bound x, last value of random walk s,
# @Returns: A simulated run of the random walk up to the first time we know an upper bound for the walk,
# the increments of the RW, the new upper bound, and the last time at which we may compute D_n
function BSAlgorithm(d::Real,et::Real,et1::Real,κ::Real,x::Real,s::Real)
  Sn = [Float64(s)]
  INCn = [0.]
  x1 = x
  t = 0

  while true
    # Step 1
    (Sn1,INCn1) = downRW(d,et,et1,κ,x1-Sn[end])
    # Step 2
    append!(Sn,Sn[end] .+ Sn1)
    append!(INCn,INCn1)
    t = t+length(Sn1)
    # Step 3 and 4
    if reflectedProcess(d,et,et1,x1-Sn[end],κ)
      x1 = Sn[end]+κ
      t = t-length(Sn1)
      break
    else
      (Sn1,INCn1) = upRW(d,et,et1,κ,x1-Sn[end])
      append!(Sn,Sn[end] .+ Sn1)
      append!(INCn,INCn1)
      t = t+length(Sn1)
    end
  end

  return (Sn[2:end],INCn[2:end],x1,t)
end

# Algorithm 3 - 1 (unconditional probabilities)
# @Input: S^+(α,ρ), δ, dega = δ*γ, dega2 = -1/(1-exp(-dega)), mom = E(S_1^γ), shift = m-k >= m*-1
function chiStable1(d::PositiveStable,δ::Real,dega::Real,dega2::Real,mom::Real,shift::Integer)
  U = rand(Uniform())
  n = shift
  S = Float64[]
  while true
    n = n+1
    aux = exp(δ*n)
    p = cdf(d,aux) # 1-pn
    aux1 = exp(-(n+1)*dega)*mom
    q = p*exp(dega2*aux1/(1-aux1)) # qn
    if U > p
      aux2 = rand(d)
      while aux2 < aux
        aux2 = rand(d)
      end
      push!(S,aux2)
      U = (U-p)/(1-p)
    elseif U < q
      aux2 = rand(d)
      while aux2 > aux
        aux2 = rand(d)
      end
      push!(S,aux2)
      return S
    else
      aux2 = rand(d)
      while aux2 > aux
        aux2 = rand(d)
      end
      push!(S,aux2)
      U = U/p
    end
  end
end

# Algorithm 3 - 2 (conditional probabilities)
# @Input: S^+(α,ρ), δ, dega = δ*γ, dega2 = -1/(1-exp(-dega)), mom = E(S_1^γ)
function chiStable2(d::PositiveStable,δ::Real,dega::Real,dega2::Real,mom::Real,shift::Integer)
  U = rand(Uniform())
  n = shift
  S = Float64[]
  while true
    n = n+1
    aux = exp(δ*n)
    aux0 = exp(δ*(n+1))
    p = cdf(d,aux)/cdf(d,aux0) # 1-pn
    aux1 = exp(-(n+1)*dega)*mom
    q = p*exp(dega2*aux1/(1 -aux1)) # qn
    if U > p
      aux2 = rand(d)
      while aux2 < aux || aux2 > aux0
        aux2 = rand(d)
      end
      push!(S,aux2)
      U = (U-p)/(1 -p)
    elseif U < q
      aux2 = rand(d)
      while aux2 > aux
        aux2 = rand(d)
      end
      push!(S,aux2)
      return S
    else
      aux2 = rand(d)
      while aux2 > aux
        aux2 = rand(d)
      end
      push!(S,aux2)
      U = U/p
    end
  end
end

import Distributions.rand

# Algorithm 9
# @Input: PerfectSupremumStable
# @Output: A random sample from a the law d, σ, and counter of missed coalescences
function rand(d::PerfectSupremumStable)
  ar = d.α*d.ρ
  et1 = 1 /(ar+d.η)
  ra = 1 /d.α
  rar = 1 /ar
  # Conditionally positive stable distribution
  dPos = PositiveStable(d.α,d.β)
  sep = d.d - d.δ
  e2 = 1 /(1 -exp(-sep))
  mom = mellin(dPos,d.γ)
  dega = d.δ*d.γ
  dega2 = -1 /(1 -exp(-dega))
  mAst = d.mAst + Int(ceil( max(0,log(mom)/dega) + rar ))
  # θ sequence
  U = Float64[]
  Λ = Float64[]
  # Line 1
  x = Inf
  σ = 0
  # Line 2
  lagU = rand(Uniform(),d.Δ)
  lagΛ = 1 .+ rand(Bernoulli(1 - d.ρ),d.Δ) .* (rand(Uniform(),d.Δ) .^ rar .- 1)
  lagS = rand(dPos,d.Δ)

  # Line 3 (first iteration)
  # Line 4 (Algorithm 3 with unconditional probabilities)
  S = rand(dPos,mAst)
  append!(S,chiStable1(dPos,d.δ,dega,dega2,mom,mAst-1))

  R = Float64[]
  C = [0.]
  F = Float64[]
  # Last value of C
  endC = 0.
  # First iteration
  # Line 5
  while (length(R) == σ) || (length(S) >= length(C))
    t = length(C)
    # Lines 6 and 7
    (C1,F1) = downRW(d.d,rar,d.η,et1,d.κ,x-endC)
    append!(F,F1)
    append!(C,endC .+ C1)
    endC = C[end]
    # Lines 8 and 9
    if reflectedProcess(d.d,d.η,et1,d.κ,x-endC)
      # Line 10
      x = endC + d.κ
      R1 = [maximum(C[t:end])]
      for i = (t-1):-1:(length(R)+1)
        push!(R1,max(R1[end],C[i]))
      end
      append!(R,R1[end:-1:1] .- C[(length(R)+1):t])
    else # Line 11
      # Lines 12 and 13
      (C1,F1) = upRW(d.d,d.η,et1,d.κ,x-endC)
      append!(F,F1)
      append!(C,endC .+ C1)
      endC = C[end]
    end # Line 14
  end # Line 15
  # Line 16
  while length(U) < length(S)
    T = rand(Poisson(d.α*(1-d.ρ)*F[length(U)+1]))
    if T == 0
      push!(U,exp(-d.α*F[length(U)+1]))
      push!(Λ,1)
    else
      push!(U,exp(-d.α*F[length(U)+1]*rand(Beta(1,T))))
      push!(Λ,exp(-F[length(U)])/U[end]^ra)
    end
  end
  # Lines 17, 18 and 19
  σ = 1
  D = exp(R[σ])*( e2*exp(-sep*(length(S)-σ-1)) + sum(S[σ:end] .*
    exp.(-d.d .* (0:(length(S)-σ))) .* (1 .- U[σ:end]) .^ ra) )
  T1 = (1 - U[σ])^ra * S[σ]
  T2 = U[σ]^ra * D
  T2 = Λ[σ] * (T1+T2)
  coal = false
  # Run chain 'forward'
  if T2 <= T1
    coal = true
    D = rand(Uniform())^rar * T1
  else
    D = T2
  end
  for i = 1:d.Δ
    T1 = (1 - lagU[i])^ra * lagS[i]
    T2 = lagU[i]^ra * D
    T2 = lagΛ[i] * (T1+T2)
    if T2 <= T1
      coal = true
      D = rand(Uniform())^rar * T1
    else
      D = T2
    end
  end
  # Line 20
  if coal # Coalescence!
    # Line 21
    return (D,σ)
  end # Line 22
  # Line 3 (second to last iterations)
  while true
    # Line 4 (Algorithm 3 with conditional probabilities)
    append!(S,chiStable2(dPos,d.δ,dega,dega2,mom,length(S)-σ-1))
    # Line 5
    while (length(R) == σ) || (length(S) >= length(C))
      t = length(C)-1
      # Lines 6 and 7
      (C1,F1) = downRW(d.d,rar,d.η,et1,d.κ,x-endC)
      append!(F,F1)
      append!(C,endC .+ C1)
      endC = C[end]
      # Lines 8 and 9
      if reflectedProcess(d.d,d.η,et1,d.κ,x-endC)
        # Line 10
        x = endC + d.κ
        R1 = [maximum(C[t:end])]
        for i = (t-1):-1:(length(R)+1)
          push!(R1,max(R1[end],C[i]))
        end
        append!(R,R1[end:-1:1] .- C[(length(R)+1):t])
      else # Line 11
        # Lines 12 and 13
        (C1,F1) = upRW(d.d,d.η,et1,d.κ,x-endC)
        append!(F,F1)
        append!(C,endC .+ C1)
        endC = C[end]
      end # Line 14
    end # Line 15
    # Line 16
    while length(U) < length(S)
      T = rand(Poisson(d.α*(1-d.ρ)*F[length(U)+1]))
      if T == 0
        push!(U,exp(-d.α*F[length(U)+1]))
        push!(Λ,1)
      else
        push!(U,exp(-d.α*F[length(U)+1]*rand(Beta(1,T))))
        push!(Λ,exp(-F[length(U)]/U[end])^ra)
      end
    end
    # Lines 17, 18 and 19
    σ += 1
    D = exp(R[σ])*( e2*exp(-sep*(length(S)-σ-1)) + sum(S[σ:end] .*
      exp.(-d.d .* (0:(length(S)-σ))) .* (1 .- U[σ:end]) .^ ra) )
    T1 = (1 - U[σ])^ra * S[σ]
    T2 = U[σ]^ra * D
    T2 = Λ[σ] * (T1+T2)
    # Line 24
    if T2 <= T1 # Coalescence!
      D = rand(Uniform())^rar * T1
      # Run chain 'forward'
      for i = (σ-1):-1:1
        T1 = (1 - U[i])^ra * S[i]
        T2 = U[i]^ra * D
        T2 = Λ[i] * (T1+T2)
        if T2 <= T1
          D = rand(Uniform())^rar * T1
        else
          D = T2
        end
      end
      for i = 1:d.Δ
        T1 = (1 - lagU[i])^ra * lagS[i]
        T2 = lagU[i]^ra * D
        T2 = lagΛ[i] * (T1+T2)
        if T2 <= T1
          D = rand(Uniform())^rar * T1
        else
          D = T2
        end
      end
      # Line 25
      return (D,σ)
    end # Line 26
  end # Line 27
end

# @Input: SupremumStable
# @Output: A random sample from a the law d
function rand(d::SupremumStable)
  return (β == -1 && α <= 1) ? 0. : rand(sampler(d))[1]
end

export SupremumStable, rand, ExactSampler, minimum, maximum, insupport, mean, params
