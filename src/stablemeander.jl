using Distributions
importall Distributions

struct StableMeander <: ContinuousUnivariateDistribution
  a::Float64
  b::Float64
  theta::Float64
  rho::Float64
  StableMeander(a,b) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 ) ?
  error("Parameters' requirements unmet:\n α∈(0,2] and β∈[-1,1])") :
  #(a==2 ? new(2,0.0,0.0,.5) : Scaled Brownian Meander, specific package currently not supported
  (b==1 && a<=1 ? StableUnilateral(a) :
  #(b==-1 && a>1 ? StablePositive(a,b):
  new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2))#))
end

import Distributions.params
function params(d::StableMeander)
  return (d.a,d.b,d.theta,d.rho)
end

import Distributions.minimum
function minimum(d::StableMeander)
  return 0
end

import Distributions.maximum
function maximum(d::StableMeander)
  return Inf
end

import Distributions.insupport
function insupport(d::StableMeander,x::Real)
  return x>=0
end

import Distributions.mean
function mean(d::StableMeander)
  return pi/(sin(pi/d.a) * gamma(d.rho + 1/d.a) * gamma(1 - d.rho))
end

struct StableMeanderPrecise <: Sampleable{Univariate,Continuous}
  a::Float64
  b::Float64
  theta::Float64
  rho::Float64
  drift::Float64
  eta::Float64
  delta::Float64
  gamma::Float64
  kappa::Float64
  epsilon::Float64
  lag::Int64
  sLag::Int64
  StableMeanderPrecise(a,b,drift,delta,gamma,kappa,epsilon,lag,sLag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2
    || 0 >= gamma || gamma >= a || 0>drift || drift>1 || 0>delta || delta >= drift || kappa < max(1,log(2)/(3*etaF(drift))) || epsilon<=0 || lag<0 || sLag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], 1>d>δ>0, α>γ>0, κ>max(1,log(2)/(3η)), ϵ>0, lag≥0, sLag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, drift, etaF(drift), delta, gamma, kappa, epsilon, lag, sLag)
  StableMeanderPrecise(a,b,epsilon, lag, sLag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0 || lag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0, lag≥0, sLag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, lag, sLag)
  StableMeanderPrecise(a,b,epsilon, lag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0 || lag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0, lag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, lag, 15)
  StableMeanderPrecise(a,b,epsilon) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, Int(floor(10+(30-log(epsilon)/log(2))*a-10*b)),15)
  StableMeanderPrecise(a,b) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 ) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1]") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3,.5^32,Int(ceil(10+60*a-10*b)),15)
end

function PreciseSampler(d::StableMeander,drift::Real,delta::Real,gamma::Real,kappa::Real,epsilon::Real,lag::Integer,sLag::Integer)
  return StableMeanderPrecise(d.a,d.b,drift,delta,gamma,kappa,epsilon,lag,sLag)
end

function PreciseSampler(d::StableMeander,epsilon::Real,lag::Integer,sLag::Integer)
  return StableMeanderPrecise(d.a,d.b,epsilon,lag,sLag)
end

function PreciseSampler(d::StableMeander,epsilon::Real,lag::Integer)
  return StableMeanderPrecise(d.a,d.b,epsilon,lag)
end

function PreciseSampler(d::StableMeander,epsilon::Real)
  return StableMeanderPrecise(d.a,d.b,epsilon)
end

function PreciseSampler(d::StableMeander)
  return StableMeanderPrecise(d.a,d.b)
end

using LambertW
# @Input: the normalized drift d
# @Returns: the value of η>0, 0 < 1/(1+η) < 1, and of κ
function etaF(drift::Real)
  return -lambertw(-drift*exp(-drift),-1)/drift-1
end

# @Input: the output of etaF() and a bound M we want to know if we exceed
# @Returns: a sample of 1{R_0 <= M}
function reflectedProcess(drift::Real,et::Real,et1::Real,M::Real)
  T = rand(Exponential())/et
  if T <= M
    return true
  end
  Sn = [0.0]
  while Sn[end] < T
    push!(Sn,Sn[end]-rand(Exponential())*et1+drift)
  end
  return maximum(Sn[1:(end-1)]) <= M
end

# @Input: the output of etaF() a bound M that it should not exceed,
# and another bound M1 we want to know if we exceed
# @Returns: a sample of 1{R_0<=M1} given R_0<M
function  reflectedProcess(drift::Real,et::Real,et1::Real,M::Real,M1::Real)
  m = exp(-M-drift)
  Sn = [0.0]
  while true
    T =-log(rand(Uniform(m,1)))/et # Exponential conditioned on <= M+d, since T>M+d always yields a rejection
    if T <= M1
      return true
    end
    while Sn[end] < T
      push!(Sn,Sn[end]-rand(Exponential())*et1+drift)
    end
    max = maximum(Sn[1:(end-1)])
    if max <= M
      return max <= M1
    else
      Sn = [0.0]
    end
  end
end

# @Input: normalized drift d, the output of etaF(), and an upper bound x
# @Returns: A simulated run of the random walk up to the hitting time of the barrier
# -2κ conditioned on never exceeding x
function downRW(drift::Real,et::Real,et1::Real,kappa::Real,x::Real)
  Sn = [0.0]
  INCn = [0.0]
  if x == Inf
    while Sn[end] > -2*kappa
      push!(INCn,rand(Exponential()))
      push!(Sn,Sn[end]-INCn[end]+drift)
    end
    return (Sn[2:end],INCn[2:end])
  else
    while true
      while Sn[end] > -2*kappa
        push!(INCn,rand(Exponential()))
        push!(Sn,Sn[end]-INCn[end]+drift)
      end
      if maximum(Sn) < x && reflectedProcess(drift,et,et1,x-Sn[end])
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.0]
        INCn = [0.0]
      end
    end
  end
end

# @Input: normalized drift d, the output of etaF(), an upper bound x0 that
# will be attained and an upper bound x1 that will not
# @Returns: A simulated run of the random walk up to the hitting time of the barrier
# κ conditioned on never exceeding x
function upRW(drift::Real,et::Real,et1::Real,kappa::Real,x::Real)
  Sn = [0.0]
  INCn = [0.0]
  if x == Inf
    while true
      while Sn[end] < kappa
        push!(INCn,rand(Exponential())*et1)
        push!(Sn,Sn[end]-INCn[end]+drift)
      end
      if rand(Uniform()) <= exp(-et*Sn[end])
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.0]
        INCn = [0.0]
      end
    end
  else
    while true
      while Sn[end] < kappa
        push!(INCn,rand(Exponential())*et1)
        push!(Sn,Sn[end]-INCn[end]+drift)
      end
      if rand(Uniform()) <= exp(-et*Sn[end]) && reflectedProcess(drift,et,et1,x-Sn[end]) # The condition maximum(Sn) < x1 is trivially satisfied by construction of κ
        return (Sn[2:end],INCn[2:end])
      else
        Sn = [0.0]
        INCn = [0.0]
      end
    end
  end
end

# Blanchet-Sigman algorithm
# @Input: normalized drift d, the output of etaValue, choice of barrier length M,
# an upper bound x that will not be attained, last value of random walk s,
# @Returns: A simulated run of the random walk up to the first time we know an upper bound for the walk,
# the increments of the RW, the new upper bound, and the last time at which we may compute D_n
function BSAlgorithm(drift::Real,et::Real,et1::Real,kappa::Real,x::Real,s::Real)
  Sn = [Float64(s)]
  INCn = [Float64(0.0)]
  x1 = x
  t = 0

  while true
    # Step 1
    (Sn1,INCn1) = downRW(drift,et,et1,kappa,x1-Sn[end])
    # Step 2
    append!(Sn,Sn[end]+Sn1)
    append!(INCn,INCn1)
    t = t+length(Sn1)
    # Step 3 and 4
    if reflectedProcess(drift,et,et1,x1-Sn[end],kappa)
      x1 = Sn[end]+kappa
      t = t-length(Sn1)
      break
    else
      (Sn1,INCn1) = upRW(drift,et,et1,kappa,x1-Sn[end])
      append!(Sn,Sn[end]+Sn1)
      append!(INCn,INCn1)
      t = t+length(Sn1)
    end
  end

  return (Sn[2:end],INCn[2:end],x1,t)
end

# Algorithm 3 - 1 (unconditional probabilities)
# @Input: dega = delta*gamma, dega2 = -1/(1-exp(-dega)), mom = E(S_1^γ), shift = m-k >= m*-1
function chiStable1(d::StablePositive,delta::Real,dega::Real,dega2::Real,mom::Real,shift::Integer)
  U = rand(Uniform())
  n = shift
  S = Float64[]
  while true
    n = n+1
    aux = exp(delta*n)
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
# @Input: dega = delta*gamma, dega2 = -1/(1-exp(-dega)), mom = E(S_1^γ)
function chiStable2(d::StablePositive,delta::Real,dega::Real,dega2::Real,mom::Real,shift::Integer)
  U = rand(Uniform())
  n = shift
  S = Float64[]
  while true
    n = n+1
    aux = exp(delta*n)
    aux0 = exp(delta*(n+1))
    p = cdf(d,aux)/cdf(d,aux0) # 1-pn
    aux1 = exp(-(n+1)*dega)*mom
    q = p*exp(dega2*aux1/(1-aux1)) # qn
    if U > p
      aux2 = rand(d)
      while aux2 < aux || aux2 > aux0
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

import Distributions.rand
# @Input: StableSupremumExact
# @Output: A random sample from a the law d, σ, and counter of missed coalescences
# Note: We follow the algorithms in the EpsStrongStable paper, with slight modifications
function rand(d::StableMeanderPrecise)
  et1 = 1/(d.eta+1)
  ra = 1/d.a
  # Theta sequence
  U = Float64[]
  # Dominating process
  D = Float64[]
  # Simple overestimate of D[1]
  DD = 0
  dPos = StablePositive(d.a,d.b)

  # Warmup
  lagU = lagS = 0
  # Epsilon after warmup
  if d.lag > 0
    lagU = rand(Beta(d.rho,1),d.lag)
    lagS = rand(dPos,d.lag)
    eps = d.epsilon / prod(lagU)^ra
  else
    eps = d.epsilon
  end

  ar = d.a*d.rho

  reDrift = d.drift/ar # Non-normalized drift
  reDelta = d.delta/ar # Non-normalized δ
  sep = reDrift - reDelta
  e2 = 1/(1-exp(-sep))
  mom = mellin(dPos,d.gamma)
  dega = reDelta*d.gamma
  dega2 = -1/(1-exp(-dega))
  mAst = d.sLag + Integer(max(0,floor(log(mom)/dega)))+1

  # Step 1 in Algorithm 2 (Algorithm 3)
  # Computing (chi_{-1},S_n)
  S = rand(dPos,mAst)
  append!(S,chiStable1(dPos,reDelta,dega,dega2,mom,mAst-2))
  chi = length(S)

  # Simple criterion to stop
  Ex = rand(Exponential())/d.eta
  DD = exp(Ex)*(e2*exp(sep*(1-chi)) + exp(2*reDrift)*sum(S[2:end].*exp.(-reDrift.*(2:length(S)))))

  n = 0

  if DD > eps
      bool = true
      if Ex > d.kappa
        Sn = [0.0]
        while Sn[end] < Ex
          push!(Sn,Sn[end]-rand(Exponential())*et1+d.drift)
        end
        bool = maximum(Sn[1:(end-1)]) <= d.kappa
      end

      t = 0
      x = Inf
      C = Float64[]
      F = Float64[]
      R = Float64[]
      if bool
          (C,F) = downRW(d.drift,d.eta,et1,d.kappa,x)
          C = C/ar
          F = F/d.rho
          while true
            (C1,F1,x,t) = BSAlgorithm(d.drift,d.eta,et1,d.kappa,x,C[end]*ar)
            t = t+length(C)
            append!(C,C1/ar)
            append!(F,F1/d.rho)
            # Compute R
            t1 = length(R)
            append!(R,zeros(t-t1))
            R[end] = maximum(C[t:(end-1)])
            for i = (t-1):-1:(t1+1)
              R[i] = max(C[i],R[i+1])
            end
            R[(t1+1):t] = R[(t1+1):t] - C[(t1+1):t]

            if length(C) >= chi
              break
            end
          end
      else
          (C,F,x,t) = BSAlgorithm(d.drift,d.eta,et1,d.kappa,x,0)
          C = C/ar
          F = F/d.rho
          R = Float64[]
          while true
            (C1,F1,x,t) = BSAlgorithm(d.drift,d.eta,et1,d.kappa,x,C[end]*ar)
            t = t+length(C)
            append!(C,C1/ar)
            append!(F,F1/d.rho)
            # Compute R
            t1 = length(R)
            append!(R,zeros(t-t1))
            R[end] = maximum(C[t:(end-1)])
            for i = (t-1):-1:(t1+1)
              R[i] = max(C[i],R[i+1])
            end
            R[(t1+1):t] = R[(t1+1):t] - C[(t1+1):t]

            if length(C) >= chi
              break
            end
          end
      end

      # Step 2 Part 2 in Algorithm 2 (computing U)
      append!(U,exp(-F))

      # Auxiliary Process
      eS = exp(-reDrift*(1:chi)) .* S .* ((1-U[1:chi]) .^ ra)

      # Step 3 in Algorithm 3
      push!(D, exp(R[1])*(e2*exp(sep*(1-chi)) + exp(2*reDrift)*sum(eS[2:end])) )
      # Step 5 in Algorithm 2
      n = 1
      while eps < exp(C[n] - n*reDrift)*D[end]
        n = n+1 # Step 6 in Algorithm 2
        # Step 7 in Algorithm 2
        append!(S,chiStable2(dPos,reDelta,dega,dega2,mom,chi-n-1))
        chi = length(S)
        # Step 8 in Algorithm 2
        while true
          (C1,F1,x,t) = BSAlgorithm(d.drift,d.eta,et1,d.kappa,x,C[end]*ar)
          t = t+length(C)
          append!(C,C1/ar)
          append!(F,F1/d.rho)
          # Compute R
          t1 = length(R)
          append!(R,zeros(t-t1))
          R[end] = maximum(C[t:(end-1)])
          for i = (t-1):-1:(t1+1)
            R[i] = max(C[i],R[i+1])
          end
          R[(t1+1):t] = R[(t1+1):t] - C[(t1+1):t]

          if length(C) >= chi
            break
          end
        end
        # Step 8 Part 2 in Algorithm 2 (computing U)
        append!(U, exp(-F[(length(U)+1):end]))

        # Auxiliary process
        append!(eS, exp(-reDrift*((length(eS)+1):length(S))) .*
          S[(length(eS)+1):end] .* ((1-U[(length(eS)+1):length(S)]) .^ ra) )
        # Step 9 in Algorithm 2
        push!(D, exp(R[n])*(e2*exp(sep*(n-chi)) + exp((n+1)*reDrift)*sum(eS[(n+1):end])) )
      end # Precision reached

      # Step 11 in Algorithm 2 would give σ=n

      # Step 2 in Algorithm 1
      X = D[end]

      # Step 3 in Algorithm 1
      for j = 1:(n-1)
        i = n-j
        X = (1-U[i])^ra*S[i] + U[i]^ra*X
      end
  else
      X = DD
  end

  if d.lag > 0
    for i = 1:d.lag
      X = (1-lagU[i])^ra * lagS[i] + lagU[i]^ra * X
    end
  end
  return (X,n)
end

# @Input: StableSupremum
# @Output: A random sample from a the law d
function rand(d::StableMeander)
  return rand(PreciseSampler(d))
end

export StableMeander, rand, PreciseSampler, minimum, maximum, insupport, mean, params
