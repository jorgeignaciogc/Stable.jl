struct ConvexMinorantMeander <: Sampleable{Multivariate,Continuous}
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
  ConvexMinorantMeander(a,b,drift,delta,gamma,kappa,epsilon,lag,sLag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2
    || 0 >= gamma || gamma >= a || 0>drift || drift>1 || 0>delta || delta >= drift || kappa < max(1,log(2)/(3*etaF(drift))) || epsilon<=0 || lag<0 || sLag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], 1>d>δ>0, α>γ>0, κ>max(1,log(2)/(3η)), ϵ>0, lag≥0, sLag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, drift, etaF(drift), delta, gamma, kappa, epsilon, lag, sLag)
  ConvexMinorantMeander(a,b,epsilon, lag, sLag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0 || lag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0, lag≥0, sLag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, lag, sLag)
  ConvexMinorantMeander(a,b,epsilon, lag) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0 || lag<0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0, lag≥0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, lag, 15)
  ConvexMinorantMeander(a,b,epsilon) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 || epsilon<=0) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1], ϵ>0") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3, epsilon, Int(floor(10+(30-log(epsilon)/log(2))*a-10*b)),15)
  ConvexMinorantMeander(a,b) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 ) ?
    error("Parameters' requirements unmet:\n α∈(0,2], β∈[-1,1]") :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2, .95, etaF(.95),.94,a*.95,max(1,log(2)/(3*.95))+3,.5^32,Int(ceil(10+60*a-10*b)),15)
end

import Distributions.rand
# @Input: StableSupremumExact
# @Output: A random sample from a the law d, σ, and counter of missed coalescences
# Note: We follow the algorithms in the EpsStrongStable paper, with slight modifications
function rand(d::ConvexMinorantMeander)
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
      X = Float64[D[end]]
      if n > 1
        ell = Float64[1.0]
        ell = append!(U[n-1:-1:1])
      else
        ell = Float64[1.0]
      end

      # Step 3 in Algorithm 1
      for j = 1:(n-1)
        i = n-j
        push!(X, (1-U[i])^ra*S[i] + U[i]^ra*X[end])
      end
  else
      X = Float64[DD]
      ell = Float64[1.0]
  end

  if d.lag > 0
    append!(ell, lagU[d.lag:-1:1])
    for i = 1:d.lag
      push!(X,(1-lagU[i])^ra * lagS[i] + lagU[i]^ra * X[end])
    end
  end
  return (X,ell)
end
