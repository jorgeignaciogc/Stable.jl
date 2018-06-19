using Distributions
using FastGaussQuadrature
importall Distributions

struct Stable <: ContinuousUnivariateDistribution
  a::Float64
  b::Float64
  theta::Float64
  rho::Float64
  Stable(a,b) = ( b<-1 || b>1 || a<=0 || a>2 ) ? error("Parameters' requirements unmet:\n α∈(0,2] and β∈[-1,1])") :
    (a==2 ? Normal(0,sqrt(2)) :
    (a==1 ? (abs(b)<1 ? Cauchy(-cos(pi*(1+b)/2),sin(pi*(1+b)/2)) :
    new(a,b,b,(b+1)/2)) :
    a<1 && abs(b)==1 ? UnilateralStable(a,b) : new(a,b,b*(a<=1 ? 1 :
    (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2)))
end

import Distributions.params
function params(d::Stable)
  return (d.a,d.b,d.theta,d.rho)
end

import Distributions.minimum
function minimum(d::Stable)
  return -Inf
end

import Distributions.maximum
function maximum(d::Stable)
  return Inf
end

import Distributions.insupport
function insupport(d::Stable,x::Real)
  return isfinite(x)
end

########################
# BASIC SIMULATION TOOLS
########################

import Distributions.rand
function rand(d::Stable) # Chambers. Mellows, and Stuck
  B = rand(Bernoulli(d.rho))
  return B==1 ? rand(PositiveStable(d.a,d.b)) : -rand(PositiveStable(d.a,-d.b))
end

function rand(d::Stable,n::Int) # Chambers. Mellows, and Stuck
  B = rand(Bernoulli(d.rho),n)
  N = Integer(sum(B))
  S1 = rand(PositiveStable(d.a,d.b),N)
  S2 = rand(PositiveStable(d.a,-d.b),n-N)
  res = zeros(0)
  j = 0
  for i = 1:n
    if B[i] == 1
      j = j+1
      push!(res,S1[j])
    else
      push!(res,-S2[i-j])
    end
  end
  return res
end
########################
# PDF, CDF
########################

function auxV2(x::AbstractArray,a::Real,th::Real)
  y = (pi/2).*x
  t = (pi/2)*a*th
  return ((sin(a.*y+t).^a)./cos(y)).^(1/(1-a)).*cos((a-1).*y+t)
end

import Distributions.pdf
function pdf(d::Stable,x::Real)
  return x==0 ? gamma(1+1/d.a)*sin(pi*d.rho)/pi : (x>0 ? pdf(StablePositive(d.a,d.b),x)*d.rho : pdf(StablePositive(d.a,-d.b),-x)*(1-d.rho))
end

function pdf(d::Stable,X::AbstractArray)
  return pdf(StablePositive(d.a,d.b),X).*d.rho + pdf(StablePositive(d.a,-d.b),-X).*(1-d.rho)
end

import Distributions.cdf
function cdf(d::Stable,x::Real)
  return d.rho + cdf(StablePositive(d.a,d.b),x)*d.rho - cdf(StablePositive(d.a,-d.b),-x)*(1-d.rho)
end

function cdf(d::Stable,X::AbstractArray)
  return d.rho + cdf(StablePositive(d.a,d.b),X).*d.rho - cdf(StablePositive(d.a,-d.b),-X).*(1-d.rho)
end

import Distributions.cf
function cf(d::Stable,x::Real)
  if x==0
    return 1
  end
  if d.a==1 && abs(d.b)==1
    return exp(d.b*im*x)
  end
  return exp(-abs(x)^d.a*exp(-pi*d.theta*d.a*sign(x)*im/2))
end

function cf(d::Stable,X::AbstractArray)
  if d.a==1 && abs(d.b)==1
    return exp((d.b*im).*X)
  end
  x = vec(X)
  s1 = exp(-pi*d.theta*d.a*im/2)
  s2 = exp(pi*d.theta*d.a*im/2)
  return reshape(exp(-abs(x).^d.a.*[x[i]>0 ? s1 : s2 for i=1:length(x)]),size(X))
end

import Distributions.mgf
function mgf(d::Stable,x::Real)
  if x==0
    return 1
  end
  if d.a==2
    return exp(x^2)
  elseif d.a==1 && abs(d.b)==1
      return exp(d.b*x)
  elseif x>0
    if d.b==-1 && d.a!=1
        return exp(sign(d.a-1)*x^d.a)
    else
      return Inf
    end
  elseif d.b==1 && d.a!=1
    return exp(sign(d.a-1)*(-x)^d.a)
  end
  return Inf
end

function mgf(d::Stable,X::AbstractArray)
  if d.a==2
    return exp(X.^2)
  elseif d.a==1 && abs(d.b)==1
    return exp(d.b.*X)
  end
  x = vec(X)
  res = exp(sign(d.a-1).*abs(x).^d.a)
  if d.b==-1
    return reshape([x[i]>=0 ? res[i] : Inf for i=1:length(x)],size(X))
  elseif d.b==1
    return reshape([x[i]<=0 ? res[i] : Inf for i=1:length(x)],size(X))
  end
  return reshape([x[i]==0 ? 1 : Inf for i=1:length(x)],size(X))
end

import Distributions.mean
function mean(d::Stable)
  if d.a<=1
    return Inf
  end
  cr=1-d.rho
  return (sin(pi*d.rho)-sin(pi*cr))/(d.a*sin(pi/d.a)*gamma(1+1/d.a))
end

import Distributions.var
function var(d::Stable)
  return Inf
end

function mellin(d::Stable,x::Complex)
  if (real(x)>=d.a && (d.a<=1 || (d.a<2 && d.b!=-1))) || real(x)<=-1
    return Inf
  end
  if (d.a>1 && d.b==-1) || d.a==2
    return gamma(1+x)/(d.a*gamma(1+x/d.a))
  elseif d.a<=1
    r=(1+d.b)/2
  else
    r=(1+d.b*(1-2/d.a))/2
  end
  return (sin(pi*r*x)*gamma(1+x))/(d.a*sin(pi*x/d.a)*gamma(1+x/d.a))
end

export Stable, rand, minimum, maximum, insupport, pdf, cdf, cf, mgf, mean, mellin, params
