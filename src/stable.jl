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
function pdf(d::Stable,X::AbstractArray)
  if x==0
    return gamma(1+1/d.a)*sin(pi*d.rho)/pi
  end
  m=10^3
  delta = 1
    #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  #

  if d.a>1
    a0 = d.a/(abs(d.a-1)*2)
    raa = 1/(d.a-1)
    a1 = x==0 ? 0.0 : abs(x)^raa
    a2 = -a1.^d.a
    s1 = d.rho
    s2 = 1-d.rho
    a01 = a0*s1
    a02 = a0*s2
    if x>0
      if x<delta
        v = 1:Int(floor(18*d.a))
        w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*gamma(v+1))
        return sum(w.*abs(x).^(v-1))
      else
        nodes, weights = gausslegendre(m)
        seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
        return a01*a1*(dot(seq1.*exp(a2.*seq1),weights))
      end
    elseif x<0
      if x>-delta
        v = 1:Int(floor(18*d.a))
        w1 = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir1.*v)./(pi.*gamma(v+1))
        return sum(w1.*abs(x).^(v-1))
      else
        nodes, weights = gausslegendre(m)
        seq2 = auxV2(nodes .*s2+s1,d.a,-d.theta)
        return a02*a1*(dot(seq2.*exp(a2.*seq2),weights))
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = nodes./(1-nodes).^3
    pir = pi*d.rho
    mat = abs(x)*nodes1
    mat = pdf(UnilateralStable(d.a),mat)
    if x>0
      fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
      return dot(fC.*mat,weights)
    else
      fC1 = nodes2.*pdf(Cauchy(cos(pir),sin(pir)),nodes1)
      return dot(fC1.*mat,weights)
    end
  end
end

function pdf(d::Stable,X::AbstractArray)
  m=10^3
  delta = 1
  x=vec(X)
  res = zeros(0)
  #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  v = 1:Int(floor(18*d.a))
  w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*gamma(v+1))
  w1 = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir1.*v)./(pi.*gamma(v+1))
  #
  aux = gamma(1+1/d.a)*sin(pi*d.rho)/pi
  if d.a>1
    a0 = d.a/((d.a-1)*2)
    raa = 1/(d.a-1)
    a1 = abs(x).^raa
    a2 = -a1.^d.a
    nodes, weights = gausslegendre(m)
    s1 = d.rho
    s2 = 1-d.rho
    seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
    seq2 = auxV2(nodes .*s2+s1,d.a,-d.theta)
    a01 = a0*s1
    a02 = a0*s2
    for i=1:length(x)
      if x[i]>0
        if x[i]<delta
          push!(res,sum(w.*abs(x[i]).^(v-1)))
        else
          push!(res,a01*a1[i]*(dot(seq1.*exp(a2[i].*seq1),weights)))
        end
      elseif x[i]<0
        if x[i]>-delta
          push!(res,sum(w1.*abs(x[i]).^(v-1)))
        else
          push!(res,a02*a1[i]*(dot(seq2.*exp(a2[i].*seq2),weights)))
        end
      else
        push!(res,aux)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = nodes./(1-nodes).^3
    pir = pi*d.rho
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    fC1 = nodes2.*pdf(Cauchy(cos(pir),sin(pir)),nodes1)
    mat = abs(x)*transpose(nodes1)
    mat = pdf(UnilateralStable(d.a),mat)
    for i=1:length(x)
      if x[i]>0
        push!(res,dot(fC.*mat[i,:],weights))
      elseif x[i]<0
        push!(res,dot(fC1.*mat[i,:],weights))
      else
        push!(res,aux)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.cdf
function cdf(d::Stable,X::AbstractArray)
  if x==0
    return 1-d.rho
  end
  m=10^3
  delta = 1
    #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  #

  if d.a>1
    raa = 1/(d.a-1)
    a1 = x==0 ? 0.0 : abs(x)^raa
    a2 = -a1.^d.a
    s1 = d.rho
    s2 = 1-d.rho
    if x>0
      if x<delta
        v = 1:Int(floor(18*d.a))
        w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*v.*gamma(v+1))
        return sum(w.*abs(x).^v)+1-d.rho
      else
        nodes, weights = gausslegendre(m)
        seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
        return 1-s1*(dot(exp(a2.*seq1),weights))/2
      end
    elseif x<0
      if x>-delta
        v = 1:Int(floor(16*d.a))
        w1 = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir1.*v)./(pi.*v.*gamma(v+1))
        return 1-d.rho-sum(w1.*abs(x).^v)
      else
        nodes, weights = gausslegendre(m)
        seq2 = auxV2(nodes .*s2+s1,d.a,-d.theta)
        return s2*(dot(exp(a2.*seq2),weights))/2
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = 1./(1-nodes).^2
    pir = pi*d.rho
    mat = abs(x)*nodes1
    mat = cdf(UnilateralStable(d.a),mat)
    if x>0
      fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
      return 1-d.rho+dot(fC.*mat,weights)
    else
      fC1 = nodes2.*pdf(Cauchy(cos(pir),sin(pir)),nodes1)
      return dot(fC1.*(1-mat),weights)
    end
  end
end

function cdf(d::Stable,X::AbstractArray)
  m=10^3
  delta = 1
  x=vec(X)
  res = zeros(0)
  #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  v = 1:Int(floor(18*d.a))
  w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*v.*gamma(v+1))
  w1 = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir1.*v)./(pi.*v.*gamma(v+1))
  #
  aux = 1-d.rho
  if d.a>1
    raa = d.a/(d.a-1)
    a2 = -abs(x).^raa
    nodes, weights = gausslegendre(m)
    s1 = d.rho
    s2 = 1-d.rho
    seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
    seq2 = auxV2(nodes .*s2+s1,d.a,-d.theta)
    s1 = s1/2
    s2 = s2/2
    for i=1:length(x)
      if x[i]>0
        if x[i]<delta
          push!(res,aux+sum(w.*abs(x[i]).^v))
        else
          push!(res,1-s1*(dot(exp(a2[i].*seq1),weights)))
        end
      elseif x[i]<0
        if x[i]>-delta
          push!(res,aux-sum(w1.*abs(x[i]).^v))
        else
          push!(res,s2*(dot(exp(a2[i].*seq2),weights)))
        end
      else
        push!(res,aux)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = 1./(1-nodes).^2
    pir = pi*d.rho
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    fC1 = nodes2.*pdf(Cauchy(cos(pir),sin(pir)),nodes1)
    mat = abs(x)*transpose(nodes1)
    mat = cdf(UnilateralStable(d.a),mat)
    for i=1:length(x)
      if x[i]>0
        push!(res,aux+dot(fC.*mat[i,:],weights))
      elseif x[i]<0
        push!(res,dot(fC1.*(1-mat[i,:]),weights))
      else
        push!(res,aux)
      end
    end
  end
  return reshape(res,size(X))
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
