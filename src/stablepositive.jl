using Distributions
using FastGaussQuadrature
importall Distributions

struct StablePositive <: ContinuousUnivariateDistribution
  a::Float64
  b::Float64
  theta::Float64
  rho::Float64
  StablePositive(a,b) = ( b<-1 || b>1 || (b==-1 && a<=1) || a<=0 || a>2 ) ?
    error("Parameters' requirements unmet:\n α∈(0,2] and β∈[-1,1])") :
    (a==2 ? new(2,0.0,0.0,.5) :
    new(a,b,b*(a<=1 ? 1 : (a-2)/a),(1+b*(a<=1 ? 1 : (a-2)/a))/2) )
end

import Distributions.params
function params(d::StablePositive)
  return (d.a,d.b,d.theta,d.rho)
end

import Distributions.minimum
function minimum(d::StablePositive)
  return 0
end

import Distributions.maximum
function maximum(d::StablePositive)
  return Inf
end

import Distributions.insupport
function insupport(d::StablePositive,x::Real)
  return x>=0
end

########################
# BASIC SIMULATION TOOLS
########################

import Distributions.rand
function rand(d::StablePositive)
  if d.a<=1
    r = (d.b+1)/2
  else
    r = (d.b*(d.a-2)/d.a+1)/2
  end
  ra = r*d.a
  u1 = rand(Uniform())*pi
  uni1 = ((((sin(r*u1)^r*sin((1-r)*u1)^(1-r))/sin(u1))^(1/(1-r)))/rand(Exponential()))^(1/r-1)
  u2 = rand(Uniform())*pi
  uni2 = ((((sin(ra*u2)^ra*sin((1-ra)*u2)^(1-ra))/sin(u2))^(1/(1-ra)))/rand(Exponential()))^(1/ra-1)
  return (uni2/uni1)^r
end

function rand(d::StablePositive,n::Int)
  if d.a<=1
    r = (d.b+1)/2
  else
    r = (d.b*(d.a-2)/d.a+1)/2
  end
  ra = r*d.a
  inva = 1/d.a
  invra = inva-r
  cr = 1-r
  cra = 1-ra
  invcr = 1/cr
  invcra = 1/cra

  if d.b==1 && d.a==1
    return ones(n)
  elseif d.b==1 && d.a<1
    u2 = rand(Uniform(0,pi),n)
    e2 = rand(Exponential(),n)
    ca = 1- d.a
    invra = inva-1
    return [((((sin(d.a*u2[i])*sin(ca*u2[i])^invra)/sin(u2[i]))^inva)/e2[i]^invra) for i=1:n]
  elseif d.b==-1 && d.a>1
    u1 = rand(Uniform(0,pi),n)
    e1 = rand(Exponential(),n)
    return [e1[i]^cr /
    (((sin(r*u1[i])^r*sin(cr*u1[i])^cr)/sin(u1[i]))) for i=1:n]
  else
    u1 = rand(Uniform(0,pi),n)
    u2 = rand(Uniform(0,pi),n)
    e1 = rand(Exponential(),n)
    e2 = rand(Exponential(),n)
    return [((((sin(ra*u2[i])^r*sin(cra*u2[i])^invra)/sin(u2[i]))^inva)/e2[i]^invra) /
    ((((sin(r*u1[i])^r*sin(cr*u1[i])^cr)/sin(u1[i])))/e1[i]^cr) for i=1:n]
  end
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
function pdf(d::StablePositive,x::Real)
  if x<0
    return 0.0
  elseif x==0
    return gamma(1+1/d.a)*sin(pi*d.rho)/(d.rho*pi)
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
    if x<delta
      v = 1:Int(floor(18*d.a))
      w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*gamma(v+1))
      return sum(w.*abs(x).^(v-1))
    else
      nodes, weights = gausslegendre(m)
      seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
      return a01*a1*(dot(seq1.*exp(a2.*seq1),weights))
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = nodes./(1-nodes).^3
    pir = pi*d.rho
    mat = abs(x)*nodes1
    mat = pdf(StableUnilateral(d.a),mat)
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    return dot(fC.*mat,weights)/d.rho
  end
end

function pdf(d::StablePositive,X::AbstractArray)
  m=10^3
  delta = 1
  x=vec(X)
  res = zeros(0)
  #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  v = 1:Int(floor(18*d.a))
  w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./((d.rho*pi).*gamma(v+1))
  w1 = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir1.*v)./((d.rho*pi).*gamma(v+1))
  #
  aux = gamma(1+1/d.a)*sin(pi*d.rho)/(d.rho*pi)
  if d.a>1
    a0 = d.a/((d.a-1)*2)
    raa = 1/(d.a-1)
    a1 = abs(x).^raa
    a2 = -a1.^d.a
    nodes, weights = gausslegendre(m)
    s1 = d.rho
    s2 = 1-d.rho
    seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
    weights = (a0*s1)*weights
    for i=1:length(x)
      if x[i]>0
        if x[i]<delta
          push!(res,sum(w.*abs(x[i]).^(v-1)))
        else
          push!(res,a1[i]*(dot(seq1.*exp(a2[i].*seq1),weights)))
        end
      elseif x[i]<0
        push!(res,0.0)
      else
        push!(res,aux)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/(d.rho*2)
    nodes1 = nodes./(1-nodes)
    nodes2 = nodes./(1-nodes).^3
    pir = pi*d.rho
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    mat = abs(x)*transpose(nodes1)
    mat = pdf(StableUnilateral(d.a),mat)
    for i=1:length(x)
      if x[i]>0
        push!(res,dot(fC.*mat[i,:],weights))
      elseif x[i]<0
        push!(res,0.0)
      else
        push!(res,aux)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.cdf
function cdf(d::StablePositive,x::Real)
  if x<=0
    return 0.0
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
    if x<delta
      v = 1:Int(floor(18*d.a))
      w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./(pi.*v.*gamma(v+1))
      return sum(w.*abs(x).^v)/d.rho
    else
      nodes, weights = gausslegendre(m)
      seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
      return 1-s1*(dot(exp(a2.*seq1),weights))/(2*d.rho)
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/2
    nodes1 = nodes./(1-nodes)
    nodes2 = 1./(1-nodes).^2
    pir = pi*d.rho
    mat = abs(x)*nodes1
    mat = cdf(StableUnilateral(d.a),mat)
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    return dot(fC.*mat,weights)/d.rho
  end
end

function cdf(d::StablePositive,X::AbstractArray)
  m=10^3
  delta = 1
  x=vec(X)
  res = zeros(0)
  #
  pir = pi*d.rho
  pir1 = pi*(1-d.rho)
  v = 1:Int(floor(18*d.a))
  w = (-1).^(v-1).*gamma(v./d.a+1).*sin(pir.*v)./((pi*d.rho).*v.*gamma(v+1))
  #
  if d.a>1
    raa = d.a/(d.a-1)
    a2 = -abs(x).^raa
    nodes, weights = gausslegendre(m)
    s1 = d.rho
    s2 = 1-d.rho
    seq1 = auxV2(nodes .*s1+s2,d.a,d.theta)
    s1 = s1/(2*d.rho)
    for i=1:length(x)
      if x[i]>0
        if x[i]<delta
          push!(res,sum(w.*abs(x[i]).^v))
        else
          push!(res,1-s1*(dot(exp(a2[i].*seq1),weights)))
        end
      else
        push!(res,0.0)
      end
    end
  else
    nodes, weights = gausslegendre(m)
    nodes = nodes/2+.5
    weights = weights/(2*d.rho)
    nodes1 = nodes./(1-nodes)
    nodes2 = 1./(1-nodes).^2
    pir = pi*d.rho
    fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes1)
    mat = abs(x)*transpose(nodes1)
    mat = cdf(StableUnilateral(d.a),mat)
    for i=1:length(x)
      if x[i]>0
        push!(res,aux+dot(fC.*mat[i,:],weights))
      else
        push!(res,0.0)
      end
    end
  end
  return reshape(res,size(X))
end


import Distributions.mgf
function mgf(d::StablePositive,x::Real)
  if x==0
    return 1
  end
  if d.a==2
    return 2*exp(x^2/4)*cdf(Normal(0,sqrt(2)),x/2)
  end
  if d.b==-1 && x>=-1
    v = 0:Int(floor(18*d.a))
    w = 1./gamma(v+1)
    return dot(w,x.^v)
  end

  nodes, weights = gausslegendre(m)
  nodes = nodes/2+.5
  weights = weights/(2*d.rho)
  nodes2 = 1./(1-nodes).^2
  nodes = nodes./(1-nodes)
  pir = pi*d.rho
  fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes)
  if d.b==-1
    mat = exp(-(abs(x).*nodes).^d.a)
    push!(res,dot(fC.*mat,weights))
  else
    for i = 1:length(x)
      if x>0
        return Inf
      else
        mat = exp(-(abs(x).*nodes).^d.a)
        return dot(fC.*mat[i,:],weights)
      end
    end
  end
end

function mgf(d::StablePositive,X::AbstractArray)
  x=vec(X)
  if d.a==2
    return 2.*exp(x.^2./4)*cdf(Normal(0,sqrt(2)),x/2)
  end
  nodes, weights = gausslegendre(m)
  nodes = nodes/2+.5
  weights = weights/(2*d.rho)
  nodes2 = 1./(1-nodes).^2
  nodes = nodes./(1-nodes)
  pir = pi*d.rho
  fC = nodes2.*pdf(Cauchy(-cos(pir),sin(pir)),nodes)
  res = zeros(0)
    mat = exp(-(abs(x)*transpose(nodes)).^d.a)
  if d.b==-1
    v = 0:Int(floor(18*d.a))
    w = 1./gamma(v+1)
    if x[i]==0
      push!(res,1)
    elseif x[i]>=-1
      push!(res,dot(w,x[i].^v))
    else
      push!(res,dot(fC.*mat[i,:],weights))
    end
  else
    for i = 1:length(x)
      if x[i]>0
        push!(res,Inf)
      elseif x[i]<0
        push!(res,dot(fC.*mat[i,:],weights))
      else
        push!(res,1)
      end
    end
  end
  return reshape(res,size(X))
end

import Distributions.mean
function mean(d::StablePositive)
  if d.a<=1
    return Inf
  end
  return sin(pi*d.rho)/(d.a*d.rho*sin(pi/d.a)*gamma(1+1/d.a))
end

import Distributions.var
function var(d::StablePositive)
  if d.a<2 && d.b!=-1
    return Inf
  end
  return 2/gamma(1+2/d.a)-1/gamma(1+1/d.a)^2
end

function mellin(d::StablePositive,x::Complex)
  if (real(x)>=d.a && (d.a<=1 || (d.a<2 && d.b!=-1))) || real(x)<=-1
    return Inf
  end
  if (d.a>1 && d.b==-1) || d.a==2
    return gamma(1+x)/gamma(1+x/d.a)
  end
  return (sin(pi*d.rho*x)*gamma(1+x))/(d.a*d.rho*sin(pi*x/d.a)*gamma(1+x/d.a))
end

function mellin(d::StablePositive,x::Real)
  if (real(x)>=d.a && (d.a<=1 || (d.a<2 && d.b!=-1))) || real(x)<=-1
    return Inf
  end
  if (d.a>1 && d.b==-1) || d.a==2
    return gamma(1+x)/gamma(1+x/d.a)
  end
  return (sin(pi*d.rho*x)*gamma(1+x))/(d.a*d.rho*sin(pi*x/d.a)*gamma(1+x/d.a))
end

export Stable, rand, minimum, maximum, insupport, pdf, cdf, cf, mgf, mean, mellin, params
