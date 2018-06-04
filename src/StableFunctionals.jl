__precompile__(true)

module StableFunctionals

  using Distributions
  using StatsBase
  using LambertW
  using FastGaussQuadrature

  importall Distributions

  export
    # distribution types
    Stable,
    StablePositive,
    StableUnilateral,
    StableSupremum,

    # methods
    rand,
    maximum,
    minimum,
    insupport,
    pdf,
    cdf,
    mgf,
    cf,
    mellin,
    mean,
    params

  # Source Files

  include("stable.jl")
  include("stablepositive.jl")
  include("stableunilateral.jl")
  include("stablesupremum.jl")

end
