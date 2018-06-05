__precompile__(true)

module StableSupremum

  using Distributions
  using StatsBase
  using LambertW
  using FastGaussQuadrature

  importall Distributions

  export
    # distribution types
    StablePositive,
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

  include("stablepositive.jl")
  include("stablesupremum.jl")

end
