__precompile__(true)

module SupStable

  using Distributions
  using StatsBase
  using LambertW
  using FastGaussQuadrature

  importall Distributions

  export
    # distribution types
    StablePositive,
    StableSupremum,
    StableUnilateral,
    Stable,

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
  include("stableunilateral.jl")
  include("stable.jl")

end
