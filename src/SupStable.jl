__precompile__(true)

module SupStable

  using Distributions
  using StatsBase
  using LambertW
  using FastGaussQuadrature

  importall Distributions

  export
    # distribution types
    PositiveStable,
    SupremumStable,
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

  include("positivestable.jl")
  include("supremumstable.jl")
  include("stable.jl")

end
