module HalftoneFunction
using QuadGK,Roots
export tryFunc1,areaBelow
export outegrand

# The halftone function h(x) is defined as follows:
# h(x) increases as x goes from 0 to 1.
# The area of {(x,y) in [0,1]×[0,1]: h(x)h(y)<z} is z whenever z is in [0,1].
# h(-x)=-h(x).
# h(1-x)=h(x).
# This module tries to compute an approximation to the halftone function.

function tryFunc1(x::AbstractFloat)
  (1-(x-1)^2)^(2/3)
end

# +-----------------+
# |      *          |
# |      *          |
# |       *         |
# |diagCross*       |
# |         | *     |
# |         |   ****|
# |         |       |
# |         |       |
# +-----------------+

function outegrand(h::Function,z::T,y::T) where T<:AbstractFloat
  diagCross=find_zero(x->h(x)-sqrt(z),[0,1])
  diagCross,find_zero(x->h(x)-z/y,[diagCross,1])
end

function areaBelow(h::Function,z::AbstractFloat)
  diagCross=find_zero(x->h(x)-sqrt(z),[0,1])
  integrand=y->find_zero(x->h(x)-z/y,[diagCross,1])
  quadgk(integrand,diagCross,1)*2+diagCross^2
end

end # module HalftoneFunction
