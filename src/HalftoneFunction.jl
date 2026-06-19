module HalftoneFunction
using QuadGK,Roots,OffsetArrays
export tryFunc1,tryFunc2,areaBelow,HalftoneApprox
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

function tryFunc2(x::AbstractFloat)
  (1-(x-1)^2)^(3/4)
end

function scale(x::AbstractFloat)
  (1-(x-1)^2)^(3/4)
end

function invScale(x::AbstractFloat)
  1-√(1-x^(4/3))
end

# upperBound and lowerBound are a halftone pair; upperBound(x)*lowerBound(y)
# (with upperBound extended to be periodic the same way as lowerBound) is an
# exact halftone grid with approximately elliptical small dots, because a
# sphere slice has area proportional to its height, regardless of where along
# the sphere it is sliced. The halftone function therefore must be between them.

function upperBound(x::AbstractFloat)
  √(1-(x-1)^2)
end

function lowerBound(x::AbstractFloat)
  sin(x*π/2)
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
  diagCross=find_zero(x->h(x)-√z,[0,1])
  diagCross,find_zero(x->h(x)-z/h(y),[0,1])
end

function areaBelow(h::Function,z::AbstractFloat)
  diagCross=find_zero(x->h(x)-√z,[0,1])
  integrand=y->find_zero(x->h(x)-z/h(y),[0,1])
  quadgk(integrand,diagCross,1)[1]*2+diagCross^2
end

function newPoints(T::DataType,n::Integer)
  ret=OffsetVector(T[],-1)
  for i in 0:n
    push!(ret,invScale(lowerBound(T(i)/n)))
  end
  ret
end

struct HalftoneApprox{T}
  points	::OffsetVector{T}
  HalftoneApprox(T::DataType,n::Integer)=new{T}(newPoints(T,n))
end

end # module HalftoneFunction
