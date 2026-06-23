module HalftoneFunction
using QuadGK,Roots,OffsetArrays,Printf
export tryFunc1,tryFunc2,areaBelow,HalftoneApprox,ht
export outegrand,htError,adjust!

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

function ht(x::T,hta::HalftoneApprox{T}) where T<:AbstractFloat
  if x<0 || x>1
    throw(DomainError(x,"The argument to ht must be between 0 and 1 inclusive."))
  end
  pos=floor(Int,x*lastindex(hta.points))
  if pos==lastindex(hta.points)
    return scale(hta.points[pos])
  else
    along=x*lastindex(hta.points)-pos
    return scale(hta.points[pos]+along*(hta.points[pos+1]-hta.points[pos]))
  end
end

function htError(n::Int,hta::HalftoneApprox)
  z=scale(hta.points[n])
  area=areaBelow(x->ht(x,hta),z)
  z-area
end

function adjust!(hta::HalftoneApprox,n::Int)
  if n<=0 || n>=lastindex(hta.points)
    return
  end
  lo=invScale(lowerBound(oftype(hta.points[n],n)/lastindex(hta.points)))
  if n+1==lastindex(hta.points)
    hi=(hta.points[end]+hta.points[n])/2
  else
    hi=invScale(2*scale(hta.points[n+1])-scale(hta.points[n+2]))
  end
  hta.points[n]=lo
  loval=htError(n,hta)
  hta.points[n]=hi
  hival=htError(n,hta)
  mid=lo # arbitrary
  midval=loval
  sgn=sign(hival-loval)
  while midval!=0 && hi-lo>eps(hta.points[end])
    @assert hival*loval<=0
    mid=lo-loval*(hi-lo)/(hival-loval) # secant rule
    hta.points[n]=mid
    midval=htError(n,hta)
    #@printf "sec: %20.17f %20.17f %20.17f\n" lo mid hi
    if sgn*midval>0
      hi=mid
      hival=midval
    else
      lo=mid
      loval=midval
    end
    mid=(lo+hi)/2 # midpoint rule
    hta.points[n]=mid
    midval=htError(n,hta)
    #@printf "mid: %20.17f %20.17f %20.17f\n" lo mid hi
    if sgn*midval>0
      hi=mid
      hival=midval
    else
      lo=mid
      loval=midval
    end
  end
end

function adjust!(hta::HalftoneApprox)
  for n in reverse(eachindex(hta.points))
    adjust!(hta,n)
    @printf "\r%d " n
  end
end

end # module HalftoneFunction
