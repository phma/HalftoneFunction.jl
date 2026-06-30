module HalftoneFunction
using QuadGK,Roots,OffsetArrays,Printf,CairoMakie
export HalftoneApprox,ht,adjust!,plotHalftoneFunction
export writeHalftone,writeHalftones,readHalftone,readHalftones,findHta
export marshal,unmarshal,readHeader,writeHeader

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

const SCALE_ID=0x0001

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

"""
    areaBelow(h::Function,z::AbstractFloat)

Compute the area of the unit square in which `h(x)*h(y)` is below `z`.
"""
function areaBelow(h::Function,z::AbstractFloat)
  diagCross=find_zero(x->h(x)-√z,[0,1])
  integrand=y->find_zero(x->h(x)-z/h(y),[0,1])
  quadgk(integrand,diagCross,1)[1]*2+diagCross^2
end

function newPoints(T::DataType,n::Integer)
  ret=OffsetVector(T[],-1)
  for i in 0:n
    push!(ret,invScale(lowerBound((T(i)/n)^2)))
  end
  ret
end

"""
    struct HalftoneApprox{T}

Holds the points of an approximation to the halftone function. Construct with
`hta=HalftoneApprox(T,n)` where `T` is a type such as `Float64` and `n` is
one more than the number of points.

See also `adjust!`.
"""
struct HalftoneApprox{T}
  points	::OffsetVector{T}
  HalftoneApprox(T::DataType,n::Integer)=new{T}(newPoints(T,n))
end

const htaList=Dict{Tuple{DataType,Int},HalftoneApprox}()

"""
    ht(x::T,hta::HalftoneApprox{T})

Compute the halftone function at `x` in [0,1], using `hta`.
"""
function ht(x::T,hta::HalftoneApprox{T}) where T<:AbstractFloat
  if x<0 || x>1
    throw(DomainError(x,"The argument to ht must be between 0 and 1 inclusive."))
  end
  pos=floor(Int,√(x*lastindex(hta.points)^2))
  if pos==lastindex(hta.points)
    return scale(hta.points[pos])
  else
    along=(x*lastindex(hta.points)^2-pos^2)/(2*pos+1)
    return scale(hta.points[pos]+along*(hta.points[pos+1]-hta.points[pos]))
  end
end

"""
    htError(n::Int,hta::HalftoneApprox)

Compute the error in the `n`th point in `hta`. Except for the next-to-last
point, this is not the error in the point's elevation, but the difference
between the elevation and the area below that elevation. For the next-to-last
point, the area is independent of the elevation and is simply the area outside
the upper right quarter circle.
"""
function htError(n::Int,hta::HalftoneApprox)
  z=scale(hta.points[n])
  area=areaBelow(x->ht(x,hta),z)
  z-area
end

"""
    adjust!(hta::HalftoneApprox,n::Int)

Adjust the `n`th point of `hta` so that its `htError` is zero. All points after
the `n`th must already be adjusted.
"""
function adjust!(hta::HalftoneApprox,n::Int)
  if n<=0 || n>=lastindex(hta.points)
    return
  end
  lo=invScale(lowerBound(oftype(hta.points[n],n)^2/lastindex(hta.points)^2))
  if n+1==lastindex(hta.points)
    hi=(hta.points[end]+hta.points[n])/2
  else
    hi=invScale((4*(n+1)*scale(hta.points[n+1])-(2n+1)*scale(hta.points[n+2]))/(2n+3))
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

"""
    adjust!(hta::HalftoneApprox)

Adjust `hta` so that `ht(x,hta)` is a good approximation to the halftone function.
Time to adjust increases as the 1.5th power of the number of points.
"""
function adjust!(hta::HalftoneApprox)
  for n in reverse(eachindex(hta.points))
    adjust!(hta,n)
    @printf "\r%d " n
  end
end

function findHta(T::DataType,n::Integer)
  if haskey(htaList,(T,n))
    htaList[(T,n)]
  else
    hta=HalftoneApprox(T,n)
    adjust!(hta)
    htaList[(T,n)]=hta
    hta
  end
end

"""
    marshal(p::AbstractFloat)

Convert `p` to a sequence of bytes, which can be converted back to `p` given
its type. `p` must be in [0,1).
"""
function marshal(p::AbstractFloat)
  ret=[0x00,0x00]
  while p>0
    p*=256
    push!(ret,floor(UInt8,p))
    p=modf(p)[1]
  end
  ret[2]=(length(ret)-2)%256
  ret[1]=(length(ret)-2)÷256
  ret
end

"""
    unmarshal(T::DataType,bytes::Vector{UInt8})

Convert `bytes` to a floating-point number of type `T`. `bytes` has had the two
length bytes stripped off when read from a file.
"""
function unmarshal(T::DataType,bytes::Vector{UInt8})
  ret=zero(T)
  for b in reverse(bytes)
    ret=(b+ret)/256
  end
  ret
end

function typeName(T::DataType)
  ret=string(T)
  if findlast('.',ret)==nothing
    ret="Base."*ret
  end
  ret
end

function nameType(name::String)
  pos=findlast('.',name)
  mod=name[1:pos-1]
  tname=name[pos+1:end]
  if mod=="Base"
    mod=Base
  elseif mod=="Core"
    mod=Core
  end
  getfield(mod,Symbol(tname))
end

function writeHeader(f::IO)
  write(f,marshal(√π/2))
  # magic number expressing how fast halftone circles grow compared to squares
  write(f,hton(SCALE_ID))
end

function readFloat(f::IO,T::DataType)
  lenBytes=read(f,2)
  if length(lenBytes)==2
    len=lenBytes[1]*256+lenBytes[2]
  else
    len=0
  end
  fpBytes=read(f,len)
  if length(lenBytes)<2 || length(fpBytes)<len
    error("end of file reading floating point number")
  end
  unmarshal(T,fpBytes)
end

"""
    readHeader(f::IO)

Read the header of a halftone function file, returning the number identifying the scale
function used when it was written.
"""
function readHeader(f::IO)
  magic=readFloat(f,Float64)
  if magic==√π/2
    ntoh(read(f,UInt16))
  else
    throw(ErrorException("invalid header"))
  end
end

function writeHalftone(f::IO,hta::HalftoneApprox)
  write(f,typeName(eltype(hta.points)))
  write(f,0x00)
  write(f,hton(UInt32(lastindex(hta.points))))
  for i in eachindex(hta.points)
    if i>0 && i<lastindex(hta.points) # these points are 0 and 1, no need to write them
      write(f,marshal(hta.points[i]))
    end
  end
end

function readHalftone(f::IO)
  if eof(f)
    return nothing
  end
  typestr=String(readuntil(f,0x00))
  T=nameType(typestr)
  if !(T<:AbstractFloat)
    throw(ErrorException("invalid type"))
  end
  len=ntoh(read(f,UInt32))
  ret=HalftoneApprox(T,len)
  for i in 1:len-1
    ret.points[i]=readFloat(f,T)
    if ret.points[i]<ret.points[i-1]
      throw(ErrorException("points out of order"))
    end
  end
  ret
end

function writeHalftones(f::IO,hta::Vector{<:HalftoneApprox})
  for h in hta
    writeHalftone(f,h)
  end
end

function readHalftones(f::IO)
  ret=HalftoneApprox[]
  while true
    ht=readHalftone(f)
    if ht==nothing
      break
    end
    push!(ret,ht)
  end
  ret
end

function writeHalftones(filename::String,hta::Vector{<:HalftoneApprox})
  f=open(filename,"w")
  writeHeader(f)
  writeHalftones(f,hta)
  close(f)
end

function readHalftones(filename::String)
  f=open(filename,"r")
  try
    scaleId=readHeader(f)
    ret=readHalftones(f)
    # TODO rescale them if scaleId!=SCALE_ID
    return ret
  catch
    rethrow()
  finally
    close(f)
  end
end

function plotHalftoneFunction()
  hta34=findHta(Float64,34)
  htf=Figure(size=(1189,841))
  htfax=Axis(htf[1,1])
  x=(0:34*21)./34/21
  y=(x->ht(x,hta34)).(x)
  lines!(htfax,x,y)
  save("halftone-function.svg",htf)
end

end # module HalftoneFunction
