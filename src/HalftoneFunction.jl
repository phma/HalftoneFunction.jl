module HalftoneFunction
using Integrals,Roots
export tryFunc1

function tryFunc1(x::AbstractFloat)
  (1-(x-1)^2)^(2/3)
end

end # module HalftoneFunction
