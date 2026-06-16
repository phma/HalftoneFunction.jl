module HalftoneFunction
using Integrals,Roots
export tryFunc1

# The halftone function h(x) is defined as follows:
# h(x) increases from 0 to 1.
# The area of {(x,y) in [0,1]×[0,1]: h(x)h(y)<z} is z whenever z is in [0,1].
# h(-x)=-h(x).
# h(1-x)=h(x).
# This module tries to compute an approximation to the halftone function.

function tryFunc1(x::AbstractFloat)
  (1-(x-1)^2)^(2/3)
end

end # module HalftoneFunction
