module SOCRATESSingleColumnForcingsNonNegLeastSquaresExt

using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using NonNegLeastSquares: NonNegLeastSquares

SSCF.Interpolation.nnls_solve(A, b; kwargs...) = NonNegLeastSquares.nonneg_lsq(A, b; kwargs...)

end # module
