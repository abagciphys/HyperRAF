module HyperRAF

using Nemo
using JRAF

sprec = 750;
ARBF = ArbField(sprec);
ACBF = AcbField(sprec);
ComplexField(prec::Int) = AcbField(prec);
CF = ComplexField(sprec);
RealField(prec::Int) = ArbField(prec);
RF = RealField(sprec);
export sprec, RF, CF

include("hyper_radials.jl")
include("octe_integrals.jl")

end
