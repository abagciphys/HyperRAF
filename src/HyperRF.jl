#module HyperRF

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

include("sto_at_integ_two_elect.jl")

#end
