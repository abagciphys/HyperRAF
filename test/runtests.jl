using HyperRF
using Test

@testset "HyperRF.jl" begin
    HypREval(HyperRadialF(OTwoEeL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(OTwoEgL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(OTwoEhL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(OTwoElL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(HyperRLTest, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(HyperRL0, 1), RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(HyperRL1, 1), RF(15//10),RF(5//10),RF(16//10),RF(6//10))
    HypREval(HyperRadialF(OTwoERL, 1), RF(15//10),RF(5//10),RF(16//10),RF(6//10))
end
