using HyperRAF
using Test

L=5
n1 = RF(15//10)
n2 = RF(16//10)
ζ1 = RF(5//10)
ζ2 = RF(6//10)
x1 = ζ1 // (ζ1 + ζ2)
x2 = ζ2 // (ζ1 + ζ2)
sprec = 750

lefthyper1 = Float64(NO(HypREval(HyperRadialF(HyperRL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
righthyper1= Float64(NO(Hypergeometric2F1(RF(1), n1 + n2 + RF(1), n1 + L + RF(2), x1, sprec)))

lefthyper2 = Float64(NO(HypREval(HyperRadialF(HyperRL, 2), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
righthyper2= Float64(NO(Hypergeometric2F1(RF(1), n1 + n2 + RF(1), n2 + L + RF(2), x2, sprec)))

@testset "HyperRF.jl" begin
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEeL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10)))) 
    == 0.6145299803434426)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEgL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -31630.170839641007)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEhL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -6.374045915920716)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoElL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 31618.898932025255)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEmL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -6.92692516554162)
    @test (Float64(NO(HypREval(HyperRadialF(HyperRLTest, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 1.2925402895942881)
    @test (Float64(NO(HypREval(HyperRadialF(HyperRL0, 1), RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 2.0790239933881467)
    @test (Float64(NO(HypREval(HyperRadialF(HyperRL1, 1), RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 1.719834620383692)
    @test (Float64(NO(HypREval(HyperRadialF(HyperRL, 1), 5, RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 1.2925402895942881)
    @test lefthyper1 == righthyper1
    @test lefthyper2 == righthyper2
end