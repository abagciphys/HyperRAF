using HyperRAF
using Test

@testset "HyperRF.jl" begin
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEeL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10)))) 
    == 0.6145299803434426)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEgL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -31630.170839641007)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEhL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -6.374045915920716)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoElL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 31618.898932025255)
    @test (Float64(NO(HypREval(HyperRadialF(OTwoEmL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == -6.92692516554162)
    @test (Float64(NO(HypREval(HyperRadialF(HyperRLTest, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
    == 0.6145299803434426)
#     @test Float64(NO(HypREval(HyperRadialF(OTwoEeL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
#     == 0.6145299803434426
#     @test Float64(NO(HypREval(HyperRadialF(OTwoEeL, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))))
#     == 0.6145299803434426
end
