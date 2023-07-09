###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                  HYPERGEOMETRIC RADIALS FUNCTIONS                                       #
###########################################################################################################
###########################################################################################################
###########################################################################################################
export HyperRadialF, HypREval, OTwoEeL, OTwoEfL, OTwoEgL, OTwoEhL, OTwoElL, OTwoEmL, HyperRLTest, HyperRL0,
HyperRL1, HyperRL
###########################################################################################################
######################################## TRANSLATION STRUCTURE ############################################
struct HyperRadialF
    HypR::Function
    i::Int
end
#################
function HypREval(HypR::HyperRadialF, L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if HypR.i == 1
        return HypR.HypR(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    elseif HypR.i == 2
        return HypR.HypR(L::Int, n2::arb, ζ2::arb, n1::arb, ζ1::arb)
    else
        error("Invalid indicator value. Expected 1 or 2.")
    end
end

function HypREval(HypR::HyperRadialF, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if HypR.i == 1
        return HypR.HypR(n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    elseif HypR.i == 2
        return HypR.HypR(n2::arb, ζ2::arb, n1::arb, ζ1::arb)
    else
        error("Invalid indicator value. Expected 1 or 2.")
    end
end
####Evaluation example
#HypREval(HyperRadialF(HyperRLTest, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
###########################################################################################################
########################################## AUXILIARY FUNCTIONS ############################################
################# FOR HYPERGEOMETRIC RADIALS FUNCTIONS
function OTwoEeL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = Gamma(n1 + n2 + RF(1))
    res2 = Power(ζ1 + ζ2, n1 + n2 + RF(1))
    res3 = n1 + L + RF(1)

    res = res1 // (res2 * res3)
end
#################
function OTwoEfL(L::Int, n1::arb, n2::arb)
    res1 = RF(pi) * Csc((-n2 + L) * RF(pi))
    res2 = Gamma(-n2 + L + RF(1))
    res3 = (n1 + L + RF(1))
    res4 = Gamma(n2 - L + RF(1))

    res = (res1 // res2) * (res3 // res4)
end
#################
function OTwoEgL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = RF(pi) * Csc((-n2 + L) * RF(pi))
    res2 = Gamma(-n2 + L + RF(1))
    res3 = Power(ζ1 // (ζ1 + ζ2), -n1 - L - RF(1))
    res4 = Power(ζ2 // (ζ1 + ζ2), -n2 + L)
    res5 = Gamma(n1 + L + 2)
    res6 = Gamma(n1 + n2 + RF(1))

    res = (res1 // res2) * res3 * res4 * (res5 // res6)
end
#################
function OTwoEhL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = (n2 + L + RF(1)) // (n1 + L + RF(1))
    res2 = Pochhammer(n1 - L, 2 * L + RF(1))
    res3 = Pochhammer(-n2 - L - RF(1), 2 * L + RF(1))
    res4 = Power(-(ζ2 // ζ1), 2 * L + RF(1))
    res5 = OTwoEfL(L, n1, n2)

    res = res1 * (res2 // res3) * res4 * res5 + RF(1)
end
#################
function OTwoElL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = (RF(pi) * Csc((-n2 + L) * RF(pi))) // (Gamma(-n2 + L + RF(1)))
    res2 = (n1 + L + RF(1)) // (Gamma(n2 - L + RF(1)))
    res3 = Power(-(ζ2 // ζ1), 2 * L + RF(1))
    res4 = (ζ1 + ζ2) // (ζ2)

    lowl = 1
    upl1 = 2 * L + 1
    upl2 = upl1 - 2

    al1 = zeros(RF, upl1 - lowl + 1)
    bl1 = zeros(RF, upl1 - lowl - 1)
    cl1 = zeros(RF, upl1 - lowl - 1)

    al2 = zeros(RF, upl1 - lowl + 1)

    al1[1] = Pochhammer(n1 - L + RF(1), 2 * L)
    al1[2] = Pochhammer(n1 - L + RF(2), 2 * L - RF(1))
    al2[2 * L + 1] = RF(1)
    al2[2 * L] = (n2 - L + RF(1))
    
    for s in lowl : upl2
        s1 = upl2 + lowl - s
        cl1[s] = RF(1) // ((n1 - L + s) * (n1 - L + s + RF(1))) 
        bl1[s] = cl1[s] * al1[s]
        al1[s + 2] = bl1[s]

        al2[s1] = (n2 + L + RF(1) - s1) * al2[s1 + 1]
    end
    res8 = RF(0)
    for s in lowl : upl1
        res5 = al1[s]
        res6 = Power(-RF(1), 2 * L + 1 - s) * al2[s]
        res7 = Power(-ζ1 // ζ2, s - 1)
        res8 += (res5 // res6) * res7
    end
    res = res1 * res2 * res3 * res4 * res8 
end
#################
function OTwoEmL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = OTwoEeL(L, n1, ζ1, n2, ζ2)
    res2 = OTwoEgL(L, n1, ζ1, n2, ζ2)
    res3 = OTwoElL(L, n1, ζ1, n2, ζ2)

    res = res1 * (res2 + res3)
end
###########################################################################################################
################################### HYPERGEOMETRIC RADIALS FUNCTIONS ######################################(43)
################# RECURRENCE RELATION TECHNIQUE
function HyperRLTest(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    ### L>1
    res1 = OneCenterTwoERζ(:test, L, n1, ζ1, n2, ζ2)
    res2 = OTwoEmL(L, n2, ζ2, n1, ζ1)
    res3 = OTwoEeL(L, n1, ζ1, n2, ζ2)
    res4 = OTwoEhL(L, n2, ζ2, n1, ζ1)

    res = (res1 + res2) // (res3 * res4)
end

function HyperRL0(n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = (n1 + RF(1)) * Power(ζ1 // (ζ1 + ζ2), -n1 - RF(1))
    res2 = Power(ζ2 // (ζ1 + ζ2), -n2)
    res3 = (Gamma(n1 + RF(1)) * Gamma(n2)) // (Gamma(n1 + n2 + RF(1)))
    res4 = Beta(n2, n1 + RF(1), ζ2 // (ζ1 + ζ2))

    res = res1 * res2 * (res3 - res4)
end
#################
function HyperRL1(n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = (n1 + RF(2)) // (n2 - RF(1))
    res2 = ζ2 // ζ1
    res3 = HyperRL0(n1, ζ1, n2, ζ2)
    res4 = (ζ1 + ζ2) // ζ1

    res = res1 * (res2 * res3 - res4)
end
#################
function HyperRL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if L == 0
        res = HyperRL0(n1, ζ1, n2, ζ2)
    elseif L == 1
        res =  HyperRL1(n1, ζ1, n2, ζ2)
    else
        lowl = 0
        upl = L

        otwoerl = zeros(RF, upl - lowl + 1)

        otwoerl[1] = HyperRL0(n1, ζ1, n2, ζ2)
        otwoerl[2] = HyperRL1(n1, ζ1, n2, ζ2)

        for s in lowl : upl - 2
            s1 = s - lowl + 1
            res1 = (n1 + s + 3) // (ζ1 * (n1 + s + 2) * (-n2 + s + 2))
            res2 = ζ2 * (n1 + s + 2)
            res3 = ζ1 * (-n2 + s + 1) - ζ2 * (n1 + s + 2)

            otwoerl[s1 + 2] = res1 * (res2  * otwoerl[s1] + res3 * otwoerl[s1 + 1])
        end
        res = otwoerl[L + 1]
    end
    return res
end