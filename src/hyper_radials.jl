###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                  HYPERGEOMETRIC RADIALS FUNCTIONS                                       #
###########################################################################################################
###########################################################################################################
###########################################################################################################
export HyperRaF, HypREval, HyperReL, HyperRfL, HyperRgL, HyperRhL, HyperRlL, HyperRmL, 
HyperRLTest, HyperRL0, HyperRL1, HyperRL
###########################################################################################################
######################################## TRANSLATION STRUCTURE ############################################
struct HyperRaF
    HyperR::Function
    i::Int
end
#################
function HypREval(HyperR::HyperRaF, L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if HyperR.i == 1
        return HyperR.HyperR(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    elseif HyperR.i == 2
        return HyperR.HyperR(L::Int, n2::arb, ζ2::arb, n1::arb, ζ1::arb)
    else
        error("Invalid indicator value. Expected 1 or 2.")
    end
end

function HypREval(HyperR::HyperRaF, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if HyperR.i == 1
        return HyperR.HyperR(n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    elseif HyperR.i == 2
        return HyperR.HyperR(n2::arb, ζ2::arb, n1::arb, ζ1::arb)
    else
        error("Invalid indicator value. Expected 1 or 2.")
    end
end
####Evaluation example
#HypREval(HyperRadialF(HyperRLTest, 1), 5,RF(15//10),RF(5//10),RF(16//10),RF(6//10))
###########################################################################################################
########################################## AUXILIARY FUNCTIONS ############################################
################# FOR HYPERGEOMETRIC RADIALS FUNCTIONS
function HyperReL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = Gamma(n1 + n2 + RF(1))
    res2 = Power(ζ1 + ζ2, n1 + n2 + RF(1))
    res3 = n1 + L + RF(1)

    res = res1 // (res2 * res3)
end
#################
function HyperRfL(L::Int, n1::arb, n2::arb)
    res1 = RF(pi) * Csc((-n2 + L) * RF(pi))
    res2 = Gamma(-n2 + L + RF(1))
    res3 = (n1 + L + RF(1))
    res4 = Gamma(n2 - L + RF(1))

    res = (res1 // res2) * (res3 // res4)
end
#################
function HyperRgL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = RF(pi) * Csc((-n2 + L) * RF(pi))
    res2 = Gamma(-n2 + L + RF(1))
    res3 = Power(ζ1 // (ζ1 + ζ2), -n1 - L - RF(1))
    res4 = Power(ζ2 // (ζ1 + ζ2), -n2 + L)
    res5 = Gamma(n1 + L + 2)
    res6 = Gamma(n1 + n2 + RF(1))

    res = (res1 // res2) * res3 * res4 * (res5 // res6)
end
#################
function HyperRhL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = (n2 + L + RF(1)) // (n1 + L + RF(1))
    res2 = Pochhammer(n1 - L, 2 * L + RF(1))
    res3 = Pochhammer(-n2 - L - RF(1), 2 * L + RF(1))
    res4 = Power(-(ζ2 // ζ1), 2 * L + RF(1))
    res5 = HyperRfL(L, n1, n2)

    res = res1 * (res2 // res3) * res4 * res5 + RF(1)
end
#################
function HyperRlL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
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
function HyperRmL(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    res1 = HyperReL(L, n1, ζ1, n2, ζ2)
    res2 = HyperRgL(L, n1, ζ1, n2, ζ2)
    res3 = HyperRlL(L, n1, ζ1, n2, ζ2)

    res = res1 * (res2 + res3)
end
###########################################################################################################
################################### HYPERGEOMETRIC RADIALS FUNCTIONS ######################################(43)
################# RECURRENCE RELATION TECHNIQUE
function HyperRLTest(L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    ### L>1
    res1 = OneCenterTwoERζ(:test, L, n1, ζ1, n2, ζ2)
    res2 = HyperRmL(L, n2, ζ2, n1, ζ1)
    res3 = HyperReL(L, n1, ζ1, n2, ζ2)
    res4 = HyperRhL(L, n2, ζ2, n1, ζ1)

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
function HyperRL(mode::Symbol, L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    if L == 0
        hyperrl = [HyperRL0(n1, ζ1, n2, ζ2)]
        otwoerl = [HyperRL0(n1, ζ1, n2, ζ2) // (n1 + 0 + 1)]
    elseif L == 1
        hyperrl =  [HyperRL1(n1, ζ1, n2, ζ2)]
        otwoerl =  [HyperRL1(n1, ζ1, n2, ζ2) // (n1 + 1 + 1)]
    else
        lowl = 0
        upl = L

        hyperrl = zeros(RF, upl - lowl + 1)
        otwoerl = zeros(RF, upl - lowl + 1)
        #otwoerl is for use in two-electron intgrals(memorizing the coefficients already)

        hyperrl[1] = HyperRL0(n1, ζ1, n2, ζ2)
        otwoerl[1] = HyperRL0(n1, ζ1, n2, ζ2) // (n1 + 0 + 1)

        hyperrl[2] = HyperRL1(n1, ζ1, n2, ζ2)
        otwoerl[2] = HyperRL1(n1, ζ1, n2, ζ2) // (n1 + 1 + 1)

        for s in lowl : upl - 2
            s1 = s - lowl + 1
            res1 = (n1 + s + 3) // (ζ1 * (n1 + s + 2) * (-n2 + s + 2))
            res2 = ζ2 * (n1 + s + 2)
            res3 = ζ1 * (-n2 + s + 1) - ζ2 * (n1 + s + 2)

            hyperrl[s1 + 2] = res1 * (res2  * hyperrl[s1] + res3 * hyperrl[s1 + 1])
            otwoerl[s1 + 2] = hyperrl[s1 + 2] // (n1 + s + 2 + 1)
        end
        if mode == :test
            return hyperrl
        elseif mode == :use
            return otwoerl
        else
            error("Invalid mode specified. Please use either :test or :use.")
        end
    end
end