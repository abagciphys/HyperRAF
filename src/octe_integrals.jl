###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                   TWO ELECTRON ATOMIC INTEGRALS                                         #
#                                  NONINTEGER SLATER-TYPE ORBITALS                                        #
###########################################################################################################
###########################################################################################################
###########################################################################################################
export OneCenterTwoERρ, OneCenterTwoERζ, OneCenterTwoER, OneCenterTwoE
###########################################################################################################
############################## ONE CENTER TWO ELECTRON BASIC INTEGRAL ####################################(43)
function OneCenterTwoERρ(mode::Symbol, L::Int, n1::arb, ρ1::arb, n2::arb, ρ2::arb)
    if mode == :test
        c1 = Gamma(n1 + n2 + 1) // Power(2 * (ρ1 + ρ2), n1 + n2 + 1)
        c2 = 1 // (n1 + L + 1)
        c3 = 1 // (n2 + L + 1)
    
        x1 = (2 * ρ1) // (2 * (ρ1 + ρ2))
        x2 = (2 * ρ2) // (2 * (ρ1 + ρ2))
    
        hyper1 = Hypergeometric2F1(RF(1), n1 + n2 + RF(1), n1 + L + RF(2), x1, sprec)
        hyper2 = Hypergeometric2F1(RF(1), n1 + n2 + RF(1), n2 + L + RF(2), x2, sprec)
    
        res1 = c2 * hyper1
        res2 = c3 * hyper2
    
        res = c1 * (res1 + res2)
    elseif mode == :use
        c1 = Gamma(n1 + n2 + 1) // Power(2 * (ρ1 + ρ2), n1 + n2 + 1)

        x1 = (2 * ρ1)
        x2 = (2 * ρ2)

        lowl = 0
        upl = L

        otwoerl = zeros(RF, upl - lowl + 1)
        otwoer2 = zeros(RF, upl - lowl + 1)

        otwoerl = HyperRL(:use, upl, n1, x1, n2, x2)
        otwoer2 = HyperRL(:use, upl, n2, x2, n1, x1)
    
        res = c1 .* (otwoerl + otwoer2)
        return res
    else
        error("Invalid mode specified. Please use either :test or :use.")
    end

end

function OneCenterTwoERζ(mode::Symbol, L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    return OneCenterTwoERρ(mode, L, n1, ζ1//2, n2, ζ2//2)
end
###########################################################################################################
################################# ONE CENTER TWO ELECTRON INTEGRALS #######################################(43)
function Slatern(n1::arb, n2::arb, ρ::arb, τ::arb)
    res1 = ((ρ*(One(τ)+τ))^(n1 + 1//2)) * ((ρ*(One(τ)-τ))^(n2 + 1//2))
    res2 = Sqrt(Gamma(2*n1 + One(n1)) * Gamma(2*n2 + One(n2)))
    res1 // res2
end
#################
function OneCenterTwoERTest(L::Int, n1::arb, ρ1::arb, τ1::arb, n2::arb, ρ2::arb, τ2::arb)
    lowl = 0
    upl = L

    res3 = zeros(RF, upl - lowl + 1)

    res1 = Slatern(n1, n2, ρ1, τ1)
    res2 = Slatern(n1, n2, ρ2, τ2)
    res3 = OneCenterTwoERρ(:use, upl, n1, ρ1, n2, ρ2)

    res =  res1 * res2 .* res3
end

function OneCenterTwoER(L::Int, 
    n1::arb, n1p::arb, ρ1::arb, τ1::arb,
    n2::arb, n2p::arb, ρ2::arb, τ2::arb
    )

    lowl = 0
    upl = L

    res3 = zeros(RF, upl - lowl + 1)

    res1 = Slatern(n1, n1p, 2 * ρ1, τ1)
    res2 = Slatern(n2, n2p, 2 * ρ2, τ2)
    res3 = OneCenterTwoERρ(:use, upl, n1 + n1p, ρ1, n2 + n2p, ρ2)
    res = res1 * res2 .* res3
end

function OneCenterTwoE(
    n1::arb, l1::Int, m1::Int, 
    n1p::arb, l1p::Int, m1p::Int,
    ρ1::arb, τ1::arb,
    n2::arb, l2::Int, m2::Int, 
    n2p::arb, l2p::Int, m2p::Int, 
    ρ2::arb, τ2::arb)

    if Abs(m1) > Abs(l1) || Abs(m1p) > Abs(l1p) || Abs(m2) > Abs(l2) || Abs(m2p) > Abs(l2p)
        error(InvalidQuantumNumbersError("Invalid quantum numbers! l must be greater than or equal to |m|."))
        return RF(0)
    else
        lowl = max(abs(l1-l1p), abs(l2-l2p))
        upl = min(l1 + l1p, l2 + l2p)

        gaunt11p = zeros(RF, upl - lowl + 1)
        gaunt22p = zeros(RF, upl - lowl + 1)
        onectwoer = zeros(CF, upl - lowl + 1)

        onectwoer = OneCenterTwoER(upl, n1, n1p, ρ1, τ1, n2, n2p, ρ2, τ2)
        onectwoe = RF(0)

        for s in lowl : upl
            s1 = s + 1
            gaunt11p[s1] = GGauntG(l1, m1, l1p, m1p, s, m1 - m1p)
            gaunt22p[s1] = GGauntG(l2, m2, l2p, m2p, s, m2 - m2p)

            onectwoe += gaunt11p[s1] * gaunt22p[s1] * onectwoer[s1]
        end 
    end
    if isnan(NO(onectwoe))
        return RF(0)
    else
        return onectwoe
    end
end