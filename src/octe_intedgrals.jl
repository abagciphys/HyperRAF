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
        c2 = 1 // (n1 + L + 1)
        c3 = 1 // (n2 + L + 1)

        x1 = (2 * ρ1)
        x2 = (2 * ρ2)
    
        hyper1 = OTwoERL(L, n1, x1, n2, x2)
        hyper2 = OTwoERL(L, n2, x2, n1, x1)
    
        res1 = c2 * hyper1
        res2 = c3 * hyper2
    
        res = c1 * (res1 + res2)
    else
        error("Invalid mode specified. Please use either :test or :use.")
    end

end

function OneCenterTwoERζ(mode::Symbol, L::Int, n1::arb, ζ1::arb, n2::arb, ζ2::arb)
    return OneCenterTwoERρ(mode, L, n1, ζ1//2, n2, ζ2//2)
end
###########################################################################################################
################################# ONE CENTER TWO ELECTRON INTEGRALS #######################################(43)
function OneCenterTwoER(L::Int, n1::arb, n2::arb, ρ1::arb, τ1::arb, ρ2::arb, τ2::arb)
    res1 = SlaterON(n1, n2, ρ1, τ1)
    res2 = SlaterON(n1, n2, ρ2, τ2)
    res3 = OneCenterTwoER(L, n1, n2, ρ1, ρ2)

    res = res1 * res2 * res3
end

function OneCenterTwoE(
    n1::arb, l1::Int, m1::Int, ρ1::arb, τ1::arb, 
    n1p::arb, l1p::Int, m1p::Int, ρ1p::arb, τ1p::arb, 
    n2::arb, l2::Int, m2::Int, ρ2::arb, τ2::arb, 
    n2p::arb, l2p::Int, m2p::Int, ρ2p::arb, τ2p::arb)

    if Abs(m1) > Abs(l1) || Abs(m1p) > Abs(l1p) || Abs(m2) > Abs(l2) || Abs(m2p) > Abs(l2p)
        res = 0
    else
        lowl = max(abs(l1-l1p), abs(l2-l2p))
        upl = min(l1 + l1p, l2 + l2p)

        gaunt11p = zeros(RF, upl - lowl + 1)
        gaunt22p = zeros(RF, upl - lowl + 1)
        onectwoer = zeros(CF, upl - lowl + 1)
        onectwoe = zeros(CF, upl - lowl + 1)

        for s in lowl : upl
            s1 = s + 1
            gaunt11p[s1] = GGauntG(l1, m1, l1p, m1p, s, m1 - m1p)
            gaunt22p[s1] = GGauntG(l2, m2, l2p, m2p, s, m2 - m2p)
            onectwoer[s1] = OneCenterTwoER(s, n1, n2, ρ1, τ1, ρ2, τ2)

            onectwoe[s1] = gaunt11p[s1] * gaunt22p[s1] * onectwoer[s1]
        end 
    end
    return onectwoe
end