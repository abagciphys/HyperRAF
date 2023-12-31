###########################################################################################################
###########################################################################################################
###########################################################################################################
#                                               BENCHMARKING                                              #
###########################################################################################################
###########################################################################################################
###########################################################################################################
export HyperRaFD, HyperRaFT
###########################################################################################################
struct HyperRaFD
    HyperD::Function
    mode::Symbol
end

function HyperRaFT(HyperD::HyperRaFD, L::Int, n1::arb, ρ1::arb, n2::arb, ρ2::arb)
    user_home_dir = homedir()
    file_time = "cpu_$(HyperD.mode)_jtime.dat"
    file_val = "cpu_$(HyperD.mode)_jval.dat"
    file_tpath = joinpath(user_home_dir, file_time)
    file_vpath = joinpath(user_home_dir, file_val)
    filet = open(file_tpath, "w")
    filev = open(file_vpath, "w")

    total_timet1 = 0.0
    total_timet2 = 0.0
    total_timeu = 0.0
    if HyperD.mode == :test
        for s in 0 : L
            for s1 in 0 : s
                time = Float32(@elapsed HyperD.HyperD(HyperD.mode, s1, n1, ρ1, n2, ρ2))
                total_timet1 += time
            end
            total_timet2 += total_timet1
            println(filet, "$s\t$total_timet1\t$total_timet2")
            val = HyperD.HyperD(HyperD.mode, s, n1, ρ1, n2, ρ2)
            println(filev, "$s\t$val")
        end
    elseif HyperD.mode == :use
        val = HyperD.HyperD(HyperD.mode, L, n1, ρ1, n2, ρ2)
        for (index1, value) in enumerate(val)
            index2 = index1 - 1
            println(filev, "$index2\t$value")
        end
        for s in 0 : L
            time = Float32(@elapsed HyperD.HyperD(HyperD.mode, s, n1, ρ1, n2, ρ2))
            println(filet, "$s\t$time\t$total_timeu")
            total_timeu += time
        end 
    end
    close(filet)
    close(filev)
end
############## Example of using
#HyperRaFT(HyperRaFD(OneCenterTwoERρ, :test), 50, RF(15//10),RF(5//10),RF(16//10),RF(6//10))
