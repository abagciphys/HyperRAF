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

function HyperRaFT(HyperD::HyperRaFD, mode::Symbol, L::Int, n1::arb, ρ1::arb, n2::arb, ρ2::arb)
    user_home_dir = homedir()
    file_name = "cpu_$(HyperD.mode)_time.dat"
    file_path = joinpath(user_home_dir, file_name)
    file = open(file_path, "w")

    total_time = 0.0
    for s in 1:L
        time = @elapsed HyperD.HyperD(:mode, L, n1, ρ1, n2, ρ2)
        total_time += time
        println(file, "$L\t$time\t$total_time")  # Write time and total_time in two separate columns
    end

    close(file)
end