using DelimitedFiles
using StaticArrays

const G = 6.67e-11
const dt = 0.1

function main()
    fname = ARGS[1]
    iterations = parse(Int64, ARGS[2])

    data = readdlm(fname)::Array{Float64, 2}
    n = size(data, 1)
    positions = [SVector{3, Float64}(data[i, 1:3]) for i in 1:n]
    velocities = [SVector{3, Float64}(data[i, 4:6]) for i in 1:n]
    masses = data[:, 7]

    start = time_ns()
    for iter in 1:iterations
        # zero(type) returns additive identity
        # zero(SVector{3, Float64}) -> static array with three zeroes
        forces_per_thread::Array{Array{SVector{3, Float64}, 1}, 1} = [fill(zero(SVector{3, Float64}), n) for _ in 1:Threads.nthreads()]
        Threads.@threads for q in 1:n
            thread_id = Threads.threadid()
            @inbounds for k in 1:n
                k <= q && continue

                diff = positions[q] - positions[k]
                dist = sqrt(sum(diff .^ 2))
                dist_cubed = dist ^ 3

                tmp = -G * masses[q] * masses[k] / dist_cubed
                force_qk = tmp * diff
                forces_per_thread[thread_id][q] += force_qk
                forces_per_thread[thread_id][k] -= force_qk
            end
        end
        # Reduction step
        sum_forces = sum(forces_per_thread)

        Threads.@threads for q in 1:n
            @inbounds positions[q] += (dt * velocities[q])
            @inbounds velocities[q] += (dt / masses[q] * sum_forces[q])
        end
    end
    println((time_ns() - start) / 1.0e9)
    # writedlm("out.txt", hcat(positions, velocities))
end

main()
