include("common.jl")
include("body.jl")
include("cell.jl")

using DelimitedFiles
using StaticArrays

function get_universe_size(bodies::Array{Body})::Tuple{Vec3, Vec3}
    universe_min = [Inf, Inf, Inf]
    universe_max = [-Inf, -Inf, -Inf]

    for body in bodies
        for c in 1:DIMS
            if body.position[c] < universe_min[c]
                universe_min[c] = body.position[c]
            end
            if body.position[c] > universe_max[c]
                universe_max[c] = body.position[c]
            end
        end
    end

    (Vec3(universe_min), Vec3(universe_max))
end

function main(args)
    fname = args[1]
    iterations = parse(Int64, args[2])

    data = readdlm(fname)::Array{Float64, 2}
    n = size(data, 1)
    positions = [Vec3(data[i, 1:3]) for i in 1:n]
    velocities = [Vec3(data[i, 4:6]) for i in 1:n]
    masses::Array{Float64} = data[:, 7]
    bodies::Array{Body} = []

    for i in 1:n
        push!(bodies, Body(i, positions[i], velocities[i], masses[i]))
    end

    start = time_ns()
    dt = 0.1
    for iter in 1:iterations
        (universe_min, universe_max) = get_universe_size(bodies)

        @time for body in bodies
            reset_force(body)
        end

        root = Cell()
        root.min_bounds = universe_min
        root.max_bounds = universe_max

        @time for body in bodies
            insert_body(root, body)
        end
        
        @time for body in bodies
            compute_force(root, body)
        end

        @time for body in bodies
            body.position[X] += dt * body.velocity[X];
            body.position[Y] += dt * body.velocity[Y];
            body.position[Z] += dt * body.velocity[Z];

            body.velocity[X] += dt / body.mass * body.force[X];
            body.velocity[Y] += dt / body.mass * body.force[Y];
            body.velocity[Z] += dt / body.mass * body.force[Z];
        end
    end
    elapsed = (time_ns() - start) / 1.0e9

    final_positions::Array{Float64, 2} = zeros(n, 3)
    final_velocities::Array{Float64, 2} = zeros(n, 3)
    for i in 1:n
        body = bodies[i]
        final_positions[i, :] = body.position
        final_velocities[i, :] = body.velocity
    end
    writedlm("out.txt", hcat(final_positions, final_velocities))
    
    elapsed
end

# println(ARGS)
# nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
# times = []
# for i in 1:nprecompilesteps
#     println(times)
#     push!(times, main(ARGS))
# end
# println(minimum(times))
