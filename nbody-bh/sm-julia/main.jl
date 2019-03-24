include("common.jl")
include("body.jl")
include("cell.jl")

using DelimitedFiles
using StaticArrays
using Base.Threads

const LEVEL = 2

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

function compute_leaf(leaf::Cell, bodies::Array{Body})
    for body in bodies
        if cell_contains_position(leaf, body.position)
            insert_body(leaf, body)
        end
    end
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
    n_leaves::Int = N_CELL_CHILDREN ^ LEVEL
    for iter in 1:iterations
        (universe_min, universe_max) = get_universe_size(bodies)

        for body in bodies
            reset_force(body)
        end

        root = Cell()
        root.min_bounds = universe_min
        root.max_bounds = universe_max

        construct_empty_tree(root, LEVEL)
        leaves = get_leaves(root, LEVEL)

        @threads for leaf in leaves
            compute_leaf(leaf, bodies)
        end

        update_empty_cells(root, LEVEL)

        @threads for body in bodies
            @inbounds compute_force(root, body)
        end

        @threads for body in bodies
            body.position = body.position .+ (dt .* body.velocity)
            body.velocity = body.velocity .+ (dt ./ body.mass .* body.force)
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

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
end
println(minimum(times))
