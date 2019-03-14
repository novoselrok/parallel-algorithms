include("common.jl")
include("body.jl")
include("worker.jl")

using DelimitedFiles
@everywhere using StaticArrays

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
        push!(bodies, Body(i, positions[i], velocities[i], masses[i], 1.0))
    end

    start_time = time_ns()

    world_size = nworkers()
    bodies_channels = [RemoteChannel(()->BodiesCh(world_size)) for _ in 1:world_size]
    bounds_channels = [RemoteChannel(()->BoundsCh(world_size)) for _ in 1:world_size]
    group_channels = [RemoteChannel(()->GroupCh(world_size)) for _ in 1:world_size]
    float_channels = [RemoteChannel(()->FloatCh(world_size)) for _ in 1:world_size]
    cell_channels = [RemoteChannel(()->CellTuplesCh(world_size)) for _ in 1:world_size]

    @sync for worker in workers()
        zb_rank = worker - 2
        ob_rank = worker - 1

        blocksize = div(n + world_size - 1, world_size)
        start = zb_rank * blocksize
        end_ = min(start + blocksize, n)
        start += 1
        subarray_size = end_ - start + 1
    
        my_bodies = bodies[start:end_]
        @spawnat worker work(my_bodies, bodies_channels, bounds_channels, group_channels, float_channels, cell_channels, iterations)
    end
    (time_ns() - start_time) / 1.0e9
end

println(ARGS)
nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    println(times)
    push!(times, main(ARGS))
end
println(minimum(times))
