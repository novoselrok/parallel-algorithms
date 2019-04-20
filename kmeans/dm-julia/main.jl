using Distributed
using DelimitedFiles
@everywhere using LinearAlgebra
@everywhere using StaticArrays
@everywhere using DistributedArrays

@everywhere import Base.+

@everywhere const POINT_SIZE = 128
# Type alias for Point type
@everywhere const Point = SVector{POINT_SIZE, Float64}

# Datatypes
@everywhere mutable struct Cluster
   size::Int64
   point_sum::Point
   mean::Point
end
@everywhere const ClustersCh = Channel{Array{Cluster}}

# Constructor
@everywhere Cluster() = Cluster(0, Point(zeros(POINT_SIZE)), Point(zeros(POINT_SIZE)))
@everywhere Cluster(point::Point) = Cluster(0, Point(zeros(POINT_SIZE)), point)

@everywhere function distance(cluster::Cluster, point::Point)
    diff = cluster.mean .- point
    norm(diff)
end

@everywhere function add_point(cluster::Cluster, point::Point)
    cluster.point_sum += point
    cluster.size += 1
end

@everywhere function calc_mean(cluster::Cluster)
    cluster.mean = cluster.point_sum ./ cluster.size
end

@everywhere function +(cluster1::Cluster, cluster2::Cluster)
    Cluster(cluster1.size + cluster2.size, cluster1.point_sum .+ cluster2.point_sum, cluster1.mean)
end

@everywhere function compute_local_labels_and_clusters(k::Int, points::Array{Point}, clusters::Array{Cluster})
    new_clusters = [Cluster() for _ in 1:k]
    labels = zeros(length(points))
    for point_idx in 1:length(points)
        point = points[point_idx]
        min_value = Inf
        min_index = 0
        for (idx, cluster) in enumerate(clusters)
            dist = distance(cluster, point)
            if dist < min_value
                min_value = dist
                min_index = idx
            end
        end
        labels[point_idx] = min_index
        add_point(new_clusters[min_index], point)
    end

    (new_clusters, labels)
end

@everywhere function work(k::Int, max_iter::Int, dpoints::DArray{Point}, initial_clusters::Array{Cluster}, clusters_channels::Array{RemoteChannel{ClustersCh}})
    points = localpart(dpoints)
    clusters = copy(initial_clusters)
    rank = myid()
    labels::Array{Int} = []
    for iter in 1:max_iter
        (my_clusters, labels) = compute_local_labels_and_clusters(k, points, clusters)

        for worker in workers()
            if worker != rank
                put!(clusters_channels[worker - 1], my_clusters)
            end
        end

        all_clusters::Array{Array{Cluster}} = [my_clusters]
        
        for _ in 1:nworkers() - 1
            push!(all_clusters, take!(clusters_channels[rank - 1]))
        end

        new_clusters = sum(all_clusters)
        [calc_mean(cluster) for cluster in new_clusters]
        clusters = new_clusters
    end
    labels
end

function main(args)
    # Configuration
    filename = args[1]
    out_filename = "out.txt"
    k = parse(Int64, args[2])
    max_iter = parse(Int64, args[3])

    # Input points
    points_read = readdlm(filename)::Array{Float64, 2}
    num_points = size(points_read, 1)
    points::Array{Point} = [Point(points_read[i, :]) for i in 1:num_points]
    arr = distribute(points)
    clusters::Array{Cluster} = []
    world_size = nworkers()

    # Algorithm
    start_time = time_ns()
    clusters_channels = [RemoteChannel(()->ClustersCh(world_size)) for _ in 1:world_size]

    @inbounds for i in 1:k
        idx = trunc(Int64, floor(rand() * num_points + 1))
        push!(clusters, Cluster(points[idx]))
    end

    tasks = []
    @sync for worker in workers()
        zb_rank = worker - 2
        task = @spawnat worker work(k, max_iter, arr, clusters, clusters_channels)
        push!(tasks, task)
    end

    labels::Array{Int64} = []
    for task in tasks
        append!(labels, fetch(task))
    end

    elapsed = (time_ns() - start_time) / 1.0e9
    writedlm("out.txt", labels)
    elapsed
end

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
end
println(minimum(times))

