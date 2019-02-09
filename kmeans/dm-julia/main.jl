using Distributed
using DelimitedFiles
@everywhere using StaticArrays
@everywhere using SharedArrays

import Base.+

@everywhere const POINT_SIZE = 100
# Type alias for Point type
@everywhere const Point = SVector{POINT_SIZE, Float64}

# Datatypes
@everywhere mutable struct Cluster
   size::Int64
   point_sum::Point
   mean::Point
end

# Constructor
@everywhere Cluster() = Cluster(0, Point(zeros(POINT_SIZE)), Point(zeros(POINT_SIZE)))
@everywhere Cluster(point::Point) = Cluster(0, Point(zeros(POINT_SIZE)), point)

@everywhere function distance(cluster::Cluster, point::Point)
    sqrt(sum((cluster.mean .- point) .^ 2))
end

@everywhere function add_point(cluster::Cluster, point::Point)
    cluster.point_sum += point
    cluster.size += 1
end

@everywhere function calc_mean(cluster::Cluster)
    cluster.mean = cluster.point_sum ./ cluster.size
end

function +(cluster1::Cluster, cluster2::Cluster)
    Cluster(cluster1.size + cluster2.size, cluster1.point_sum .+ cluster2.point_sum, cluster1.mean)
end

@everywhere function compute_local_labels_and_clusters(k::Int, points::SharedArray{Point}, labels::SharedArray{Int64}, clusters::Array{Cluster})
    local_indices = localindices(points)
    new_clusters = [Cluster() for _ in 1:k]

    for point_idx in local_indices
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

    new_clusters
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
    points::SharedArray{Point} = SharedArray{Point}([Point(points_read[i, :]) for i in 1:num_points])
    labels::SharedArray{Int64} = SharedArray{Int64}(zeros(num_points))
    clusters::Array{Cluster} = []

    # Algorithm
    start = time_ns()

    @inbounds for i in 1:k
        idx = convert(Int64, floor(rand() * num_points + 1))
        push!(clusters, Cluster(points[idx]))
    end

    @inbounds for iter in 1:max_iter
        tasks = @sync [@spawnat worker compute_local_labels_and_clusters(k, points, labels, clusters) for worker in workers()]
        computed_new_clusters = [fetch(task) for task in tasks]

        new_clusters = sum(computed_new_clusters)
        [calc_mean(cluster) for cluster in new_clusters]
        clusters = new_clusters
    end

    elapsed = (time_ns() - start) / 1.0e9
    writedlm("out.txt", labels)
    elapsed
end

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
end
println(minimum(times))
