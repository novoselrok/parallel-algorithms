using DelimitedFiles
using StaticArrays
using Base.Threads

import Base.+

const POINT_SIZE = 100
# Type alias for Point type
const Point = SVector{POINT_SIZE, Float64}

# Datatypes
mutable struct Cluster
   size::Int64
   point_sum::Point
   mean::Point
end

# Constructor
Cluster() = Cluster(0, Point(zeros(POINT_SIZE)), Point(zeros(POINT_SIZE)))
Cluster(point::Point) = Cluster(0, Point(zeros(POINT_SIZE)), point)

function distance(cluster::Cluster, point::Point)
    sqrt(sum((cluster.mean - point) .^ 2))
end

function add_point(cluster::Cluster, point::Point)
    cluster.point_sum += point
    cluster.size += 1
end

function calc_mean(cluster::Cluster)
    cluster.mean = cluster.point_sum ./ cluster.size
end

function reset_cluster(cluster::Cluster)
    cluster.size = 0
    cluster.point_sum = Point(zeros(POINT_SIZE))
    cluster.mean = Point(zeros(POINT_SIZE))
end

function +(cluster1::Cluster, cluster2::Cluster)
    Cluster(cluster1.size + cluster2.size, cluster1.point_sum .+ cluster2.point_sum, cluster1.mean)
end

function compute_label_for_point(index::Int, points::Array{Point}, clusters::Array{Cluster}, labels::Array{Int}, clusters_per_thread::Array{Array{Cluster,1},1})
    thread_id = threadid()
    point = points[index]
    min_value = Inf
    min_index = 0
    for (idx, cluster) in enumerate(clusters)
        dist = distance(cluster, point)
        if dist < min_value
            min_value = dist
            min_index = idx
        end
    end
    labels[index] = min_index
    add_point(clusters_per_thread[thread_id][min_index], point)
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
    labels = convert(Array{Int64, 1}, zeros(num_points))
    clusters::Array{Cluster} = []

    # Algorithm
    start = time_ns()

    @inbounds for i in 1:k
        idx = convert(Int64, floor(rand() * num_points + 1))
        push!(clusters, Cluster(points[idx]))
    end

    for iter in 1:max_iter
        clusters_per_thread = [[Cluster() for _ in 1:k] for _ in 1:nthreads()]

        @threads for i in 1:num_points
            compute_label_for_point(i, points, clusters, labels, clusters_per_thread)
        end
        new_clusters = sum(clusters_per_thread)

        [calc_mean(cluster) for cluster in new_clusters]
        clusters = new_clusters
    end
    
    println((time_ns() - start) / 1.0e9)
    writedlm("out.txt", labels)
end

main(ARGS)
