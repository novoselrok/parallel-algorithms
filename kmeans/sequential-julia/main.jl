using DelimitedFiles
using StaticArrays
using LinearAlgebra

const POINT_SIZE = 128
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
    diff = cluster.mean .- point
    norm(diff)
end

function add_point(cluster::Cluster, point::Point)
    cluster.point_sum += point
    cluster.size += 1
end

function calc_mean(cluster::Cluster)
    cluster.mean = cluster.point_sum ./ cluster.size
end

function compute_label_for_point(index::Int, points::Array{Point}, clusters::Array{Cluster}, labels::Array{Int}, new_clusters::Array{Cluster})
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
    add_point(new_clusters[min_index], point)
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
    points::Array{Point, 1} = [Point(points_read[i, :]) for i in 1:num_points]
    labels = convert(Array{Int64, 1}, zeros(num_points))
    clusters::Array{Cluster} = []

    # Algorithm
    start = time_ns()

    @inbounds for i in 1:k
        idx = trunc(Int64, floor(rand() * num_points + 1))
        push!(clusters, Cluster(points[idx]))
    end

    @inbounds for iter in 1:max_iter
        new_clusters = [Cluster() for _ in 1:k]
        for i in 1:num_points
            compute_label_for_point(i, points, clusters, labels, new_clusters)
        end
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
