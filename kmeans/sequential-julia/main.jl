using DelimitedFiles
using StaticArrays

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

function main(args)
    # Configuration
    filename = args[1]
    out_filename = "out.txt"
    k = parse(Int64, args[2])
    max_iter = parse(Int64, args[3])

    # Input points
    points = readdlm(filename)::Array{Float64, 2}
    num_points = size(points, 1)
    points = [Point(points[i, :]) for i in 1:num_points]
    labels = convert(Array{Int64, 1}, zeros(num_points))
    
    # Algorithm
    start = time_ns()
    clusters = [Cluster() for i in 1:k]

    @inbounds for i in 1:k
        idx = convert(Int64, floor(rand() * num_points + 1))
        clusters[i].mean = copy(points[idx])
    end

    @inbounds for iter in 1:max_iter
        new_clusters = [Cluster() for _ in 1:k]

        for i in 1:num_points
            distances = [distance(cluster, points[i]) for cluster in clusters]
            _, min_index = findmin(distances)
            labels[i] = min_index
            add_point(new_clusters[min_index], points[i])
        end
        [calc_mean(cluster) for cluster in new_clusters]
        clusters = new_clusters
    end
    println((time_ns() - start) / 1.0e9)

    writedlm("out.txt", labels)
end

main(ARGS)
