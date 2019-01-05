using DelimitedFiles
using StaticArrays

const POINT_SIZE = 10
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

function compute_label_for_point(index::Int, points::Array{Point}, clusters::Array{Cluster}, labels::Array{Int}, new_clusters::Array{Cluster})
    point = points[index]
    distances = [distance(cluster, point) for cluster in clusters]
    _, min_index = findmin(distances)
    labels[index] = min_index
    add_point(new_clusters[min_index], point)
end

function kmeanspp_init(points::Array{Point}, k::Int)
    clusters::Array{Cluster} = []
    push!(clusters, Cluster(points[1]))

    for _ in 2:k
        distances = [minimum([distance(cluster, point) for cluster in clusters]) for point in points]
        distances_sum = sum(distances)
        cumulative_probs = cumsum(distances ./ distances_sum)
        random = rand()
        new_idx = 0
        for (i, prob) in enumerate(cumulative_probs)
            if random < prob
                new_idx = i
                break
            end
        end
        push!(clusters, Cluster(points[new_idx]))
    end
    clusters
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
    # Algorithm
    start = time_ns()
    clusters = kmeanspp_init(points, k)

    @inbounds for iter in 1:max_iter
        new_clusters = [Cluster() for _ in 1:k]
        for i in 1:num_points
            compute_label_for_point(i, points, clusters, labels, new_clusters)
        end
        [calc_mean(cluster) for cluster in new_clusters]
        clusters = new_clusters
    end
    println((time_ns() - start) / 1.0e9)

    writedlm("out.txt", labels)
end

main(ARGS)
