using Distributed
using DelimitedFiles
@everywhere using SharedArrays

const OVERSAMPLING_FACTOR = 128
const BINS_TYPE = Array{Array{Int,1},1}

### QuickSort
@everywhere function partition(arr::Array{Int}, left::Int, right::Int)
    pivot = arr[right]
    while left < right
        while arr[left] < pivot 
            left += 1
        end
        while arr[right] > pivot
            right -= 1
        end
        if left <= right
            temp = arr[left]
            arr[left] = arr[right]
            arr[right] = temp
        end
    end
    left
end

@everywhere function myqsort(arr::Array{Int}, left::Int, right::Int)
    if left >= right
        return
    end
    pivot_idx = partition(arr, left, right)
    myqsort(arr, left, pivot_idx - 1)
    myqsort(arr, pivot_idx, right)
end

@everywhere function myqsort(arr::Array{Int})
    @inbounds myqsort(arr, 1, length(arr))
end
### 

function get_sample_keys(arr::SharedArray{Int}, m::Int)
    sampled_keys = zeros(m * OVERSAMPLING_FACTOR)
    for i in 1:length(sampled_keys)
        sampled_keys[i] = arr[i]
    end
    sample_keys::Array{Int} = []
    for i in 1:m-1
        push!(sample_keys, sampled_keys[i * OVERSAMPLING_FACTOR])
    end
    myqsort(sample_keys)
    sample_keys
end

@everywhere function binary_search(arr::Array{Int}, key::Int)
    left = 0
    right = length(arr)
    while left < right
        middle = div(left + right, 2)
        if arr[middle + 1] >= key
            right = middle
        else
            left = middle + 1
        end
    end
    left + 1
end

@everywhere function map_keys_to_bins(arr::AbstractArray{Int,1}, sample_keys::Array{Int}, m::Int)
    nkeys = length(arr)
    index::Array{Int} = zeros(nkeys)
    for i in 1:nkeys
        index[i] = binary_search(sample_keys, arr[i])
    end
    index
end

@everywhere function compute_bin_array(arr::SharedArray{Int}, sample_keys::Array{Int}, m::Int)
    local_indices = localindices(arr)
    subarray_size = length(local_indices)
    bins = [[] for _ in 1:m]

    index = map_keys_to_bins(@view(arr[local_indices]), sample_keys, m)
    for (bin_idx, key_idx) in zip(index, local_indices)
        push!(bins[bin_idx], arr[key_idx])
    end
    bins
end

@everywhere function bin_subsort(arr::SharedArray{Int}, sample_keys::Array{Int,1}, m::Int, channels::Array{RemoteChannel{Channel{Array{Int64,1}}},1})
    workerid = myid() - 1
    bins = compute_bin_array(arr, sample_keys, m)
    
    for (idx, channel) in enumerate(channels)
        if idx == workerid
            continue
        end
        put!(channel, bins[idx])
    end
    subarray::Array{Int,1} = bins[workerid]
    received_arrays = m - 1
    while received_arrays > 0
        received_array = take!(channels[workerid])
        subarray = vcat(subarray, received_array)
        received_arrays -= 1
    end
    myqsort(subarray)
    subarray
end

function main(args)
    arr = convert(Array{Int, 1}, readdlm(args[1])[:, 1])
    arr = SharedArray{Int}(arr)
    nbins = nworkers()
    
    start = time_ns()
    sample_keys = get_sample_keys(arr, nbins)
    channels = [RemoteChannel(()->Channel{Array{Int,1}}(nbins)) for _ in 1:nworkers()]

    sorted_array = @sync @distributed (vcat) for i in 1:nbins
        bin_subsort(arr, sample_keys, nbins, channels)
    end

    elapsed = (time_ns() - start) / 1.0e9

    for i in 1:length(sorted_array) - 1
        if sorted_array[i] > sorted_array[i + 1]
            println("Array not sorted.")
            exit(1)
        end
    end
    elapsed
end

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
    println(times)
    sleep(1)
end
println(minimum(times))
