using Base.Threads
using DelimitedFiles
using Random

const OVERSAMPLING_FACTOR = 64
const BINS_TYPE = Array{Array{Array{Int,1},1},1}

### QuickSort
function partition(arr::Array{Int}, left::Int, right::Int)
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

function myqsort(arr::Array{Int}, left::Int, right::Int)
    if left >= right
        return
    end
    pivot_idx = partition(arr, left, right)
    myqsort(arr, left, pivot_idx - 1)
    myqsort(arr, pivot_idx, right)
end

function myqsort(arr::Array{Int})
    @inbounds myqsort(arr, 1, length(arr))
end
### 

function get_sample_keys(arr::Array{Int}, m::Int)
    # sampled_keys = rand(arr, m * OVERSAMPLING_FACTOR)
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

function binary_search(arr::Array{Int}, key::Int)
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

function map_keys_to_bins(arr::AbstractArray{Int,1}, sample_keys::Array{Int}, m::Int)
    nkeys = length(arr)
    index::Array{Int} = zeros(nkeys)
    for i in 1:nkeys
        index[i] = binary_search(sample_keys, arr[i])
    end
    index
end

function compute_bin_array(thread_id::Int, blocksize::Int, arr::Array{Int}, sample_keys::Array{Int}, bins::BINS_TYPE, m::Int)
    nkeys = length(arr)
    start = (thread_id - 1) * blocksize
    end_ = min(start + blocksize, nkeys)
    start += 1
    subarray_size = end_ - start + 1

    index = map_keys_to_bins(@view(arr[start:end_]), sample_keys, m)
    for i in 1:subarray_size
        bin_idx = index[i]
        push!(bins[thread_id][bin_idx], arr[start + i - 1])
    end
end

function bin(arr::Array{Int}, m::Int)::BINS_TYPE
    bins::BINS_TYPE = []
    for _ in 1:m
        row = [[] for _ in 1:m]
        push!(bins, row)
    end
    nkeys = length(arr)
    sample_keys = get_sample_keys(arr, m)
    blocksize = div(nkeys + m - 1, m)
    @threads for i in 1:m
        compute_bin_array(threadid(), blocksize, arr, sample_keys, bins, m)
    end
    bins
end

function compute_sorted_subarray(thread_id::Int, subarray_sizes::Array{Int}, subarray_sizes_cumsum::Array{Int}, sorted_subarrays::Array{Int}, bins::BINS_TYPE)
    sorted_subarray::Array{Int,1} = zeros(subarray_sizes[thread_id])
    offset = 1
    for row in bins
        bin = row[thread_id]
        nbinkeys = length(bin)
        sorted_subarray[offset:offset+nbinkeys-1] .= bin
        offset += nbinkeys
    end
    myqsort(sorted_subarray)
    end_ = subarray_sizes_cumsum[thread_id]
    start = end_ - subarray_sizes[thread_id] + 1
    sorted_subarrays[start:end_] .= sorted_subarray
end

function subsort(bins::BINS_TYPE, nkeys::Int, m::Int)
    tally::Array{Int} = zeros(m, m)
    for i in 1:m
        for j in 1:m
            tally[i, j] = length(bins[i][j])
        end
    end

    colsum = vec(sum(tally, dims=1))
    colcumsum = cumsum(colsum)

    sorted_subarrays::Array{Int} = zeros(nkeys)
    @threads for i in 1:m
        compute_sorted_subarray(threadid(), colsum, colcumsum, sorted_subarrays, bins)
    end
    sorted_subarrays
end

function main(args)
    arr = convert(Array{Int, 1}, readdlm(args[1])[:, 1])
    nbins = nthreads()
    
    start = time_ns()
    
    @inbounds bins = bin(arr, nbins)
    @inbounds sorted_array = subsort(bins, length(arr), nbins)
    
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
end
println(minimum(times))
