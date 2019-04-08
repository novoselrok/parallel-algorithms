using Distributed
using DelimitedFiles
@everywhere using DistributedArrays

const OVERSAMPLING_FACTOR = 128
const BINS_TYPE = Array{Array{Int,1},1}
@everywhere const IntArrayCh = Channel{Array{Int64,1}}

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

function get_sample_keys(arr::Array{Int}, m::Int)
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

@everywhere function compute_bin_array(arr::Array{Int}, sample_keys::Array{Int}, m::Int)
    bins::Array{Array{Int}} = [[] for _ in 1:m]

    index = map_keys_to_bins(arr, sample_keys, m)
    for (key_idx, bin_idx) in enumerate(index)
        push!(bins[bin_idx], arr[key_idx])
    end
    bins
end

@everywhere function bin_subsort(darr::DArray{Int}, sample_keys::Array{Int,1}, m::Int, channels::Array{RemoteChannel{IntArrayCh}})
    arr = localpart(darr)
    ob_rank = myid() - 1
    bins = compute_bin_array(arr, sample_keys, m)
    
    for (idx, channel) in enumerate(channels)
        if idx == ob_rank
            continue
        end
        put!(channel, bins[idx])
    end
    subarray::Array{Int,1} = bins[ob_rank]
    received_arrays = m - 1
    while received_arrays > 0
        received_array = take!(channels[ob_rank])
        subarray = vcat(subarray, received_array)
        received_arrays -= 1
    end
    myqsort(subarray)
    subarray
end

const REPEAT = 100
const MY_RAND_MAX = ((1 << 31) - 1)

function get_random_number(seed::Int)
    (seed * 1103515245 + 12345) & MY_RAND_MAX
end

function init_random_array(n::Int, initial_seed::Int)
    arr::Array{Int} = zeros(n)
    random_num = get_random_number(initial_seed)
    for i in 1:n
        arr[i] = random_num
        random_num = get_random_number(random_num)
    end
    arr
end

function main(args)
    times = []
    n = parse(Int64, args[1])

    for iter in 1:REPEAT
        randomarr = init_random_array(n, iter)
        arr = distribute(randomarr)
        nbins = nworkers()
        
        start_time = time_ns()
        sample_keys = get_sample_keys(randomarr, nbins)
        channels = [RemoteChannel(()->IntArrayCh(nbins)) for _ in 1:nworkers()]

        tasks = []
        @sync for worker in workers()
            zb_rank = worker - 2

            blocksize = div(n + nbins - 1, nbins)
            start = zb_rank * blocksize
            end_ = min(start + blocksize, n)
            start += 1
            subarray_size = end_ - start + 1

            task = @spawnat worker bin_subsort(arr, sample_keys, nbins, channels)
            push!(tasks, task)
        end

        sorted_array::Array{Int64} = []
        for task in tasks
            append!(sorted_array, fetch(task))
        end

        elapsed = (time_ns() - start_time) / 1.0e9

        for i in 1:length(sorted_array) - 1
            if sorted_array[i] > sorted_array[i + 1]
                println("Array not sorted.")
                exit(1)
            end
        end
        push!(times, elapsed)
        close(arr)
    end
    sum(times) / REPEAT
end

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
    sleep(1)
end
println(minimum(times))
