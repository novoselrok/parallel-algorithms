using DelimitedFiles

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
        arr = init_random_array(n, iter)
        start = time_ns()
        myqsort(arr)
        elapsed = ((time_ns() - start) / 1.0e9)
    
        for i in 1:length(arr) - 1
            if arr[i] > arr[i + 1]
                println("Array not sorted.")
                exit(1)
            end
        end
        push!(times, elapsed)
    end
    sum(times) / REPEAT
end

nprecompilesteps = haskey(ENV, "JL_NRETRIES") ? parse(Int, ENV["JL_NRETRIES"]) : 0
times = []
for i in 1:nprecompilesteps
    push!(times, main(ARGS))
end
println(minimum(times))
