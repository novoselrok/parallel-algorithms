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

function main(args, verbose::Bool)
    # arr = [10, 18, 16, 14, 0, 17, 11, 2, 3, 9, 5, 7, 4, 19, 6, 15, 8, 1, 13, 12]
    arr = convert(Array{Int, 1}, readdlm(args[1])[:, 1])

    start = time_ns()
    myqsort(arr)
    verbose && println((time_ns() - start) / 1.0e9)

    for i in 1:length(arr) - 1
        if arr[i] > arr[i + 1]
            println("Array not sorted.")
            exit(1)
        end
    end
end

nprecompilesteps = 3
verbose = false
for i in 1:nprecompilesteps
    main(ARGS, verbose)
end
verbose = true
main(ARGS, verbose)
