using DelimitedFiles

function main(args)
    # arr = [10, 18, 16, 14, 0, 17, 11, 2, 3, 9, 5, 7, 4, 19, 6, 15, 8, 1, 13, 12]
    arr = convert(Array{Int, 1}, readdlm(args[1])[:, 1])

    start = time_ns()
    sort!(arr)
    println((time_ns() - start) / 1.0e9)

    for i in 1:length(arr) - 1
        if arr[i] > arr[i + 1]
            println("Array not sorted.")
            exit(1)
        end
    end

end

main(ARGS)
