using Base.Threads
function compute_pi(n::Int)
    within_circle = 0
    for _ in 1:n
        x = rand() * 2 - 1
        y = rand() * 2 - 1
        r2 = x * x + y * y
        if r2 < 1.0
            within_circle += 1
        end
    end
    within_circle / n * 4.0
end
const N = 10000000
results = zeros(nthreads())
@threads for _ in 1:nthreads()
    results[threadid()] = compute_pi(ceil(Int, N / nthreads()))
end
println(sum(results) / nthreads())

