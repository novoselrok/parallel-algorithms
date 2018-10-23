using Distributed
n_tasks = 8
addprocs(n_tasks)

@everywhere using DataStructures
@everywhere using DelimitedFiles
@everywhere using Profile

graph = convert(Array{Int, 2}, readdlm("../data/mat"))
n_cities = size(graph, 1)
@everywhere home = 1

@everywhere mutable struct Tour
    cities::Array{Int, 1}
    cities_set::Set{Int}
    cost::Int
end

function initial_tour(city::Int)
    cities = Array{Int, 1}()
    cities_set = Set{Int}()
    pushfirst!(cities, city)
    push!(cities_set, city)
    Tour(cities, cities_set, 0)
end

@everywhere function add_city(t::Tour, city::Int)
    @inbounds t.cost += graph[t.cities[1], city]
    pushfirst!(t.cities, city)
    push!(t.cities_set, city)
end

@everywhere function total_cost(t::Tour)
    @inbounds t.cost + graph[t.cities[1], home]
end

@everywhere Base.copy(t::Tour) = Tour(copy(t.cities), copy(t.cities_set), t.cost)

best_tours = RemoteChannel(() -> Channel{Tour}(1))
best_tour = initial_tour(home)
best_tour.cost = typemax(Int)

function passobj(src::Int, target::AbstractVector{Int}, nm::Symbol;
    from_mod=Main, to_mod=Main)
    r = RemoteChannel(src)
    @spawnat(src, put!(r, getfield(from_mod, nm)))
    @sync for to in target
        @spawnat(to, Core.eval(to_mod, Expr(:(=), nm, fetch(r))))
    end
    nothing
end

macro passobj(src::Int, target, val, from_mod=:Main, tomod=:Main)
    quote
        passobj($(esc(src)), $(esc(target)), $(QuoteNode(val)); from_mod=$from_mod, to_mod=$tomod)
    end
end


function passobj(src::Int, target::Int, nm::Symbol; from_mod=Main, to_mod=Main)
    passobj(src, [target], nm; from_mod=from_mod, to_mod=to_mod)
end


function passobj(src::Int, target, nms::Vector{Symbol}; from_mod=Main, to_mod=Main)
    for nm in nms
        passobj(src, target, nm; from_mod=from_mod, to_mod=to_mod)
    end
end

@everywhere function is_best_tour(t::Tour)
    total_cost(t) < total_cost(best_tour)
end

function update_best_tour()
    println("updating")
    while true
        t = take!(best_tours)
        if is_best_tour(t)
            add_city(t, home)
            global best_tour = copy(t)
            @passobj 1 workers() best_tour
            println(best_tour)
        end
    end
    println("done updating")
end

@everywhere function par_tree_search(s::Stack{Tour}, graph::Array{Int, 2}, n_cities::Int, best_tours::RemoteChannel)
    while !isempty(s)
        t = pop!(s)
        # push!(time_elapsed, @elapsed (is_best = remotecall_fetch((t) -> is_best_tour(t), 1, t)))
        is_best = is_best_tour(t)
        if length(t.cities_set) == n_cities && is_best
            put!(best_tours, t)
        else
            for city in n_cities:-1:2
                if in(city, t.cities_set) || @inbounds t.cost + graph[t.cities[1], city] >= best_tour.cost
                    continue
                end
                new_tour = copy(t)
                add_city(new_tour, city)
                push!(s, new_tour)
            end
        end
    end
end

function main()
    @passobj 1 workers() best_tour

    q = Array{Tour, 1}()
    push!(q, initial_tour(home))

    while length(q) < n_tasks
        t = popfirst!(q)
        for city in 2:n_cities
            if in(city, t.cities_set)
                continue
            end

            new_tour = copy(t)
            add_city(new_tour, city)
            push!(q, new_tour)
        end
    end

    stack_per_task = [Stack{Tour}() for _ in 1:n_tasks]
    tours_per_task = convert(Int, ceil((length(q) / n_tasks)))
    # idx = 1
    for i in 0:n_tasks-1
        fst = i * tours_per_task + 1
        lst = i * tours_per_task + tours_per_task + 1

        if lst > length(q)
            sub = q[fst:end]
        else
            sub = q[fst:lst]
        end

        for t in sub
            push!(stack_per_task[i + 1], t)
        end
    end

    updater = @async update_best_tour()

    @time @sync for s in stack_per_task
        @spawn par_tree_search(s, graph, n_cities, best_tours)
    end

    # fetch(updater)
    println("done")
end

main()
