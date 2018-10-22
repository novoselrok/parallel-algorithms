use List;
use Types;
use Time;

config const filename = "input.txt";
config const n_tasks = 4;
var home = 0;
var tasks_running: atomic int;
var n_cities: int;
var graph_domain: domain(2);
var graph: [graph_domain] int;

record Tour {
    var cities: domain(int);
    var cities_ordered: list(int) = new list(int);
    var cost: int = 0;
    var cost_to_home: int = 0;
    var last_city: int = 0;

    proc init() {
        cities_ordered = new list(int);
    }

    proc init(city: int) {
        cities = {city};
        cities_ordered = new list(int);
        cities_ordered.push_front(city);
    }

    proc add(city: int) {
        // var last_city = cities.pop_front();
        cost += graph[last_city, city];
        // cities.push_front(last_city);
        cities_ordered.push_front(city);
        cities += city;
        last_city = city;
        cost_to_home = graph[last_city, home];
    }

    proc visited(city: int) {
        return cities.member(city);
    }

    proc set_cost(new_cost: int) {
        cost = new_cost;
    }

    proc total_cost() {
        return cost + cost_to_home;
    }

    proc length() {
        return cities.size;
    }
}

var best_tour: Tour;
best_tour.set_cost(max(int));
var best_tour$: sync bool;

proc is_best_tour(tour: Tour) {
    return tour.total_cost() < best_tour.total_cost();
}

proc update_best_tour(ref tour: Tour) {
    if (is_best_tour(tour)) {
        tour.add(home);
        best_tour = tour;
        writeln(" update tour ", best_tour);
    }
}

proc par_tree_search(rank: int, ref stack: list(Tour)) {
    tasks_running.add(1);

    while (stack.length > 0) {
        // if should_split?
        // writeln(rank, " loop start ", stack.length);
        var tour: Tour = stack.pop_front();
        if (tour.length() == n_cities && is_best_tour(tour)) {
            // writeln(rank, " waiting");
            best_tour$ = true; // fill/write best_tour$ i.e. lock
            update_best_tour(tour);
            var val = best_tour$; // empty/read best_tour$ i.e. unlock
        } else {
            // writeln(rank, " adding cities");
            for city in 1..n_cities - 1 by -1 {
                if (tour.visited(city) || tour.cost + graph[tour.last_city, city] >= best_tour.cost) {
                    continue;
                }
                var new_tour = tour; // copy
                new_tour.add(city);
                stack.push_front(new_tour);
            }
        }
        // writeln(rank, " loop done ", stack.length);
    }
    tasks_running.sub(1);
}


proc main() {
    tasks_running.write(0);

    var f = open(filename, iomode.r);
    var reader = f.reader();
    n_cities = reader.read(int);
    graph_domain = {0..#n_cities, 0..#n_cities};
    reader.read(graph);
    writeln(graph);
    f.close();
    reader.close();
    writeln("Done reading...");

    var initial_tours_q = new list(Tour);
    // 0 indicates the starting city
    initial_tours_q.append(new Tour(home));

    while (initial_tours_q.length < n_tasks) {
        var tour = initial_tours_q.pop_front();
        for city in 1..n_cities - 1 {
            if (tour.visited(city)) {
                continue;
            }
            var new_tour = tour; // copy
            new_tour.add(city);
            initial_tours_q.append(new_tour);
        }
    }

    var stack_per_task: [0..#n_tasks] list(Tour);

    var counter = 0;
    for tour in initial_tours_q {
        stack_per_task[counter % n_tasks].append(tour);
        counter += 1;
    }

    sync {
        coforall i in 0..#n_tasks {
            // writeln(i, " ", stack_per_task[i]);
            par_tree_search(i, stack_per_task[i]);
        }
    }
    
    writeln(best_tour);
}
