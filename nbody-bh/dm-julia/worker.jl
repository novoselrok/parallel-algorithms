include("common.jl")
include("body.jl")
include("cell.jl")

using Distributed
@everywhere using DelimitedFiles

@everywhere const BISECTION_MAX_ITER = 200
@everywhere const BISECTION_TOL = 1e-10
@everywhere const dt = 0.1

@everywhere const Bounds = Tuple{Vec3, Vec3}
@everywhere const Vec3Ch = Channel{Vec3}
@everywhere const FloatCh = Channel{Tuple{Float64, Float64}}
@everywhere const BodiesCh = Channel{Array{Body}}
@everywhere const CellTuplesCh = Channel{Array{CellTuple}}
@everywhere const BoundsCh = Channel{Bounds}
@everywhere const GroupPartnerInfo = Tuple{Int, Int, RemoteChannel{FloatCh}}
@everywhere const GroupCh = Channel{GroupPartnerInfo}
@everywhere const PartnerInfo = Tuple{Int, Bool}

@everywhere function get_workers_info()
    (myid() - 2, myid() - 1, nworkers())
end

@everywhere function get_universe_size(bodies::Array{Body})::Bounds
    universe_min = [Inf, Inf, Inf]
    universe_max = [-Inf, -Inf, -Inf]

    for body in bodies
        for c in 1:DIMS
            if body.position[c] < universe_min[c]
                universe_min[c] = body.position[c]
            end
            if body.position[c] > universe_max[c]
                universe_max[c] = body.position[c]
            end
        end
    end

    (Vec3(universe_min), Vec3(universe_max))
end

@everywhere function determine_group_partners(
    group::Int,
    group_channels::Array{RemoteChannel{GroupCh}},
    float_channels::Array{RemoteChannel{FloatCh}}
    )::Array{RemoteChannel{FloatCh}}

    (zb_rank, ob_rank, world_size) = get_workers_info()

    for (idx, channel) in enumerate(group_channels)
        if idx != ob_rank
            put!(channel, (group, ob_rank, float_channels[ob_rank]))
        end
    end

    group_partners::Array{GroupPartnerInfo} = []
    for _ in 1:world_size - 1
        info = take!(group_channels[ob_rank])
        if info[1] == group
            push!(group_partners, info)
        end
    end

    map(x -> x[3], sort(group_partners, by = x -> x[2]))
end

@everywhere function is_above_split(rank::Int, n_procs_left::Int)::Bool
    rel_bit = trunc(Int, log2(n_procs_left))
    is_above = (rank >> (rel_bit - 1)) & 1 # get last bit
    is_above == 1
end

@everywhere function get_partner_rank(rank::Int, n_procs_left::Int)
    rel_bit = trunc(Int, log2(n_procs_left))
    xor(rank, (1 << (rel_bit - 1))) # rank xor (2^(rel_bit-1))
end

@everywhere function frac_weight_below(
    bodies::Array{Body},
    split::Float64,
    coord::Int,
    my_channel::RemoteChannel{FloatCh},
    group_partners::Array{RemoteChannel{FloatCh}})

    work_below = 0.0
    work_above = 0.0

    for body in bodies
        if body.position[coord] > split
            work_above += body.work
        else
            work_below += body.work
        end
    end

    all_work_below = work_below
    all_work_above = work_above

    for partner in group_partners
        put!(partner, (work_below, work_above))
    end

    for i in 1:length(group_partners)
        (other_work_below, other_work_above) = take!(my_channel)
        all_work_below += other_work_below
        all_work_above += other_work_above
    end

    fraction = all_work_below / (all_work_below + all_work_above)
    fraction - 0.5
end

@everywhere function bisection(
    min::Float64,
    max::Float64,
    bodies::Array{Body},
    coord::Int,
    my_channel::RemoteChannel{FloatCh},
    group_partners::Array{RemoteChannel{FloatCh}})

    fmin = frac_weight_below(bodies, min, coord, my_channel, group_partners)
    mid = 0.0

    iter = 0
    while iter < BISECTION_MAX_ITER && abs((max - min) / 2) > BISECTION_TOL
        mid = (min + max) / 2
        fmid = frac_weight_below(bodies, mid, coord, my_channel, group_partners)
        if abs(fmid) < BISECTION_TOL
            break
        elseif fmin * fmid > 0
            min = mid
        else
            max = mid
        end

        iter += 1
    end

    mid
end

@everywhere function orb(
    bodies::Array{Body},
    universe_min::Vec3,
    universe_max::Vec3,
    bodies_channels::Array{RemoteChannel{BodiesCh}},
    group_channels::Array{RemoteChannel{GroupCh}},
    float_channels::Array{RemoteChannel{FloatCh}})

    (zb_rank, ob_rank, world_size) = get_workers_info()
    my_bodies::Array{Body} = []
    other_bodies::Array{Body} = []
    my_bounds::Array{Bounds} = []
    other_bounds::Array{Bounds} = []
    partners::Array{PartnerInfo} = []

    my_min = copy(universe_min)
    my_max = copy(universe_max)

    group = 0
    above_split = false

    n_splits = trunc(Int, log2(world_size))
    for i in 0:n_splits-1
        n_procs_left = trunc(Int, world_size / (2^i));

        group = group << 1
        if above_split
            group = group | 1
        end

        group_partners = determine_group_partners(group, group_channels, float_channels)

        coord = (i % DIMS) + 1

        split = bisection(my_min[coord], my_max[coord], bodies, coord, float_channels[ob_rank], group_partners)

        above_split = is_above_split(zb_rank, n_procs_left)

        other_min = copy(my_min)
        other_max = copy(my_max)

        if above_split
            my_min = setindex(my_min, split, coord)
            other_max = setindex(other_max, split, coord)
        else
            my_max = setindex(my_max, split, coord)
            other_min = setindex(other_min, split, coord)
        end

        push!(my_bounds, (my_min, my_max))
        push!(other_bounds, (other_min, other_max))

        my_bodies = []
        other_bodies = []

        for body in bodies
            if (body.position[coord] - split > 0) == above_split
                push!(my_bodies, body)
            else
                push!(other_bodies, body)
            end
        end

        # partner rank used for indexing into the body channels
        partner_rank = get_partner_rank(zb_rank, n_procs_left) + 1
        push!(partners, (partner_rank, above_split))
        put!(bodies_channels[partner_rank], other_bodies)
        recv_bodies = take!(bodies_channels[ob_rank])

        bodies = vcat(my_bodies, recv_bodies)
        my_min = copy(my_min)
        my_max = copy(my_max)
    end

    (bodies, my_bounds, other_bounds, partners)
end

@everywhere function work(
    input_bodies::Array{Body},
    bodies_channels::Array{RemoteChannel{BodiesCh}},
    bounds_channels::Array{RemoteChannel{BoundsCh}},
    group_channels::Array{RemoteChannel{GroupCh}},
    float_channels::Array{RemoteChannel{FloatCh}},
    cell_channels::Array{RemoteChannel{CellTuplesCh}},
    iterations::Int)
    # zero based rank, used for ORB calculation
    # one based rank, used for indexing
    (zb_rank, ob_rank, world_size) = get_workers_info()

    bodies::Array{Body} = copy(input_bodies)

    for iter in 1:iterations
        (universe_min, universe_max) = get_universe_size(bodies)

        for body in bodies
            reset_force(body)
        end

        # Calculate universe bounds
        for (idx, channel) in enumerate(bounds_channels)
            if idx != ob_rank
                put!(channel, (universe_min, universe_max))
            end
        end

        for _ in 1:world_size - 1
            (other_universe_min, other_universe_max) = take!(bounds_channels[ob_rank])
            universe_min = min.(universe_min, other_universe_min)
            universe_max = max.(universe_max, other_universe_max)
        end

        # ORB
        (bodies, my_bounds, other_bounds, partners) = orb(bodies, universe_min, universe_max, bodies_channels, group_channels, float_channels)

        root = Cell()
        root.min_bounds = universe_min
        root.max_bounds = universe_max
        root.cell_opened = true
        
        for (min_bounds, max_bounds) in my_bounds
            insert_empty_cell(root, min_bounds, max_bounds)
        end

        for body in bodies
            insert_body(root, body)
        end

        for i in 1:length(my_bounds)
            (my_min, my_max) = my_bounds[i]
            (other_min, other_max) = other_bounds[i]

            cells_to_send = get_cells_to_send(root, Cell(), other_min, other_max, i - 1, 0, Array{Cell}([]))
            packed = pack_cells(cells_to_send)

            (partner, above_split) = partners[i]
            put!(cell_channels[partner], packed)
            recv_cells = take!(cell_channels[ob_rank])
            root_cells = reconstruct_received_cells(recv_cells)

            for root_cell in root_cells
                insert_cell(root, root_cell)
            end
        end

        for body in bodies
            elapsed_force_time = @elapsed compute_force(root, body)
            body.work = elapsed_force_time
        end

        for body in bodies
            body.position = body.position .+ (dt .* body.velocity)
            body.velocity = body.velocity .+ (dt ./ body.mass .* body.force)
        end
    end

    n = length(bodies)
    ids::Array{Float64, 2} = zeros(n, 1)
    final_positions::Array{Float64, 2} = zeros(n, 3)
    final_velocities::Array{Float64, 2} = zeros(n, 3)
    for i in 1:n
        body = bodies[i]
        ids[i, 1] = body.id
        final_positions[i, :] = body.position
        final_velocities[i, :] = body.velocity
    end
    writedlm(string("out", zb_rank, ".txt"), hcat(ids, final_positions, final_velocities))

end
