include("common.jl")
include("body.jl")

@everywhere import Base: @propagate_inbounds

@everywhere const N_CELL_CHILDREN = 8

@everywhere mutable struct Cell
    cell_opened::Bool # If true it means it has children
    children::Array{Cell, 1}
    body_present::Bool
    body::Body
    mass::Float64
    cm::Vec3
    min_bounds::Vec3
    max_bounds::Vec3
    parent_idx::Int
    array_idx::Int
end

@everywhere Cell() = Cell(
    false,
    Array{Cell, 1}(undef, N_CELL_CHILDREN),
    false,
    Body(),
    0.0, 
    Vec3(zeros(DIMS)),
    Vec3(zeros(DIMS)),
    Vec3(zeros(DIMS)),
    -1,
    -1
)

@everywhere const CellTuple = Tuple{Int, Vec3, Vec3, Vec3, Float64}

@everywhere function is_external(cell::Cell)
    # !cell.cell_opened
    !isassigned(cell.children, 1)
end

@everywhere function cell_contains_position(cell::Cell, position::Vec3)
    if position[X] < cell.min_bounds[X] || position[X] > cell.max_bounds[X]
        return false
    end
    if position[Y] < cell.min_bounds[Y] || position[Y] > cell.max_bounds[Y]
        return false
    end
    if position[Z] < cell.min_bounds[Z] || position[Z] > cell.max_bounds[Z]
        return false
    end
    true
end

@everywhere function pack_cells(cells::Array{Cell})::Array{CellTuple}
    map((cell) -> (cell.parent_idx, cell.cm, cell.min_bounds, cell.max_bounds, cell.mass), cells)
end

@everywhere function reconstruct_received_cells(recv_cells::Array{CellTuple})
    all_cells::Array{Cell} = []
    root_cells::Array{Cell} = []

    for (parent_idx, cm, min_bounds, max_bounds, mass) in recv_cells
        cell = Cell()
        cell.cm = copy(cm)
        cell.min_bounds = copy(min_bounds)
        cell.max_bounds = copy(max_bounds)
        cell.mass = mass

        if parent_idx == -1 
            push!(root_cells, cell)
        else
            parent = all_cells[parent_idx]
            for i in 1:N_CELL_CHILDREN
                if !isassigned(parent.children, i)
                    parent.children[i] = cell
                    break
                end
            end
        end
        push!(all_cells, cell)
    end
    root_cells
end

@everywhere function get_cells_to_send(
    cell::Cell,
    parent::Cell,
    min_bounds::Vec3,
    max_bounds::Vec3,
    min_depth::Int,
    depth::Int,
    cells_to_send::Array{Cell})

    if depth > min_depth
        cell.parent_idx = depth - min_depth == 1 ? -1 : parent.array_idx
        cell.array_idx = length(cells_to_send) + 1
        push!(cells_to_send, cell)
    end

    size = sum(cell.max_bounds .- cell.min_bounds)
    bounds_center = (max_bounds .- min_bounds) ./ 2
    d = distance(cell.cm, bounds_center)

    if size / d >= THETA
        for i in 1:N_CELL_CHILDREN
            if isassigned(cell.children, i) && cell.children[i].mass != 0.0
                get_cells_to_send(cell.children[i], cell, min_bounds, max_bounds, min_depth, depth + 1, cells_to_send)
            else
                break
            end
        end
    end
    cells_to_send
end

@everywhere function cell_cointains_bounds(cell::Cell, min_bounds::Vec3, max_bounds::Vec3)
    for c in 1:DIMS
        if cell.min_bounds[c] > min_bounds[c] || cell.max_bounds[c] < max_bounds[c]
            return false
        end
    end
    true
end

@everywhere function insert_empty_cell(cell::Cell, min_bounds::Vec3, max_bounds::Vec3)
    for i in 1:N_CELL_CHILDREN
        if isassigned(cell.children, i) && cell_cointains_bounds(cell.children[i], min_bounds, max_bounds)
            insert_empty_cell(cell.children[i], min_bounds, max_bounds)
            return
        elseif !isassigned(cell.children, i)
            child = Cell()
            cell.children[i] = child
            child.min_bounds = min_bounds
            child.max_bounds = max_bounds
            child.cell_opened = true
            return
        end
    end
end

@everywhere function insert_cell(cell::Cell, insert::Cell)
    if cell.mass != 0.0 || insert.mass != 0.0
        cell.cm = (cell.mass .* cell.cm .+ insert.mass .* insert.cm) ./ (cell.mass + insert.mass)
        cell.mass += insert.mass
    end

    for i in 1:N_CELL_CHILDREN
        if isassigned(cell.children, i) && cell_cointains_bounds(cell.children[i], insert.min_bounds, insert.max_bounds)
            insert_cell(cell.children[i], insert)
            return
        elseif !isassigned(cell.children, i)
            cell.children[i] = insert
            return
        end
    end
end

@everywhere function insert_body(cell::Cell, body::Body)
    if !cell.body_present && is_external(cell)
        cell.body = body
        cell.mass = body.mass
        cell.cm = copy(body.position)
        cell.body_present = true
    else
        if is_external(cell)
            half_sides = (cell.max_bounds .- cell.min_bounds) ./ 2
            for i in 1:N_CELL_CHILDREN
                child = Cell()
                cell.children[i] = child

                shifts = Vec3([((i - 1) >> (c - 1)) & 1 for c in 1:DIMS])
                child.min_bounds = cell.min_bounds .+ (shifts .* half_sides)
                child.max_bounds = cell.max_bounds .- ((1 .- shifts) .* half_sides)

                # Insert old body
                if (cell.body_present && cell_contains_position(child, cell.body.position))
                    insert_body(child, cell.body)
                    cell.body_present = false
                end
            end
            cell.cell_opened = true
        end

        # Insert new body
        for i in 1:N_CELL_CHILDREN
            if isassigned(cell.children, i) && cell_contains_position(cell.children[i], body.position)
                insert_body(cell.children[i], body)
                break
            end
        end

        # Update mass and center of mass of the cell
        cell.cm = (cell.mass .* cell.cm .+ body.mass .* body.position) ./ (cell.mass + body.mass)
        cell.mass += body.mass
    end
end

@everywhere function add_force(body::Body, position::Vec3, mass::Float64)::Nothing
    dx = position[X] - body.position[X]
    dy = position[Y] - body.position[Y]
    dz = position[Z] - body.position[Z]

    dist = sqrt(dx*dx + dy*dy + dz*dz)
    dist_cubed = dist * dist * dist
    f = (-G * body.mass * mass) / dist_cubed

    body.force = body.force + Vec3(f * dx, f * dy, f * dz)
    nothing
end

@everywhere @propagate_inbounds function compute_force(cell::Cell, body::Body)::Nothing
    cell_body::Body = cell.body
    if (cell.body_present && cell_body.id === body.id) || cell.mass === 0.0
        return
    end
    
    if is_external(cell) && cell.body_present
        @inbounds add_force(body, cell_body.position, cell_body.mass)
    else
        size = sum(cell.max_bounds .- cell.min_bounds)
        d = distance(cell.cm, body.position)
        
        if size / d < THETA
            # Treat cell as a body
            @inbounds add_force(body, cell.cm, cell.mass)
        else
            for i in 1:N_CELL_CHILDREN
                if isassigned(cell.children, i)
                    @inbounds compute_force(cell.children[i], body)
                end
            end
        end
    end
end
