include("common.jl")
include("body.jl")

import Base: @propagate_inbounds

const N_CELL_CHILDREN = 8

mutable struct Cell
    cell_opened::Bool # If true it means it has children
    children::Array{Cell, 1}
    body_present::Bool
    body::Body
    mass::Float64
    cm::Vec3
    min_bounds::Vec3
    max_bounds::Vec3
end

Cell() = Cell(
    false,
    Array{Cell, 1}(undef, N_CELL_CHILDREN),
    false,
    Body(),
    0.0, 
    Vec3(zeros(DIMS)),
    Vec3(zeros(DIMS)),
    Vec3(zeros(DIMS))
)

function is_external(cell::Cell)
    !cell.cell_opened
end

function init_cell_bounds_from_parent(parent::Cell, child::Cell, half_sides::Vec3, child_idx::Int)
    shifts = Vec3([((child_idx - 1) >> (c - 1)) & 1 for c in 1:DIMS])
    child.min_bounds = parent.min_bounds .+ (shifts .* half_sides)
    child.max_bounds = parent.max_bounds .- ((1 .- shifts) .* half_sides)
end

construct_empty_tree(cell::Cell, max_level::Int) = construct_empty_tree(cell, 0, max_level)

function construct_empty_tree(cell::Cell, current_level::Int, max_level::Int)
    current_level === max_level && return
    cell.cell_opened = true
    half_sides = (cell.max_bounds .- cell.min_bounds) ./ 2
    for i in 1:N_CELL_CHILDREN
        child = Cell()
        cell.children[i] = child
        init_cell_bounds_from_parent(cell, child, half_sides, i)
        construct_empty_tree(child, current_level + 1, max_level)
    end
end

get_leaves(cell::Cell, max_level::Int) = get_leaves(cell, Array{Cell}([]), 0, max_level)

function get_leaves(cell::Cell, leaves::Array{Cell}, current_level::Int, max_level::Int)::Array{Cell}
    if current_level === max_level - 1
        for child in cell.children
            push!(leaves, child)
        end
    else
        for child in cell.children
            get_leaves(child, leaves, current_level + 1, max_level)
        end
    end
    leaves
end

update_empty_cells(cell::Cell, max_level::Int) = update_empty_cells(cell, 0, max_level)

function update_empty_cells(cell::Cell, current_level::Int, max_level::Int)
    current_level === max_level && return

    for child in cell.children
        update_empty_cells(child, current_level + 1, max_level)
    end

    mass = 0.0
    for child in cell.children
        mass += child.mass

        combined_mass = cell.mass + child.mass
        combined_mass === 0.0 && continue

        cell.cm = (cell.mass .* cell.cm .+ child.mass .* child.cm) ./ combined_mass
    end
    cell.mass = mass
end

function cell_contains_position(cell::Cell, position::Vec3)
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

function insert_body(cell::Cell, body::Body)
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
                init_cell_bounds_from_parent(cell, child, half_sides, i)

                # Insert old body
                if (cell.body_present && cell_contains_position(child, cell.body.position))
                    insert_body(child, cell.body)
                    cell.body_present = false
                end
            end
            cell.cell_opened = true
        end

        # Insert new body
        for child in cell.children
            if cell_contains_position(child, body.position)
                insert_body(child, body)
                break
            end
        end

        # Update mass and center of mass of the cell
        cell.cm = (cell.mass .* cell.cm .+ body.mass .* body.position) ./ (cell.mass + body.mass)
        cell.mass += body.mass
    end
end

function add_force(body::Body, position::Vec3, mass::Float64)::Nothing
    dx = position[X] - body.position[X]
    dy = position[Y] - body.position[Y]
    dz = position[Z] - body.position[Z]

    dist = sqrt(dx*dx + dy*dy + dz*dz)
    dist_cubed = dist * dist * dist
    f = (-G * body.mass * mass) / dist_cubed

    body.force = body.force + Vec3(f * dx, f * dy, f * dz)
    nothing
end

@propagate_inbounds function compute_force(cell::Cell, body::Body)::Nothing
    cell_body::Body = cell.body
    if (cell.body_present && cell_body.id === body.id) || cell.mass === 0.0
        return
    end
    
    if is_external(cell)
        @inbounds add_force(body, cell_body.position, cell_body.mass)
    else
        size = sum(cell.max_bounds .- cell.min_bounds)
        d = distance(cell.cm, body.position)
        
        if size / d < THETA
            # Treat cell as a body
            @inbounds add_force(body, cell.cm, cell.mass)
        else
            for child in cell.children
                @inbounds compute_force(child::Cell, body::Body)
            end
        end
    end
end
