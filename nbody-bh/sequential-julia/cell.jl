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
        cell.body_present =  true
    else
        if is_external(cell)
            for i in 1:N_CELL_CHILDREN
                child = Cell()
                cell.children[i] = child

                # Find the bounds of the child
                for c in 1:DIMS
                    half_side = (cell.max_bounds[c] - cell.min_bounds[c]) / 2
                    shift::Int = ((i - 1) >> (c - 1)) & 1
                    child.min_bounds[c] = cell.min_bounds[c] + (shift * half_side)
                    child.max_bounds[c] = cell.max_bounds[c] - ((1 - shift) * half_side)
                end

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
        for c in 1:DIMS
            cell.cm[c] = (cell.mass * cell.cm[c] + body.mass * body.position[c]) / (cell.mass + body.mass)
        end
        cell.mass += body.mass
    end
end

function add_force(body::Body, position::Vec3, mass::Float64)::Nothing
    dx::Float64 = position[X] - body.position[X]
    dy::Float64 = position[Y] - body.position[Y]
    dz::Float64 = position[Z] - body.position[Z]

    dist::Float64 = sqrt(dx*dx + dy*dy + dz*dz)
    dist_cubed::Float64 = dist * dist * dist
    f::Float64 = (-G * body.mass * mass) / dist_cubed

    body.force[X] += f * dx;
    body.force[Y] += f * dy;
    body.force[Z] += f * dz;
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
        dx::Float64 = cell.max_bounds[X] - cell.min_bounds[X]
        dy::Float64 = cell.max_bounds[Y] - cell.min_bounds[Y]
        dz::Float64 = cell.max_bounds[Z] - cell.min_bounds[Z]
        size::Float64 = dx + dy + dz
        d::Float64 = distance(cell.cm, body.position)
        
        if size / d < THETA
            # Treat cell as a body
            @inbounds add_force(body, cell.cm, cell.mass)
        else
            for child in cell.children
                compute_force(child::Cell, body::Body)
            end
        end
    end
end
