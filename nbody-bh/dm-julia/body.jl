include("common.jl")

@everywhere mutable struct Body
    id::Int
    force::Vec3
    position::Vec3
    velocity::Vec3
    mass::Float64
    work::Float64

    Body() = new()
    Body(id, force, position, velocity, mass, work) = new(id, force, position, velocity, mass, work)
end

@everywhere Body(id::Int, position::Vec3, velocity::Vec3, mass::Float64, work::Float64) = Body(
    id,
    Vec3(zeros(DIMS)),
    position,
    velocity,
    mass,
    work
)

@everywhere function reset_force(body::Body)
    body.force = Vec3(0.0, 0.0, 0.0)
end
