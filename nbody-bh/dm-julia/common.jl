using Distributed
@everywhere using StaticArrays

@everywhere const X = 1
@everywhere const Y = 2
@everywhere const Z = 3
@everywhere const DIMS = 3
@everywhere const THETA = 0.5
@everywhere const G = 6.67e-11

# Custom types
@everywhere const Vec3 = SVector{DIMS, Float64}

@everywhere function distance(vec1::Vec3, vec2::Vec3)::Float64
    dx::Float64 = vec1[X] - vec2[X]
    dy::Float64 = vec1[Y] - vec2[Y]
    dz::Float64 = vec1[Z] - vec2[Z]

    sqrt(dx*dx + dy*dy + dz*dz)
end
