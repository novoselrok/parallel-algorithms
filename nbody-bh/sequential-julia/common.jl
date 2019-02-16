using StaticArrays

const X = 1
const Y = 2
const Z = 3

const DIMS = 3

const Vec3 = MVector{DIMS, Float64}

const THETA = 0.5

const G = 6.67e-11

function distance(vec1::Vec3, vec2::Vec3)::Float64
    dx::Float64 = vec1[X] - vec2[X]
    dy::Float64 = vec1[Y] - vec2[Y]
    dz::Float64 = vec1[Z] - vec2[Z]

    sqrt(dx*dx + dy*dy + dz*dz)
end
