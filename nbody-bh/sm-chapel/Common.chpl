module Common {
    type Vec3 = (real, real, real);

    const X = 1;
    const Y = 2;
    const Z = 3;
    const DIMS = 3;
    const THETA = 0.5;
    const G = 6.67e-11;

    proc distance(vec1: Vec3, vec2: Vec3) {
        var dx = vec1[X] - vec2[X];
        var dy = vec1[Y] - vec2[Y];
        var dz = vec1[Z] - vec2[Z];

        return sqrt(dx*dx + dy*dy + dz*dz);
    }
}