module Common {
    type Vec3 = (real, real, real);

    param X = 1;
    param Y = 2;
    param Z = 3;
    param DIMS = 3;
    param THETA = 0.5;
    param G = 6.67e-11;

    proc distance(vec1: Vec3, vec2: Vec3) {
        var dx = vec1[X] - vec2[X];
        var dy = vec1[Y] - vec2[Y];
        var dz = vec1[Z] - vec2[Z];

        return sqrt(dx*dx + dy*dy + dz*dz);
    }
}