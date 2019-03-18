module BodyMod {
    use Common;

    record Body {
        var id: int;
        var force: Vec3;
        var position: Vec3;
        var velocity: Vec3;
        var mass: real;
        var work: real;

        proc resetForce() {
            force = (0.0, 0.0, 0.0);
        }
    }
}
