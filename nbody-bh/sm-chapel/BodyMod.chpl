module BodyMod {
    use Common;

    type UBody = unmanaged Body;
    
    class Body {
        var id: int;
        var force: Vec3;
        var position: Vec3;
        var velocity: Vec3;
        var mass: real;

        proc resetForce() {
            force = (0.0, 0.0, 0.0);
        }
    }
}
