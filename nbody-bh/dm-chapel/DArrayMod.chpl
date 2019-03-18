module DArrayMod {

    use BodyMod;

    type UBodiesArray = unmanaged BodiesArray;
    class BodiesArray {
        var elements: [0..-1] Body;
        proc push(el: Body) {
            this.elements.push_back(el);
        }
    }
}
