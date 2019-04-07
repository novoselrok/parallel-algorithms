use Random;
proc computePi(n: int) {
    var withinCircle = 0;
    var randStream = new owned RandomStream(real);
    for _i in 0..#n {
        var x = randStream.getNext() * 2 - 1;
        var y = randStream.getNext() * 2 - 1;
        var r2 = x * x + y * y;
        if r2 < 1.0 {
            withinCircle += 1;
        }
    }
    return withinCircle / n:real * 4.0;
}
param N = 100000;
var result = 0.0;
forall _i in 0..#here.maxTaskPar with (+ reduce result) {
    result = computePi(Math.ceil(N / here.maxTaskPar):int);
}
writeln(result / here.maxTaskPar);

