use Math;
use Time;
use BlockDist;
use CommDiagnostics;

use Common;
use BodyMod;
use CellMod;
use Worker;
use DArrayMod;

config const filename = "../../nbody/data/small.chpl.txt";
config const iterations = 1;
config const nBodies = 16;

proc main() {
    var bodySpace = {0..#nBodies};
    var bodyDomain = bodySpace dmapped Block(boundingBox=bodySpace);
    var bodies: [bodyDomain] Body;

    var f = open(filename, iomode.r);
    var reader = f.reader();
    for i in bodyDomain {
        bodies[i] = new Body();
        bodies[i].id = i;
        bodies[i].position = reader.read(Vec3);
        bodies[i].velocity = reader.read(Vec3);
        bodies[i].mass = reader.read(real);
        bodies[i].work = 1.0;
    }
    f.close();
    reader.close();

    var watch: Timer;
    watch.start();
    var dt = 0.1;

    var localeSpace = {0..#numLocales};
    var localeDomain = localeSpace dmapped Block(boundingBox=localeSpace);

    var bounds: [localeDomain] Bounds;
    var bounds$: [localeDomain] sync int;

    var workSum: [localeDomain] Work;
    var workSum$: [localeDomain] sync int;

    var otherBodies: [localeDomain] owned BodiesArray;
    var otherBodies$: [localeDomain] sync int;

    var groups: [localeDomain] int;
    var groups$: [localeDomain] sync int;

    var cells: [localeDomain] owned CellTupleArray;
    var cells$: [localeDomain] sync int;

    // startVerboseComm();
    coforall L in Locales do on L {
        work(
            bodies.localSlice(bodies.localSubdomain()),
            bounds,
            bounds$,
            workSum,
            workSum$,
            otherBodies,
            otherBodies$,
            groups,
            groups$,
            cells,
            cells$,
            iterations
        );
    }
    // stopVerboseComm();

    writeln(watch.elapsed());

    // var outf = open("out.txt", iomode.cw);
    // var writer = outf.writer();

    // for body in bodies {
    //     writer.writeln("%er %er %er %er %er %er".format(body.position[X], body.position[Y], body.position[Z], body.velocity[X], body.velocity[Y], body.velocity[Z]));
    // }

    // writer.close();
    // outf.close();
}
