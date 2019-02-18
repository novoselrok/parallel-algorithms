use Math;
use Time;

use Common;
use BodyMod;
use CellMod;

config const filename = "../data/small.chpl.txt";
config const iterations = 1;
config const nBodies = 16;

proc getUniverseSize(bodies: [] UBody) {
    var universeMin = [Math.INFINITY, Math.INFINITY, Math.INFINITY];
    var universeMax = [-Math.INFINITY, -Math.INFINITY, -Math.INFINITY];

    for body in bodies {
        for c in 1..DIMS {
            if body.position[c] < universeMin[c] {
                universeMin[c] = body.position[c];
            }
            if body.position[c] > universeMax[c] {
                universeMax[c] = body.position[c];
            }
        }
    }

    return (
        (universeMin[X], universeMin[Y], universeMin[Z]), 
        (universeMax[X], universeMax[Y], universeMax[Z])
    );
}

proc main() {
    var bodyDomain = {0..#nBodies};
    var bodies: [bodyDomain] UBody;

    var f = open(filename, iomode.r);
    var reader = f.reader();
    for i in bodyDomain {
        var body = new UBody();
        bodies[i] = body;
        body.id = i;
        body.position = reader.read(Vec3);
        body.velocity = reader.read(Vec3);
        body.mass = reader.read(real);
    }
    f.close();
    reader.close();

    var watch: Timer;
    watch.start();
    var dt = 0.1;
    for iteration in 0..#iterations {
        var (universeMin, universeMax) = getUniverseSize(bodies);

        for body in bodies {
            body.resetForce();
        }

        var root: unmanaged Cell = new unmanaged Cell();
        root.minBounds = universeMin;
        root.maxBounds = universeMax;

        for body in bodies {
            insertBody(root, body);
        }
        writeln(watch.elapsed());

        for body in bodies {
            computeForce(root, body);
        }
        writeln(watch.elapsed());

        for body in bodies {
            body.position = body.position + (dt * body.velocity);
            body.velocity = body.velocity + (dt / body.mass * body.force);
        }
    }
    writeln(watch.elapsed());

    var outf = open("out.txt", iomode.cw);
    var writer = outf.writer();

    for body in bodies {
        writer.writeln("%er %er %er %er %er %er".format(body.position[X], body.position[Y], body.position[Z], body.velocity[X], body.velocity[Y], body.velocity[Z]));
    }

    writer.close();
    outf.close();
}
