use IO;
use Math;

writeln(here.maxTaskPar);

config const filename = "input.txt";
config const iterations = 100;
config const out_filename = "out.txt";
const X = 0;
const Y = 1;
const Z = 2;
const G = 6.67e-11;
config const dt = 0.1;

// Read input file, initialize bodies
var f = open(filename, iomode.r);
var reader = f.reader();

var n_bodies = reader.read(int);
const pDomain = {0..#n_bodies};

var forces: [0..#n_bodies, 0..#3] real;
var velocities: [0..#n_bodies, 0..#3] real;
var positions: [0..#n_bodies, 0..#3] real;
var masses: [0..#n_bodies] real;

for i in 0..#n_bodies {
    positions[i, X] = reader.read(real);
    positions[i, Y] = reader.read(real);
    positions[i, Z] = reader.read(real);

    velocities[i, X] = reader.read(real);
    velocities[i, Y] = reader.read(real);
    velocities[i, Z] = reader.read(real);

    masses[i] = reader.read(real);
}

f.close();
reader.close();

for i in 0..#iterations {
    // Reset forces
    forces = 0.0;

    for q in 0..#n_bodies {
        for k in 0..#n_bodies {
            if k <= q {
                continue;
            }
            var x_diff = positions[q, X] - positions[k, X];
            var y_diff = positions[q, Y] - positions[k, Y];
            var z_diff = positions[q, Z] - positions[k, Z];
            var dist = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
            var dist_cubed = dist * dist * dist;
            
            var tmp = -G * masses[q] * masses[k] / dist_cubed;
            var force_qk_x = tmp * x_diff;
            var force_qk_y = tmp * y_diff;
            var force_qk_z = tmp * z_diff;

            forces[q, X] += force_qk_x;
            forces[q, Y] += force_qk_y;
            forces[q, Z] += force_qk_z;
            forces[k, X] -= force_qk_x;
            forces[k, Y] -= force_qk_y;
            forces[k, Z] -= force_qk_z;
        }
    }

    for q in 0..#n_bodies {
        positions[q, X] += dt * velocities[q, X];
        positions[q, Y] += dt * velocities[q, Y];
        positions[q, Z] += dt * velocities[q, Z];
        velocities[q, X] += dt / masses[q] * forces[q, X];
        velocities[q, Y] += dt / masses[q] * forces[q, Y];
        velocities[q, Z] += dt / masses[q] * forces[q, Z];
    }
}

var outf = open(out_filename, iomode.cw);
var writer = outf.writer();

for q in 0..#n_bodies {
    writer.writeln("%er %er %er %er %er %er".format(positions[q, X], positions[q, Y], positions[q, Z], velocities[q, X], velocities[q, Y], velocities[q, Z]));
}

writer.close();
outf.close();
