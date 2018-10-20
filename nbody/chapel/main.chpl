use Time;

config const filename = "input.txt";
config const iterations = 1;
config const out_filename = "out.txt";
const X = 1;
const Y = 2;
const Z = 3;
const G = 6.67e-11;
config const dt = 0.1;

// Read input file, initialize bodies                                       
var f = open(filename, iomode.r);
var reader = f.reader();

var n_bodies = reader.read(int);
const pDomain = {0..#n_bodies};

type vec3 = 3*real;

var forces: [pDomain] vec3;
var velocities: [pDomain] vec3;
var positions: [pDomain] vec3;
var masses: [pDomain] real;

for i in pDomain {
  positions[i] = reader.read(real, real, real);
  velocities[i] = reader.read(real, real, real);
  masses[i] = reader.read(real);
}

f.close();
reader.close();
writeln("Done reading...");

var watch: Timer;
watch.start();
for i in 0..#iterations {
  // Reset forces                                                           
  forces = (0.0, 0.0, 0.0);
  forall q in pDomain with (+ reduce forces) {
    for k in pDomain {
      if k <= q {
        continue;
      }
      var diff = positions[q] - positions[k];
      var dist = sqrt(diff[X]**2 + diff[Y]**2 + diff[Z]**2);
      var dist_cubed = dist**3;
      var tmp = -G * masses[q] * masses[k] / dist_cubed;
      var force_qk = tmp * diff;

      forces[q] += force_qk;
      forces[k] -= force_qk;
    }
  }
  
  forall q in pDomain {
    positions[q] += dt * velocities[q];
    velocities[q] += dt / masses[q] * forces[q];
  }
}
writeln("The simulation took ", watch.elapsed(), " seconds");

var outf = open(out_filename, iomode.cw);
var writer = outf.writer();

for q in pDomain {
  writer.writeln("%er %er %er %er %er %er".format(positions[q][X], positions[q][Y], positions[q][Z], velocities[q][X], velocities[q][Y], velocities[q][Z]));
}

writer.close();
outf.close();
