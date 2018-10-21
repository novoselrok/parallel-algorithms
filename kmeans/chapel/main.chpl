use Time;

/* Config variables */
config const filename = "../data/input.txt";
config const out_filename = "out.txt";
config const k = 2;
config const max_iter = 10;
config const rows = 0;
config const cols = 0;

/* Read input file */
var f = open(filename, iomode.r);
var reader = f.reader();

/* Define data structures */
type Point = [0..#cols] real;
// type Point = 16 * real;
const pDomain = 0..#rows;
const clustersDomain = 0..#k;
var data: [pDomain] Point;
reader.read(data);

// for i in pDomain {
//     data[i] = reader.read(real, real, real, real, real, real, real, real, real, real, real, real, real, real, real, real);
// }

f.close();
reader.close();

record Cluster {
    var n_points: int;
    var acc: Point;
    var mean: Point;

    proc init() {
        n_points = 0;
        acc = 0.0;
        // acc = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    proc add_point(point: Point) {
        n_points += 1;
        acc += point;
    }

    proc add_cluster(other_cluster: Cluster) {
        n_points += other_cluster.n_points;
        acc += other_cluster.acc;
    }

    proc dist(point: Point) {
        return sqrt(+ reduce ((mean - point) ** 2));
    }

    proc calc_mean() {
        mean = acc / n_points;
    }
}

proc +(ref c1: Cluster, c2: Cluster) {
    c1.add_cluster(c2);
    return c1;
}

var labels: [pDomain] int;
var clusters: [clustersDomain] Cluster;

writeln("-> Starting simulation");

var watch: Timer;
watch.start();

// Initial clusters
forall i in pDomain with (+ reduce clusters) {
    labels[i] = i % k;
    clusters[labels[i]].add_point(data[i]);
}

writeln("-> Init done");

for iteration in 0..#max_iter {
    writeln(iteration);
    [cluster in clusters] cluster.calc_mean();

    var new_clusters: [clustersDomain] Cluster;
    forall i in pDomain with (+ reduce new_clusters) {
        ref point = data[i];
        var dists = [cluster in clusters] cluster.dist(point);
        var (_, min_index) = minloc reduce zip(dists, dists.domain);
        labels[i] = min_index;
        new_clusters[min_index].add_point(point);
    }

    clusters = new_clusters;

    // [cluster in new_clusters] cluster.calc_mean();
    // [cluster in new_clusters] writeln(cluster);
}

writeln("-> The simulation took ", watch.elapsed(), " seconds");

// Output labels
var outf = open(out_filename, iomode.cw);
var writer = outf.writer();

for lbl in labels {
  writer.writeln(lbl);
}

writer.close();
outf.close();
