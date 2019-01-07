use Random;
use Time;

/* Config variables */
config const filename = "../data/test_10_10_2.chpl.txt";
config const outFilename = "out.txt";
config const k = 2;
config const maxIter = 10;
config const numPoints = 10;

/* Datatypes */
type Point = 100 * real;

record Cluster {
    var size: int;
    var pointSum: Point;
    var mean: Point;

    proc init() {}

    proc init(point: Point) {
        mean = point;
    }

    proc distance(p: Point) {
        var squaredDiff = (mean - p) ** 2;
        var sum = 0.0;
        for el in squaredDiff {
            sum += el;
        }
        return sqrt(sum);
    }

    // Used when reducing clusters
    proc addCluster(other: Cluster) {
        size += other.size;
        pointSum += other.pointSum;
    }

    proc addPoint(ref p: Point) {
        pointSum += p;
        size += 1;
    }

    proc setMean(ref p: Point) {
        mean = p;
    }

    proc calcMean() {
        mean = pointSum / size;
    }
}

// Used when reducing clusters
proc +(ref c1: Cluster, c2: Cluster) {
    c1.addCluster(c2);
    return c1;
}

proc main() {
    var pointsDomain = {0..#numPoints};
    var points: [pointsDomain] Point;
    var clustersDomain = {0..#k};
    var labels: [pointsDomain] int;
    var clusters: [clustersDomain] Cluster;
    
    /* Read input file */
    var f = open(filename, iomode.r);
    var reader = f.reader();
    for i in pointsDomain {
        points[i] = reader.read(Point);
    }
    f.close();
    reader.close();

    var watch: Timer;
    watch.start();
    /* Algorithm */
    var randStream = new owned RandomStream(real);
    for i in clustersDomain {
        var idx: int = (randStream.getNext() * numPoints):int;
        clusters[i] = new Cluster(points[idx]);
    }

    for iteration in 0..#maxIter {
        var newClusters: [clustersDomain] Cluster;
        forall i in pointsDomain with (+ reduce newClusters) {
            // You have to be careful with tuples since they create copies
            ref point = points[i];

            var minIndex = 0;
            var minValue = Math.INFINITY;
            for (idx, cluster) in zip(clusters.domain, clusters) {
                var distance = cluster.distance(point);
                if (distance < minValue) {
                    minValue = distance;
                    minIndex = idx;
                }
            }
            
            labels[i] = minIndex;
            newClusters[minIndex].addPoint(point);
        }

        for cluster in newClusters {
            cluster.calcMean();
        }
        clusters = newClusters;
    }
    writeln(watch.elapsed());

    // Output labels
    var outf = open(outFilename, iomode.cw);
    var writer = outf.writer();
    for lbl in labels {
        writer.writeln(lbl);
    }
    writer.close();
    outf.close();
}
