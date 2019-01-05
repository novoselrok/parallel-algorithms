use Random;
use Time;

/* Config variables */
config const filename = "../data/test_10_10_2.chpl.txt";
config const outFilename = "out.txt";
config const k = 2;
config const maxIter = 10;
config const numPoints = 10;

/* Datatypes */
type Point = 10 * real;

record Cluster {
    var size: int;
    var pointSum: Point;
    var mean: Point;

    proc init() {}

    proc init(point: Point) {
        mean = point;
    }

    proc distance(p: Point) {
        return sqrt(+ reduce (mean - p) ** 2);
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

proc kmeansppInit(points: [?] Point, k: int) {
    var clusters: [0..#1] Cluster = [new Cluster(points[0])];
    var randStream = new owned RandomStream(real);

    for clusterIdx in 1..k - 1 {
        var distances = [point in points] (min reduce [cluster in clusters] cluster.distance(point));
        var distancesSum = + reduce distances;
        distances = distances / distancesSum;

        for i in 1..distances.size - 1 {
            distances[i] += distances[i - 1];
        }

        var random = randStream.getNext();
        var newIdx = 0;
        for (i, prob) in zip(distances.domain, distances) {
            if (random < prob) {
                newIdx = i;
                break;
            }
        }
        clusters.push_back(new Cluster(points[newIdx]));
    }

    return clusters;
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
    var clusters = kmeansppInit(points, k);

    for iteration in 0..#maxIter {
        writeln(iteration);
        var newClusters: [clustersDomain] Cluster;
        forall i in pointsDomain with (+ reduce newClusters) {
            // You have to be careful with records since they create copies
            ref point = points[i];
            var distances = [cluster in clusters] cluster.distance(point);

            var (_, minIndex) = minloc reduce zip(distances, distances.domain);
            labels[i] = minIndex;
            newClusters[minIndex].addPoint(point);
        }

        [cluster in newClusters] cluster.calcMean();
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
