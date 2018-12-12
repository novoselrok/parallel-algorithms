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

    proc distance(p: Point) {
        return sqrt(+ reduce (mean - p) ** 2);
    }

    proc addPoint(p: Point) {
        pointSum += p;
        size += 1;
    }

    proc setMean(p: Point) {
        mean = p;
    }

    proc calcMean() {
        mean = pointSum / size;
    }
}

proc main() {
    var pointsDomain = {0..#numPoints};
    var points: [pointsDomain] Point;
    var clustersDomain = {0..#k};
    var clusters: [clustersDomain] Cluster;
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
    var randStream = new owned RandomStream(real);
    for i in clustersDomain {
        var idx: int = (randStream.getNext() * numPoints):int;
        clusters[i].setMean(points[idx]);
    }

    for iteration in 0..#maxIter {
        var newClusters: [clustersDomain] Cluster;
        for i in pointsDomain {
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
