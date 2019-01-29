use Random;
use Time;
use BlockDist;
use CyclicDist;
use ReplicatedDist;
use CommDiagnostics;

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
    // startVerboseComm();
    var pointsSpace = {0..#numPoints};
    var pointsDomain = pointsSpace dmapped Block(boundingBox=pointsSpace);
    var points: [pointsDomain] Point;
    var labels: [pointsDomain] int;
    
    var clustersSpace = {0..#k};
    var clustersDomain = clustersSpace dmapped Replicated();
    var clusters: [clustersDomain]  Cluster;
    
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

    for L in Locales {
        clusters.replicand(L) = clusters;
    }

    for iteration in 0..#maxIter {
        var newClustersPerLocale: [clustersDomain] Cluster;
        forall (point, label_) in zip(points, labels) {
            var minIndex = 0;
            var minValue = Math.INFINITY;
            for (idx, cluster) in zip(clusters.domain, clusters) {
                var distance = cluster.distance(point);
                if (distance < minValue) {
                    minValue = distance;
                    minIndex = idx;
                }
            }
            
            label_ = minIndex;
            newClustersPerLocale[minIndex].addPoint(point);
        }
        var newClusters = + reduce [L in Locales] newClustersPerLocale.replicand(L);
        for cluster in newClusters {
            cluster.calcMean();
        }
        clusters = newClusters;
        
        for L in Locales {
            clusters.replicand(L) = clusters;
        }
    }
    // stopVerboseComm();
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
