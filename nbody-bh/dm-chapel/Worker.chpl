module Worker {
    use Math;
    use Sort;
    use IO;
    use Time;

    use Common;
    use BodyMod;
    use CellMod;
    use DArrayMod;

    type Bounds = (Vec3, Vec3);
    type PartnerInfo = (int, bool);
    type Work = (real, real);

    param BISECTION_MAX_ITER = 200;
    param BISECTION_TOL = 1e-10;

    proc getUniverseSize(bodies: [] Body) {
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

    proc determineGroupPartners(group: int, groups: [] int, groups$: [] sync int) {
        var myGroupPartners: [0..-1] int;
        var rank = here.id;

        myGroupPartners.push_back(rank);
        for i in 0..#groups.size-1 {
            var intendedReceiver = (rank + i + 1) % groups.size;
            var partner = ((rank - (i + 1)) + groups.size) % groups.size;

            // Wait until the previous value is read
            if (groups$[rank].isFull) {
                groups$[rank];
            }

            // Write values
            groups[rank] = group;
            groups$[rank] = intendedReceiver;

            // Am I the intended receiver?
            while (groups$[partner].readFF() != rank) {}

            // Read partner value
            var partnerGroup = groups[partner];
            if partnerGroup == group {
                myGroupPartners.push_back(partner);
            }

            groups$[partner]; // empty
            // Reset write, which blocks until variable is empty
            groups$[rank] = -1;
        }

        sort(myGroupPartners);
        return myGroupPartners;
    }

    proc fracWeightBelow(
        bodies: [] Body,
        split: real,
        coord: int,
        workSum: [] Work,
        workSum$: [] sync int,
        groupPartners: [] int
    ) {
        var rank = here.id;
        var (_, idxWithinGroup) = groupPartners.find(rank);
        var workBelow = 0.0;
        var workAbove = 0.0;

        for body in bodies {
            if body.position[coord] > split {
                workAbove += body.work;
            } else {
                workBelow += body.work;
            }
        }

        var allWorkBelow = workBelow;
        var allWorkAbove = workAbove;

        for i in 0..#groupPartners.size-1 {
            var intendedReceiverIdx = (idxWithinGroup + i + 1) % groupPartners.size;
            var intendedReceiver = groupPartners[intendedReceiverIdx];
            var partnerIdx = ((idxWithinGroup - (i + 1)) + groupPartners.size) % groupPartners.size;
            var partner = groupPartners[partnerIdx];

            // Wait until the previous value is read
            if (workSum$[rank].isFull) {
                workSum$[rank];
            }

            // Write values
            workSum[rank] = (workBelow, workAbove);
            workSum$[rank] = intendedReceiver;

            // Am I the intended receiver?
            while (workSum$[partner].readFF() != rank) {}

            // Read partner value
            var (otherWorkBelow, otherWorkAbove) = workSum[partner];
            allWorkBelow += otherWorkBelow;
            allWorkAbove += otherWorkAbove;

            workSum$[partner]; // empty
            // Reset write, which blocks until variable is empty
            workSum$[rank] = -1;
        }

        var fraction = allWorkBelow / (allWorkBelow + allWorkAbove);
        return fraction - 0.5;
    }

    proc bisection(
        _min: real,
        _max: real,
        bodies: [] Body,
        coord: int,
        workSum: [] Work,
        workSum$: [] sync int,
        groupPartners: [] int
    ) {
        var min = _min;
        var max = _max;
        var fmin = fracWeightBelow(bodies, min, coord, workSum, workSum$, groupPartners);
        var mid = 0.0;

        var iteration = 0;
        while iteration < BISECTION_MAX_ITER && abs((max - min) / 2) > BISECTION_TOL {
            mid = (min + max) / 2;
            var fmid = fracWeightBelow(bodies, mid, coord, workSum, workSum$, groupPartners);

            if abs(fmid) < BISECTION_TOL {
                break;
            } else if (fmin * fmid > 0) {
                min = mid;
            } else {
                max = mid;
            }

            iteration += 1;
        }

        return mid;
    }

    proc setIndexVec3(vec: Vec3, idx: int, value: real) {
        var (x, y, z) = vec;

        if idx == 1 {
            return (value, y, z);
        } else if idx == 2 {
            return (x, value, z);
        } else {
            return (x, y, value);
        }
    }

    proc isAboveSplit(rank: int, nProcsLeft: int) {
        var relBit = Math.log2(nProcsLeft):int;
        var isAbove = (rank >> (relBit - 1)) & 1;
        return isAbove == 1;
    }

    proc getPartnerRank(rank: int, nProcsLeft: int) {
        var relBit = Math.log2(nProcsLeft):int;
        return rank ^ (1 << (relBit - 1));
    }

    proc orb(
        _bodies: [] Body,
        universeMin: Vec3,
        universeMax: Vec3,
        workSum: [] Work,
        workSum$: [] sync int,
        otherBodiesWrappers: [] owned BodiesArray,
        otherBodiesWrappers$: [] sync int,
        groups: [] int,
        groups$: [] sync int
    ) {
        var rank = here.id;
        var myMin = universeMin;
        var myMax = universeMax;

        var bodies: [0..-1] Body;
        bodies.push_back(_bodies);

        var myBounds: [0..-1] Bounds;
        var otherBounds: [0..-1] Bounds;
        var partners: [0..-1] PartnerInfo;

        var group = 0;
        var aboveSplit = false;

        var nSplits = Math.log2(numLocales):int;
        for i in 0..#nSplits {
            var nProcsLeft = (numLocales / (2**i)):int;

            group = group << 1;
            if aboveSplit {
                group = group | 1;
            }

            var groupPartners = determineGroupPartners(group, groups, groups$);

            var coord = (i % DIMS) + 1;

            var split = bisection(myMin[coord], myMax[coord], bodies, coord, workSum, workSum$, groupPartners);

            aboveSplit = isAboveSplit(rank, nProcsLeft);

            var otherMin = myMin;
            var otherMax = myMax;

            if aboveSplit {
                myMin = setIndexVec3(myMin, coord, split);
                otherMax = setIndexVec3(otherMax, coord, split);
            } else {
                myMax = setIndexVec3(myMax, coord, split);
                otherMin = setIndexVec3(otherMin, coord, split);
            }

            myBounds.push_back((myMin, myMax));
            otherBounds.push_back((otherMin, otherMax));

            var myBodies: [0..-1] Body;
            var otherBodies: [0..-1] Body;

            for body in bodies {
                if (body.position[coord] - split > 0) == aboveSplit {
                    myBodies.push_back(body);
                } else {
                    otherBodies.push_back(body);
                }
            }

            var partner = getPartnerRank(rank, nProcsLeft);
            partners.push_back((partner, aboveSplit));

            if (otherBodiesWrappers$[rank].isFull) {
                otherBodiesWrappers$[rank];
            }

            var otherBodiesWrapper = new owned BodiesArray();
            otherBodiesWrapper.elements.push_back(otherBodies);
            otherBodiesWrappers[rank] = otherBodiesWrapper;
            otherBodiesWrappers$[rank] = partner;

            // Am I the intended receiver?
            while (otherBodiesWrappers$[partner].readFF() != rank) {}

            // Read the bodies
            bodies.clear();
            bodies.push_back(myBodies);
            var partnersOtherBodies = otherBodiesWrappers[partner];
            bodies.push_back(partnersOtherBodies.elements);

            otherBodiesWrappers$[partner]; // empty
            // Reset write, which blocks until variable is empty
            otherBodiesWrappers$[rank] = -1;
        }
        return (bodies, myBounds, otherBounds, partners);
    }

    proc work(
        inputBodies: [] Body, 
        bounds: [] Bounds,
        bounds$: [] sync int,
        workSum: [] Work,
        workSum$: [] sync int,
        otherBodies: [] owned BodiesArray,
        otherBodies$: [] sync int,
        groups: [] int,
        groups$: [] sync int,
        cells: [] owned CellTupleArray,
        cells$: [] sync int,
        iterations: int
    ) {
        var bodies: [0..-1] Body;
        var rank = here.id;
        var dt = 0.1;

        for inputBody in inputBodies {
            bodies.push_back(inputBody);
        }

        for iteration in 0..#iterations {
            try! stdout.flush();
            var myUniverseBounds = getUniverseSize(bodies);
            var (universeMin, universeMax) = myUniverseBounds;

            for i in 0..#bounds.size-1 {
                var intendedReceiver = (rank + i + 1) % numLocales;
                var partner = ((rank - (i + 1)) + numLocales) % numLocales;

                // Wait until the previous value is read
                if (bounds$[rank].isFull) {
                    bounds$[rank];
                }

                // Read my value
                bounds[rank] = myUniverseBounds;
                bounds$[rank] = intendedReceiver;

                // Am I the intended receiver?
                while (bounds$[partner].readFF() != rank) {}

                // Read partner value
                var (minOtherBounds, maxOtherBounds) = bounds[partner];
                universeMin = (
                    min(universeMin[X], minOtherBounds[X]),
                    min(universeMin[Y], minOtherBounds[Y]),
                    min(universeMin[Z], minOtherBounds[Z])
                );
                universeMax = (
                    max(universeMax[X], maxOtherBounds[X]),
                    max(universeMax[Y], maxOtherBounds[Y]),
                    max(universeMax[Z], maxOtherBounds[Z])
                );

                bounds$[partner]; // empty
                // Reset write, which blocks until bounds$[rank] is empty
                bounds$[rank] = -1;
            }

            var (orbBodies, myBounds, otherBounds, partners) = orb(
                bodies,
                universeMin,
                universeMax,
                workSum,
                workSum$,
                otherBodies,
                otherBodies$,
                groups,
                groups$
            );
            bodies.clear();
            bodies.push_back(orbBodies);

            for body in bodies {
                body.resetForce();
            }

            var root: unmanaged Cell = new unmanaged Cell();
            root.minBounds = universeMin;
            root.maxBounds = universeMax;

            for (minBounds, maxBounds) in myBounds {
                insertEmptyCell(root, minBounds, maxBounds);
            }

            for body in bodies {
                insertBody(root, body);
            }

            for i in 0..#myBounds.size {
                var (myMin, myMax) = myBounds[i];
                var (otherMin, otherMax) = otherBounds[i];
                
                var cellsToSend: [0..-1] UCell;  
                getCellsToSend(root, nil, otherMin, otherMax, i, 0, cellsToSend);

                var packed = packCells(cellsToSend);
                var (partner, aboveSplit) = partners[i];
                
                // Wait until the previous value is read
                if (cells$[rank].isFull) {
                    cells$[rank];
                }

                var cellTuplesWrapper = new owned CellTupleArray();
                cellTuplesWrapper.elements.push_back(packed);
                cells[rank] = cellTuplesWrapper;
                cells$[rank] = partner;

                // Am I the intended receiver?
                while (cells$[partner].readFF() != rank) {}

                var recvCells: [0..-1] CellTuple;
                recvCells.push_back(cells[partner].elements);
                var rootCells = reconstructReceivedCells(recvCells);
                for rootCell in rootCells {
                    insertCell(root, rootCell);
                }

                cells$[partner]; // empty
                // Reset write, which blocks until variable is empty
                cells$[rank] = -1;
            }

            for body in bodies {
                var watch: Timer;
                watch.start();
                computeForce(root, body);
                body.work = watch.elapsed();
            }

            for body in bodies {
                body.position = body.position + (dt * body.velocity);
                body.velocity = body.velocity + (dt / body.mass * body.force);
            }

            freeCell(root);
        }

        // var filename = "out" + here.id + ".txt";
        // try {
        //     var outf = open(filename, iomode.cw);
        //     var writer = outf.writer();

        //     for body in bodies {
        //         writer.writeln("%i %er %er %er %er %er %er".format(body.id, body.position[X], body.position[Y], body.position[Z], body.velocity[X], body.velocity[Y], body.velocity[Z]));
        //     }

        //     writer.close();
        //     outf.close();
        // } catch e {
        //     writeln(e);
        // }
    }
}
