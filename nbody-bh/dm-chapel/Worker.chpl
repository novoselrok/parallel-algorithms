module Worker {
    use Math;
    use Sort;
    use IO;

    use Common;
    use BodyMod;
    use CellMod;
    use DArrayMod;

    type Bounds = (Vec3, Vec3);
    type PartnerInfo = (int, bool);

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

    proc determineGroupPartners(group: int, groupPartners$: [] sync int) {
        var myGroupPartners: [0..-1] int;
        var rank = here.id;

        myGroupPartners.push_back(rank);
        for i in 0..#groupPartners$.size-1 {
            groupPartners$[rank] = group;

            var readFrom = (rank + i + 1) % numLocales;
            var groupPartner = groupPartners$[readFrom];
            if readFrom != rank && groupPartner == group {
                myGroupPartners.push_back(readFrom);
            }
        }

        sort(myGroupPartners);
        return myGroupPartners;
    }

    proc fracWeightBelow(
        bodies: [] Body,
        split: real,
        coord: int,
        workBelow$: [] sync real,
        workAbove$: [] sync real,
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
            workBelow$[rank] = workBelow;
            workAbove$[rank] = workAbove;

            var readFromIdx = (idxWithinGroup + i + 1) % groupPartners.size;
            var readFrom = groupPartners[readFromIdx];
            if readFrom != rank {
                allWorkAbove += workAbove$[readFrom];
                allWorkBelow += workBelow$[readFrom];
            }
        }

        var fraction = allWorkBelow / (allWorkBelow + allWorkAbove);
        return fraction - 0.5;
    }

    proc bisection(
        _min: real,
        _max: real,
        bodies: [] Body,
        coord: int,
        workBelow$: [] sync real,
        workAbove$: [] sync real,
        groupPartners: [] int
    ) {
        var min = _min;
        var max = _max;
        var fmin = fracWeightBelow(bodies, min, coord, workBelow$, workAbove$, groupPartners);
        var mid = 0.0;
        // writeln("fmin: ", fmin);
        // try! stdout.flush();
        var iteration = 0;
        while iteration < BISECTION_MAX_ITER && abs((max - min) / 2) > BISECTION_TOL {
            // writeln("iteration: ", iteration);
            // try! stdout.flush();
            mid = (min + max) / 2;
            var fmid = fracWeightBelow(bodies, mid, coord, workBelow$, workAbove$, groupPartners);
            // writeln("fmid: ", fmid);
            // try! stdout.flush();
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
        workBelow$: [] sync real,
        workAbove$: [] sync real,
        otherBodiesWrappers: [] owned BodiesArray,
        otherBodiesWrappers$: [] sync bool,
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
        // writeln("nsplits ", nSplits);
        // try! stdout.flush();
        for i in 0..#nSplits {
            var nProcsLeft = (numLocales / (2**i)):int;

            group = group << 1;
            if aboveSplit {
                group = group | 1;
            }

            var groupPartners = determineGroupPartners(group, groups$);

            // writeln(rank, " group partners ", groupPartners);
            // try! stdout.flush();

            var coord = (i % DIMS) + 1;

            var split = bisection(myMin[coord], myMax[coord], bodies, coord, workBelow$, workAbove$, groupPartners);

            // writeln(split);
            // try! stdout.flush();
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

            var partnerRank = getPartnerRank(rank, nProcsLeft);
            writeln(rank, " partnerrank ", partnerRank);
            try! stdout.flush();
            partners.push_back((partnerRank, aboveSplit));

            var otherBodiesWrapper = new owned BodiesArray();
            otherBodiesWrapper.elements.push_back(otherBodies);

            otherBodiesWrappers[rank] = otherBodiesWrapper;
            otherBodiesWrappers$[rank] = true;

            // Wait until its full
            otherBodiesWrappers$[partnerRank];
            // Read the bodies
            bodies.clear();
            bodies.push_back(myBodies);
            var partnersOtherBodies = otherBodiesWrappers[partnerRank];
            bodies.push_back(partnersOtherBodies.elements);
        }
        writeln(rank, " ", bodies.size);
        return (bodies, myBounds, otherBounds, partners);
    }

    proc work(
        inputBodies: [] Body, 
        bounds: [] Bounds,
        bounds$: [] sync bool,
        workBelow$: [] sync real,
        workAbove$: [] sync real,
        otherBodies: [] owned BodiesArray,
        otherBodies$: [] sync bool,
        groups$: [] sync int,
        iterations: int
    ) {
        var bodies: [0..-1] Body;
        var rank = here.id;

        for inputBody in inputBodies {
            bodies.push_back(inputBody);
        }

        for iteration in 0..#iterations {
            writeln(iteration);
            try! stdout.flush();
            var myUniverseBounds = getUniverseSize(bodies);

            var universeMin = myUniverseBounds[1];
            var universeMax = myUniverseBounds[2];
            for i in 0..#bounds.size-1 {
                bounds[rank] = myUniverseBounds;
                bounds$[rank] = true;

                var readFrom = (rank + i + 1) % numLocales;
                if (readFrom != rank) {
                    // Blocks until readFrom locale writes to it
                    bounds$[readFrom];
                    var (minOtherBounds, maxOtherBounds) = bounds[readFrom];
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
                }
            }
            // writeln(rank, " universe ", universeMin, " ", universeMax);
            // try! stdout.flush();

            var (orbBodies, myBounds, otherBounds, partners) = orb(
                bodies,
                universeMin,
                universeMax,
                workBelow$,
                workAbove$,
                otherBodies,
                otherBodies$,
                groups$
            );
            bodies.clear();
            bodies.push_back(orbBodies);
            writeln(rank, " ", bodies.size);
        }
        // writeln("length: ", bodies.size);
    }
}
