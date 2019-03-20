module CellMod {
    use Common;
    use BodyMod;

    const N_CELL_CHILDREN = 8;
    type UCell = unmanaged Cell;

    class Cell {
        var children: [0..#N_CELL_CHILDREN] unmanaged Cell;
        var minBounds: Vec3;
        var maxBounds: Vec3;
        var body: Body = nil;
        var bodyPresent: bool = false;
        var mass: real;
        var cm: Vec3;
        var parentIdx = -1;
        var arrayIdx = -1;

        proc isExternal() {
            return children[0] == nil;
        }

        proc containsPosition(position: Vec3) {
            if position[X] < this.minBounds[X] || position[X] > this.maxBounds[X] {
                return false;
            }
            if position[Y] < this.minBounds[Y] || position[Y] > this.maxBounds[Y] {
                return false;
            }
            if position[Z] < this.minBounds[Z] || position[Z] > this.maxBounds[Z] {
                return false;
            }
            return true;
        }
    }

    type CellTuple = (int, Vec3, Vec3, Vec3, real);

    proc freeCell(cell: UCell) {
        if cell == nil {
            return;
        }

        for child in cell.children {
            freeCell(child);
        }

        delete cell;
    }
    
    proc cellContainsBounds(cell: UCell, minBounds: Vec3, maxBounds: Vec3) {
        for c in 1..DIMS {
            if cell.minBounds[c] > minBounds[c] || cell.maxBounds[c] < maxBounds[c] {
                return false;
            }
        }
        return true;
    }
    
    proc getCellsToSend(
        cell: UCell,
        parent: UCell,
        minBounds: Vec3,
        maxBounds: Vec3,
        minDepth: int,
        depth: int,
        cellsToSend: [] UCell
    ) {
        if depth > minDepth {
            cell.parentIdx = if depth - minDepth == 1 then -1 else parent.arrayIdx;
            cell.arrayIdx = cellsToSend.size;
            cellsToSend.push_back(cell);
        }

        var size = + reduce (cell.maxBounds - cell.minBounds);
        var boundsCenter = (maxBounds - minBounds) / 2.0;
        var d = distance(cell.cm, boundsCenter);

        if size / d >= THETA {
            for i in 0..#N_CELL_CHILDREN {
                if cell.children[i] != nil && cell.children[i].mass != 0.0 {
                    getCellsToSend(cell.children[i], cell, minBounds, maxBounds, minDepth, depth + 1, cellsToSend);
                } else {
                    break;
                }
            }
        }
    }

    proc reconstructReceivedCells(recvCells: [] CellTuple) {
        var allCells: [0..-1] UCell;
        var rootCells: [0..-1] UCell;

        for (parentIdx, cm, minBounds, maxBounds, mass) in recvCells {
            var cell: UCell = new UCell();
            cell.cm = cm;
            cell.minBounds = minBounds;
            cell.maxBounds = maxBounds;
            cell.mass = mass;

            if parentIdx == -1 {
                rootCells.push_back(cell);
            } else {
                var parent = allCells[parentIdx];
                for i in 0..#N_CELL_CHILDREN {
                    if cell.children[i] == nil {
                        parent.children[i] = cell;
                        break;
                    }
                }
            }
            allCells.push_back(cell);
        }
        return rootCells;
    }

    proc packCells(cells: [] UCell) {
        var packed: [0..-1] CellTuple;
        for cell in cells {
            packed.push_back((cell.parentIdx, cell.cm, cell.minBounds, cell.maxBounds, cell.mass));
        }
        return packed;
    }

    proc insertEmptyCell(cell: UCell, minBounds: Vec3, maxBounds: Vec3) {
        for i in 0..#N_CELL_CHILDREN {
            if cell.children[i] != nil && cellContainsBounds(cell.children[i], minBounds, maxBounds) {
                insertEmptyCell(cell.children[i], minBounds, maxBounds);
                return;
            } else if (cell.children[i] == nil) {
                var child: UCell = new UCell();
                cell.children[i] = child;
                child.minBounds = minBounds;
                child.maxBounds = maxBounds;
                return;
            }
        }
    }

    proc insertCell(cell: UCell, insert: UCell) {
        if cell.mass != 0.0 || insert.mass != 0.0 {
            cell.cm = (cell.mass * cell.cm + insert.mass * insert.cm) / (cell.mass + insert.mass);
            cell.mass += insert.mass;
        }

        for i in 0..#N_CELL_CHILDREN {
            if cell.children[i] != nil && cellContainsBounds(cell.children[i], insert.minBounds, insert.maxBounds) {
                insertCell(cell.children[i], insert);
                return;
            } else if cell.children[i] == nil {
                cell.children[i] = insert;
                return;
            }
        }
    }

    proc insertBody(cell: UCell, body: Body) {
        if !cell.bodyPresent && cell.isExternal() {
            cell.body = body;
            cell.mass = body.mass;
            cell.cm = body.position;
            cell.bodyPresent = true;
        } else {
            if cell.isExternal() {
                var halfSides: Vec3 = (cell.maxBounds - cell.minBounds) / 2;
                for i in 0..#N_CELL_CHILDREN {
                    var child: UCell = new UCell();
                    cell.children[i] = child;

                    var shifts = ((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1);
                    child.minBounds = cell.minBounds + (shifts * halfSides);
                    child.maxBounds = cell.maxBounds - ((1 - shifts) * halfSides);
                    if cell.bodyPresent && child.containsPosition(cell.body.position) {
                        insertBody(child, cell.body);
                        cell.bodyPresent = false;
                    }
                }
            }

            for child in cell.children {
                if child.containsPosition(body.position) {
                    insertBody(child, body);
                    break;
                }
            }

            cell.cm = (cell.mass * cell.cm + body.mass * body.position) / (cell.mass + body.mass);
            cell.mass += body.mass;
        }
    }

    proc addForce(inout body: Body, position: Vec3, mass: real) {
        var dx = position[X] - body.position[X];
        var dy = position[Y] - body.position[Y];
        var dz = position[Z] - body.position[Z];

        var dist = sqrt(dx*dx + dy*dy + dz*dz);
        var dist_cubed = dist * dist * dist;
        var f = (-G * body.mass * mass) / dist_cubed;
        var tmp = body.force + (f * dx, f * dy, f * dz);
        body.force = tmp;
    }

    proc computeForce(cell: UCell, inout body: Body) {
        if cell.mass == 0.0 {
           return;
        }
        if cell.bodyPresent {
            if cell.body.id == body.id {
                return;
            }
        }

        if cell.isExternal() && cell.bodyPresent {
            addForce(body, cell.body.position, cell.body.mass);
        } else {
            var diff = cell.maxBounds - cell.minBounds;
            var size = diff[X] + diff[Y] + diff[Z];

            var d = distance(cell.cm, body.position);
            if (size / d < THETA) {
                addForce(body, cell.cm, cell.mass);
            } else {
                for child in cell.children {
                    if child != nil {
                        computeForce(child, body);
                    }
                }
            }
        }
    }
}
