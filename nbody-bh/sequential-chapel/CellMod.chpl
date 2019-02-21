module CellMod {
    use Common;
    use BodyMod;

    const N_CELL_CHILDREN = 8;
    type UCell = unmanaged Cell;

    class Cell {
        var children: [0..#N_CELL_CHILDREN] unmanaged Cell;
        var minBounds: Vec3;
        var maxBounds: Vec3;
        var body: UBody = nil;
        var mass: real;
        var cm: Vec3;

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

    proc freeCell(cell: UCell) {
        if cell == nil {
            return;
        }

        for child in cell.children {
            freeCell(child);
        }

        delete cell;
    }

    proc insertBody(cell: UCell, body: UBody) {
        if cell.body == nil && cell.isExternal() {
            cell.body = body;
            cell.mass = body.mass;
            cell.cm = body.position;
        } else {
            if cell.isExternal() {
                var halfSides: Vec3 = (cell.maxBounds - cell.minBounds) / 2;
                for i in 0..#N_CELL_CHILDREN {
                    var child: UCell = new UCell();
                    cell.children[i] = child;

                    var shifts = ((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1);
                    child.minBounds = cell.minBounds + (shifts * halfSides);
                    child.maxBounds = cell.maxBounds - ((1 - shifts) * halfSides);
                    if cell.body != nil && child.containsPosition(cell.body.position) {
                        insertBody(child, cell.body);
                        cell.body = nil;
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

    proc addForce(body: UBody, position: Vec3, mass: real) {
        var dx = position[X] - body.position[X];
        var dy = position[Y] - body.position[Y];
        var dz = position[Z] - body.position[Z];

        var dist = sqrt(dx*dx + dy*dy + dz*dz);
        var dist_cubed = dist * dist * dist;
        var f = (-G * body.mass * mass) / dist_cubed;

        body.force = body.force + (f * dx, f * dy, f * dz);
    }

    proc computeForce(cell: UCell, body: UBody) {
        if cell.mass == 0.0 {
           return;
        }
        if cell.body != nil {
            if cell.body.id == body.id {
                return;
            }
        }

        if cell.isExternal() {
            addForce(body, cell.body.position, cell.body.mass);
        } else {
            var diff = cell.maxBounds - cell.minBounds;
            var size = diff[X] + diff[Y] + diff[Z];

            var d = distance(cell.cm, body.position);
            if (size / d < THETA) {
                addForce(body, cell.cm, cell.mass);
            } else {
                for child in cell.children {
                    computeForce(child, body);
                }
            }
        }
    }
}
