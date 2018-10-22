module DS {
    record Stack {
        type el_type;
        var stack_domain: domain(1);
        var stack: [stack_domain] el_type;
        var idx = 0;

        proc init(type el_type, capacity: int) {
            this.el_type = el_type;
            this.stack_domain = {0..#capacity};
        }

        proc push(el: el_type) {
            this.stack[this.idx] = el;
            this.idx += 1;
        }

        proc pop() {
            if (this.is_empty()) {
                writeln("Cannot pop empty stack.");
            }
            this.idx -= 1;
            return this.stack[this.idx];
        }

        proc peek() {
            if (this.is_empty()) {
                writeln("Cannot pek empty stack.");
            }
            return this.stack[this.idx - 1];
        }

        proc is_empty() {
            return this.idx == 0;
        }
    }
}