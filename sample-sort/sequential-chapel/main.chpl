use Sort;
use Time;

config const filename = "../data/arr10.txt";
config const nkeys = 10;

proc main() {
    var arr: [{0..#nkeys}] int;

    var f = open(filename, iomode.r);
    var reader = f.reader();
    reader.read(arr);
    f.close();
    reader.close();

    var watch: Timer;
    watch.start();
    sort(arr);
    writeln(watch.elapsed());
}
