use Sort;
use Time;
use Random;

config const filename = "../data/arr10.txt";
config const nkeys = 10;

proc partition(arr: [] int, left_: int, right_: int) {
    var left = left_;
    var right = right_;
    var pivot = arr[right];
    while (left < right) {
        while (arr[left] < pivot) {
            left += 1;
        }
        while (arr[right] > pivot) {
            right -= 1;
        }
        if (left <= right) {
            var temp = arr[left];
            arr[left] = arr[right];
            arr[right] = temp;
        }
    }
    return left; // pivot index
}

proc myqsort(arr: [] int, left: int, right: int) {
    if (left >= right) {
        return;
    }

    var pivotIndex = partition(arr, left, right);
    myqsort(arr, left, pivotIndex - 1);  // sort left of pivot
    myqsort(arr, pivotIndex, right);  // sort right of pivot
}

proc myqsort(arr: [] int) {
    myqsort(arr, 0, arr.size - 1);
}

proc main() {
    var arr: [{0..#nkeys}] int;

    var f = open(filename, iomode.r);
    var reader = f.reader();
    reader.read(arr);
    f.close();
    reader.close();

    var watch: Timer;
    watch.start();
    myqsort(arr);
    writeln(watch.elapsed());
}
