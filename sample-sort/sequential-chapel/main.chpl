use Sort;
use Time;
use Random;

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

param REPEAT = 100;
param MY_RAND_MAX = (1 << 31) - 1;

proc getRandomNumber(seed: int) {
    return (seed * 1103515245 + 12345) & MY_RAND_MAX;
}

proc initRandomArray(n: int, initialSeed: int) {
    var arr: [{0..#n}] int;
    var randomNum = getRandomNumber(initialSeed);
    for i in 0..#n {
        arr[i] = randomNum;
        randomNum = getRandomNumber(randomNum);
    }
    return arr;
}

proc main() {
    var times: [{0..#REPEAT}] real;

    for i in 0..#REPEAT {
        var arr = initRandomArray(nkeys, i + 1);
        var watch: Timer;
        watch.start();
        myqsort(arr);
        times[i] = watch.elapsed();
    }

    writeln((+ reduce times) / REPEAT);
}

