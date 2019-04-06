use Sort;
use Random;
use Time;

config const nkeys = 10;

const OVERSAMPLING_FACTOR = 128;

/// QuickSort
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
///

class BinArray {
    // Weird init for empty array
    var arr: [0..-1] int;

    proc push(el: int) {
        arr.push_back(el);
    }
}

proc getSampleKeys(arr: [] int, m: int) {
    var sampleKeys: [0..#m-1] int;
    var sampledKeys: [0..#m*OVERSAMPLING_FACTOR] int;

    // var randStream = new owned RandomStream(real, 1234);
    for i in sampledKeys.domain {
        // var idx: int = (randStream.getNext() * nkeys):int;
        sampledKeys[i] = arr[i];
    }

    for i in sampleKeys.domain {
        sampleKeys[i] = sampledKeys[(i + 1) * OVERSAMPLING_FACTOR];
    }

    myqsort(sampleKeys);

    return sampleKeys;
}

proc binarySearch(arr: [] int, el: int) {
    var left = 0;
    var right = arr.size;

    while (left < right) {
        var middle = ((left + right) / 2):int;
        if (arr[middle] >= el) {
            right = middle;
        } else {
            left = middle + 1;            
        }
    }
    return left;
}

proc mapKeysToBins(arr: [] int, sampleKeys: [] int, m: int) {
    var index_: [0..#arr.size] int;
    for i in index_.domain {
        var binIdx = binarySearch(sampleKeys, arr[i]);
        index_[i] = binIdx;
    }
    return index_;
}

proc bin(arr: [] int, m: int, bins: [{0..#m}][{0..#m}] unmanaged BinArray) {
    var sampleKeys = getSampleKeys(arr, m);
    var blockSize = ((nkeys + m - 1) / m):int;

    forall threadId in 0..#m {
        var start = threadId * blockSize;
        var end = min(start + blockSize, nkeys);
        var subarraySize = end - start;
        var subarrayDomain = 0..#subarraySize;
        var subarray = arr[start..end-1].reindex(subarrayDomain);

        var index_ = mapKeysToBins(subarray, sampleKeys, m);
        for i in subarrayDomain {
            var binIdx = index_[i];
            bins[threadId][binIdx].push(arr[start + i]);
        }
    }

    return bins;
}

proc subsort(m: int, bins: [{0..#m}][{0..#m}] unmanaged BinArray) {
    var sortedArray: [0..#nkeys] int;
    var colSum: [0..#m] int;
    var prefixColSum: [0..#m] int;

    for i in 0..#m {
        for j in 0..#m {
            var binSize = bins[j][i].arr.size;
            colSum[i] += binSize;
        }
    }

    for i in 1..#m-1 {
        prefixColSum[i] = prefixColSum[i - 1] + colSum[i - 1];
    }

    forall threadId in 0..#m {
        var subarray: [0..#colSum[threadId]] int;

        var offset = 0;
        for i in 0..#m {
            var binArray = bins[i][threadId];
            var size = binArray.arr.size;
            subarray[offset..offset+size-1] = binArray.arr;
            offset += size;
        }
        myqsort(subarray);

        var start = prefixColSum[threadId];
        var end = start + colSum[threadId];
        sortedArray[start..end-1] = subarray;
    }
    return sortedArray;
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
        var nbins = here.maxTaskPar;
        var watch: Timer;
        watch.start();

        // Init bins
        var bins: [{0..#nbins}][{0..#nbins}] unmanaged BinArray;
        for i in 0..#nbins {
            for j in 0..#nbins {
                bins[i][j] = new unmanaged BinArray();
            }
        }

        bin(arr, nbins, bins);

        var sortedArray = subsort(nbins, bins);

        times[i] = watch.elapsed();

        if !isSorted(sortedArray) {
            writeln("Array not sorted!");
        }

        for i in 0..#nbins {
            for j in 0..#nbins {
                delete bins[i][j];
            }
        }
    }
    writeln((+ reduce times) / REPEAT);
}

