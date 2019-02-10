use Sort;
use Random;
use Time;
use BlockDist;
use ReplicatedDist;
use CommDiagnostics;

config const filename = "../data/arr10.txt";
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

proc getSampleKeys(sampleKeys: [] int, arr: [] int, m: int) {
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
    for (binIdx, key) in zip(index_, arr) {
        binIdx = binarySearch(sampleKeys, key);
    }
    return index_;
}

proc computeBins(arr: [] int, sampleKeys: [] int, m: int) {
    var bins: [{0..#m}] unmanaged BinArray;
    for i in bins.domain {
        bins[i] = new unmanaged BinArray();
    }

    var localSubdomain = arr.localSubdomain();
    var subarray = arr.localSlice(localSubdomain);
    var index_ = mapKeysToBins(subarray, sampleKeys, m);
    for (binIdx, key) in zip(index_, subarray) {
        bins[binIdx].push(key);
    }
    return bins;
}

proc main() {
    var keysSpace = {0..#nkeys};
    var keysDomain = keysSpace dmapped Block(boundingBox=keysSpace);
    var arr: [keysDomain] int;
    var nbins = numLocales;

    coforall L in Locales do on L {
        var f = open(filename, iomode.r);
        var reader = f.reader();
        var localSubdomain = arr.localSubdomain();
        var i = 0;
        while i < localSubdomain.low {
            reader.read(int);
            i += 1;
        }
        reader.read(arr[localSubdomain]);
        reader.close();
        f.close();
    }

    var binsLocaleView = {0..0, 0..#nbins};
    var binsLocales: [binsLocaleView] locale = reshape(Locales, binsLocaleView);

    var binsSpace = {0..#nbins, 0..#nbins};
    var binsDomain = binsSpace dmapped Block(boundingBox=binsSpace, targetLocales=binsLocales);
    var bins: [binsDomain] unmanaged BinArray;

    var watch: Timer;
    watch.start();

    var sampleKeysDomain = {0..#nbins-1} dmapped Replicated();
    var sampleKeys: [sampleKeysDomain] int;
    getSampleKeys(sampleKeys, arr, nbins);

    for L in Locales {
        sampleKeys.replicand(L) = sampleKeys;
    }

    coforall L in Locales do on L {
        var binsRow = computeBins(arr, sampleKeys, nbins);
        for i in binsRow.domain {
            bins[here.id, i] = binsRow[i];
        }
    }

    var sortedArray: [0..-1] int;
    var sortedSubarrays: [0..#numLocales] unmanaged BinArray;
    coforall L in Locales do on L {
        var subarray: [0..-1] int;
        for i in 0..#nbins {
            subarray.push_back(bins[i, here.id].arr);
        }
        myqsort(subarray);
        sortedSubarrays[here.id] = new unmanaged BinArray();
        sortedSubarrays[here.id].arr.push_back(subarray);
    }

    for sub in sortedSubarrays {
        sortedArray.push_back(sub.arr);
    }

    writeln(watch.elapsed());
    if !isSorted(sortedArray) {
        writeln("Array not sorted!");
    }
}
