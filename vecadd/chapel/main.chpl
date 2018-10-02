use Random;
use Time;

const seed = 42;
config var n = 10;

var A: [1..n] real;
var B: [1..n] real;
var alpha = 0.9;
fillRandom(A, seed);
fillRandom(B, seed);

const startTime = getCurrentTime();
var C = A + alpha * B;
const stopTime = getCurrentTime();
writeln("Elapsed time was: ", stopTime - startTime);
