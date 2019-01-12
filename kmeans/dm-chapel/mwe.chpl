var arr: [{0..#numLocales}] int;

coforall L in Locales with (+ reduce arr) {
    on L {
        arr = here.id;
    }
}

writeln(arr);
