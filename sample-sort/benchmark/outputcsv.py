import sys
import json
import numpy as np

langs = ['c', 'chapel', 'julia']
threads = [1, 2, 4, 8, 16]

def main():
    args = sys.argv[1:]
    print(args)
    results = json.load(open('results-5e7.json', encoding='utf-8'))[args[0]]
    if args[1] == 'seq':
        seq_res = results['sequential']
        print(','.join(langs))
        print(','.join([str(np.array(seq_res[lang]).min()) for lang in langs]))
    elif args[1] == 'sm':
        sm_res = results['sm']
        print(','.join(langs))
        for thread in threads:
            print(','.join([str(np.array(sm_res[lang][str(thread)]).min()) for lang in langs]))

if __name__ == "__main__":
    main()