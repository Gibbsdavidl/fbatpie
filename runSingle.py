
#!/usr/bin/python

import sys, getopt
from src import fbat

def main(argv):
    tpedfile = ''
    tfamfile = ''
    offset=0
    freq=0
    verbose="not"
    try:
        opts, args = getopt.getopt(argv,"h",["tped=","tfam=","offset=","freq="])
    except getopt.GetoptError:
        print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset <o> --freq <cutoff> --verbose v'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset <o> --freq <cutoff> --verbose v'
            sys.exit()
        elif opt in ("--tped"):
            tpedfile = arg
        elif opt in ("--tfam"):
            tfamfile = arg
        elif opt in ("--offset"):
            offset = arg
        elif opt in ("--freq"):
            freq = arg
    if tfamfile == '' or tpedfile == '':
        print 'runSingle.py --tped <tped file> --tfam <tfam file> --offset o --freq f --verbose v'
        sys.exit()
            
    f = fbat.FBAT()
    f.load(tfamfile, tpedfile)
    f.setOffset(offset)
    f.setVerbose(verbose)
    f.setFreqCutoff(freq)
    f.single()

if __name__ == "__main__":
   main(sys.argv[1:])
   
