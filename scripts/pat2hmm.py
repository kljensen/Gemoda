#!/usr/bin/env python
import sys, string, re, os, tempfile
from optparse import OptionParser

def getCommandLine():
    usage = "usage: %prog -nN -iFILE -oOUTFILE\n\n" + \
            "Converts a pattern from Gemoda into a hidden markov\n" + \
            "model using the Hmmer package available from\n" + \
	    "http://hmmer.wustl.edu/\n"
    parser = OptionParser(usage=usage)

    # get the pattern number
    parser.add_option("-n", "--number", type="int", dest="patNo",
                        metavar="N", help="number of pattern you want",
                        default=-1)

    # get the input file
    parser.add_option("-i", "--infile", type="string", dest="inFileName",
                        metavar="FILE", help="name of gemoda pattern file")

    # get the output file
    parser.add_option("-o", "--outfile", type="string", dest="outFileName",
                        metavar="OUTFILE", help="name of output file")

    # parse the command line
    (options, args) = parser.parse_args()

    # do some validity checking
    if options.patNo == -1 or options.inFileName == None or options.outFileName == None:
        parser.print_help()
        sys.exit()

    return (options.patNo, options.inFileName, options.outFileName)

def openFiles(nameArray):

    fh = []
    for file in nameArray:
        try:
            fh.append(open(file, "r"))
        except IOError, (errno, strerror):
            print "Error %s: %s" % (errno, strerror)
            sys.exit()

    return fh

def writeAlnFile(fh, n):

    #  A regular expression object that indicates the
    #  start of a new pattern
    patRegex = re.compile("^pattern (?P<patNo>[0-9]+)")
    instanceRegex = re.compile("^\s+(?P<seqNo>[0-9]+)\s+(?P<posNo>[0-9]+)\s+(?P<string>[A-Z]+)")


    #  Find the pattern we want
    i=0
    for line in fh:
        m = patRegex.match(line)
        if m != None and int(m.group('patNo')) == n:
            break
        i+=1

    i=0
    alnOutput = ""
    for line in fh:
        m = instanceRegex.match(line)
        if m == None:
            break
        else:
            alnOutput = alnOutput + "seq%s_%s\t\t%s\n" % (m.group('seqNo'), m.group('posNo'), m.group('string'))
        i+=1

    alnOutput = string.expandtabs(alnOutput, 2)

    (fileno,alnFile) = tempfile.mkstemp(prefix=sys.argv[0])
    try:
        alnFileH = open(alnFile, "w")
    except IOError, (errno, strerror):
        print "Error %s: %s" % (errno, strerror)
        sys.exit()

    alnFileH.write(alnOutput)
    alnFileH.close()
    return (alnFile)

def makeHmm(alnFile, n, hmmFile):
    os.system("hmmbuild -F -n \"gemoda-pattern-no-%d\" %s %s" % (n, hmmFile, alnFile))


# Main
def main():

    # get stuff from the command line
    (n, fin, fout) = getCommandLine()

    # open the input and output files
    [fh] = openFiles([fin])

    # write a temporary aln file to pass
    # to hmmbuild
    alnFile = writeAlnFile(fh, n)

    fh.close()

    # get the profile from hmmbuild
    makeHmm(alnFile, n, fout)

    os.remove(alnFile)

    sys.stderr.write("\n\n\tMade hmm file: %s" % fout)
    sys.stderr.write("\n\n\tTo search with hmmsearch use:")
    sys.stderr.write("\n\thmmsearch %s <your FastA file>" % (fout))
    sys.stderr.write("\n\n\tTo make a multiple sequence alignment using hmmalign:")
    sys.stderr.write("\n\thmmalign %s <your FastA file>" % (fout))
    sys.stderr.write("\n\n\t(See the Hmmer documentation for help.)\n\n")

if __name__ == "__main__":
  sys.exit(main())


