#!/usr/bin/env python
import sys, getopt, string, re

def usage():
    print "%s: PDB_FILE [OTHER_PDB_FILES]" % (sys.argv[0])

def printFasta(fhout,label,coords):
    fhout.write("%s\n" % label )
    for row in coords:
        fhout.write("%s\n" % string.join(row, " ") )

def processPdbFile(fh, defaultName):
    headers = []
    atoms = []

    isHeader = re.compile("HEADER")
    isAtom = re.compile("ATOM")
    isCA = re.compile("CA")
    #isTitle = re.compile("TITLE")

    aaDict ={   'ALA':'A',
                'CYS':'C',
                'ASP':'D',
                'GLU':'E',
                'PHE':'F',
                'GLY':'G',
                'HIS':'H',
                'ILE':'I',
                'LYS':'K',
                'LEU':'L',
                'MET':'M',
                'ASN':'N',
                'PRO':'P',
                'GLN':'Q',
                'ARG':'R',
                'SER':'S',
                'THR':'T',
                'VAL':'V',
                'TRP':'W',
                'TYR':'Y'}


    id=""
    lastChain = ""
    seenChains = -1
    chainSizes = []
    chain = ""
    coords = [[], [], []]
    for line in fh:
        line = string.strip(line)

        if isHeader.match(line):
            id = string.split(line)[-1]
            headers.append(id)

        elif isAtom.match(line):
            atomFields = string.split(line)
            atomType = atomFields[2]
            if isCA.match(atomType) == None:
                continue

            chain = atomFields[4]

            if chain != lastChain:
                chainSizes.append(0)
                seenChains += 1
                lastChain = chain

                if seenChains is not 0:
                    if id is "": id = defaultName
                    printFasta(sys.stdout, ">%s-chain_%s" % (id, chain), coords)
                    del coords
                    coords = [[], [], []]


            residue = atomFields[3]
            if aaDict.has_key(residue):
                residue =  aaDict[residue]
            else:
                residue = 'X'

            [x, y, z] = atomFields[6:9]

            coords[0].append(x)
            coords[1].append(y)
            coords[2].append(z)

            chainSizes[seenChains-1] += 1

    printFasta(sys.stdout, ">%s-chain_%s" % (id, chain), coords)


def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(2)


    i=0
    for pdbFile in sys.argv[1:]:

        try:
            fh = open(pdbFile, "r")
        except IOError, (errno, strerror):
            print "Error %s: %s" % (errno, strerror)
            sys.exit()

        processPdbFile(fh, "file%d" % (i))
        i+=1



if __name__ == "__main__":
    sys.exit(main())
