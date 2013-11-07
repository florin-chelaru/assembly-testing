import math
import random

# Produces random assemblies, and reads with duplicates in certain regions of
# those assemblies. Our error detection should then detect duplicates in those
# regions of the assembly in which duplicates were introduced
class CorrectnessTests:

    def __init__(self):
        self.assembly = []
        self.reads = []
        self.readnamestart = 10000000

    # Generates contigcount number of contigs of an assembly, where each contig
    # is of length length
    # returns: an array of contigs, each a string of a, t, c, and g's.
    def gen_assembly(self, contigcount, length):
        self.assembly = []
        for s in range(contigcount):
            self.assembly.append("".join([random.choice(['a', 't', 'c', 'g']) for i in range(length)]))

        return self.assembly

    # Generates reads corresponding to the generated assembly.
    # readlen: the length of reads to produce
    #   note - some reads will be shorter if contig length is not a multiple of readlen
    # repeatcoords: list of tuples designating repeat locations, (contig, startpos)
    # returns: an array of strings of reads
    def gen_reads(self, readlen, repeatcoords):
        self.reads = []
        for contig in self.assembly:
            readcount = int(math.ceil(1.0 * len(contig) / readlen))
            for pos in [i * readlen for i in range(readcount)]:
                self.reads.append(contig[pos:min(pos + readlen, len(contig))])

        # insert repeats
        for repeat in repeatcoords:
            contig = self.assembly[repeat[0]]
            start = repeat[1]
            stop = min(start + readlen, len(contig))
            self.reads.append(contig[start:stop])

        return self.reads

    # Writes the last assembly generated to the given file
    def write_assembly(self, filename):
        f = open(filename, 'w')
        for i in range(len(self.assembly)):
            f.write(">" + str(i + 1) + "\n")
            f.write(self.assembly[i] + "\n")
        f.close()

    # Writes the last set of reads generated to the given file
    def write_reads(self, filename):
        f = open(filename, 'w')
        for i in range(len(self.reads)):
            f.write(">" + str(i + self.readnamestart) + "\n")
            f.write(self.reads[i] + "\n")
        f.close()

if __name__ == "__main__":
    tests = CorrectnessTests()
    tests.gen_assembly(4, 500)
    tests.gen_reads(100, [(0, 20), (3, 400)])
    tests.write_assembly("resources/assembly")
    tests.write_reads("resources/reads")
