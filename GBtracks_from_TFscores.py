#!/usr/bin/python

import os, sys, string, argparse

#===============================================================================
#  Command Line Arguments

parser = argparse.ArgumentParser(description = 'Generating Genome Browser Tracks')
parser.add_argument('-t', metavar = 'modeltype', 
                    help = 'Type of scores, either PWM (log ratio), or SVR (normalized between 0 and 1)', 
                    dest = 'runtype' ,
                    choices = ['SVR','PWM'] ,
                    required=True)
parser.add_argument('-s', metavar = 'SequenceFile', 
                    help = 'Sequence file where the first three columns are "Chromosome Name", "Start Position", and "Stop Position" (i.e. .bed format), and may or may not also have the sequences', 
                    dest = 'seqfile' , 
                    required=True)
parser.add_argument('-o', metavar = 'OutFile', 
                    help = 'Optional, the name of the output file', 
                    dest = 'outfile')
args = parser.parse_args()
seqfile = args.seqfile
runtype = args.runtype
if args.outfile: outfile = args.outfile
else: outfile = os.path.splitext(seqfile)[0] + '_browser-track.bed'


def read_data(filename):
    ''' Creates an array for every cell in a tab separated text file'''
    data = []
    if os.path.isfile(filename): f = open(filename, 'r') #opens the file as "f"
    else:
        sys.exit("The file", filename, "does not exist!")
    for line in f: #for each line in the file
        l = string.split(line.strip(), '\t') #removes any carriage returns, then splits the tab separated line into columns
        data.append(l) #Add each line to the array
    f.close() #closing the file
    return data



def get_genome_browser_tracks(fullscoresfile):
    '''Takes a file generated from a "peak_full_scores..." module, and generates a custom trac for the genome browser in np format.
    Columns are chromosome, start position, end position, name, score.'''
    data = read_data(fullscoresfile)
    
    # Getting these from the name of the file
    if args.runtype == 'SVR': info = 'E2F1 SVR Model'
    elif args.runtype == 'PWM': info = 'E2F Transfac PWM'

    for peak in data:
        numscores = len(peak)-4 #first 4 are not scores, but everything else is
        chrom,start,stop,seq=peak[:4]
        f = open(outfile, 'w')
        print "writing results to", outfile
        print >>f, 'browser position '+chrom+':'+start+'-'+stop
        print >>f, '#browser hide all\nbrowser pack refGene encodeRegions\nbrowser full altGraph'
        print >>f, 'track type=bedGraph name="'+outfile+'" description="'+info+'" visibility=full color=200,100,0 altColor=0,100,200 priority=20'
        start=int(start)
        stop=int(stop)
        for i in range(numscores):
            #if float(peak[4+i]) > 0: print chrom, int(start)+i, int(start)+i+36, peak[4+i] # for testing
            if args.runtype == 'SVR':
                newstart = (int(start)+i)+17
                if float(peak[4+i]) > 0.207: print >>f, chrom, newstart, newstart+4,i, float(peak[4+i])*100 #, seq[i:i+36]
            elif args.runtype == 'PWM':
                newstart = (int(start)+i)+4
                if float(peak[4+i]) > 3: print >>f, chrom, newstart, newstart+1,i, float(peak[4+i])*10 #, seq[i:i+36]
        f.close()
        
get_genome_browser_tracks(seqfile)

