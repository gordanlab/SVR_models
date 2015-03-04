#!/usr/bin/python

import argparse, string, sys, os, time, itertools, re, operator, random
from subprocess import *

print '''\nRunning the program to get scores for genomic sequences using a SVR or PWM model.
Check example files for proper input file formats.  
Generally speaking, input files should not contain any headers.

Note!  If the program needs to get genomic sequences, more than 3GB of memory will be required.'''


#===============================================================================
#  Command Line Arguments

parser = argparse.ArgumentParser(description = 'E2F Binding Model')
parser.add_argument('-t', metavar = 'RunType', 
                    help = 'Type of run, either using PWM, or SVR for finding scores', 
                    dest = 'runtype' ,
                    choices = ['SVR','PWM'] ,
                    required=True)
parser.add_argument('-g', metavar = 'GenomeFile', 
                    help = 'Genome File - Fasta format, containing the sequence for each chromosome as a separate entry', 
                    dest = 'genomefile')
parser.add_argument('-s', metavar = 'SequenceFile', 
                    help = 'Sequence file where the first three columns are "Chromosome Name", "Start Position", and "Stop Position" (i.e. .bed format), and may or may not also have the sequences', 
                    dest = 'seqfile' , 
                    required=True)
parser.add_argument('-n', metavar = 'NegativeSequenceFile', 
                    help = "Sequence file for negative control sequences for ROC curves. If this flag is used, ROC data and AUC will be calculated. If this file contains more sequences than the input sequence file, it will be normalized. This file could be, for example, DNase hypersensitive sites that don't overlap the ChIP peaks given as the input sequence file" , 
                    dest = 'negseqfile')
parser.add_argument('-m', metavar = 'ModelFile', 
                    help = 'The E2F .model file generated from LibSVM, or the PWM model file (in probabilities or log ratios)' , 
                    dest = 'modelfile' , 
                    required=True)
parser.add_argument('-o', metavar = 'OutFilePrefix', 
                    help = 'Optional, the prefix that all output files will be based on (do not include file extension)', 
                    dest = 'outprefix')
parser.add_argument('--direction',
                    help =  'Optional, specify the direction of the sequences to get scores for; "fwd" for scoring the forward (sense) strand, "rev" for scoring the reverse (anti-sense) strand, and "best" (default) to use the best score from both directions for each position.' ,
                    choices = ['fwd','rev','best'] )

args = parser.parse_args()
genomefile = args.genomefile
seqfile = args.seqfile
modelfile = args.modelfile
runtype = args.runtype
if args.outprefix: outprefix = args.outprefix
else: outprefix = os.path.splitext(seqfile)[0] + '_' + modelfile + '_' + runtype + 'predict'
if args.negseqfile: 
    negseqfile = args.negseqfile
    print "\nNegative sequence file has been specified, so ROC curve and AUC will be calculated."

    
#===============================================================================


'''Preparation Stuff============================================================'''



infofilename = outprefix + '_info.txt'
f_info = open(infofilename, 'a')
print "\nWriting info about this run to", infofilename
print >>f_info, '\n================================================================================\n' , \
'\nRun started at', time.strftime("%m\%d\%Y %H:%M") ,\
'\nRun type is', runtype ,\
'\nOutput files will be given the prefix:', outprefix,\
'\nSequence file being used:', seqfile,\
'\nModel file being used (SVR or PWM file):', modelfile
if args.genomefile:
    print >>f_info, '\nGenome file (if needed)', genomefile
if args.negseqfile: 
    negseqfile = args.negseqfile
    print >>f_info, '\nNegative sequence file for ROC:', negseqfile 
    

''' Parameters =========================================================='''

### These parameters are querried by the program

### Files with large numbers of sequences need to be broken up into chunks for LIBSVM to run efficently
chunksize = 1000

if runtype == 'SVR':
    ### Length of the sequences used for building the model: 36 is the default value
    length = 36 
    
    ### Defining what kind of features we want. I.e. '1' for 1mers, or 123. Best value for E2Fs is '3'.
    rawkmer = '3'
    kmers = []
    for item in str(rawkmer): kmers.append(int(item))
    kinfo = ''
    for x in kmers: kinfo = kinfo + str(x) + " + " 
    kinfo = kinfo[:-3] + " mer features"

    ### Sequences that don't match the criteria (SVR) should be assigned a low score 
    badscore = 0

### Defining whether we're scoring the forward sequences ('fwd'), reverse ('rev'), or the best of both sequences ('best')
if args.direction: seqtype = args.direction
else: seqtype = 'best'
    

'''Defining the Modules ========================================================='''
    
def read_data(filename):
    ''' Creates an array for every cell in a tab separated text file'''
    data = []
    try: f = open(filename, 'r') #opens the file as "f"
    except IOError: 
        print "Could not open the file:", filename
        sys.exit()
    for line in f: #for each line in the file
        l = string.split(line.strip(), '\t') #removes any carriage returns, then splits the tab separated line into columns
        data.append(l) #Add each line to the array
    f.close() #closing the file
    return data

def read_pwm(pwmfile):
    '''function that reads a pwm from a text file and converts the information into a table, (uses the top one if multiple).
    The first column contains the corresponding nucleotide (A, C, G, T) - the other columns contain probabilities in the 
    float format. '''
    import math, copy
    
    pwmdata_1 = read_data(pwmfile)
    
    ### Removing all the blank lines in list (or, rather, blank lists in the array), and getting only the first PWM if multiples
    pwmdata_2 = [ x for x in pwmdata_1 if x != ['']]
    #for line in pwmdata_2: print line
    
    pwmdata = []
    counter = 0
    for i in range(len(pwmdata_2)):
        if counter < 4 and any(string in str(pwmdata_2[i][0]) for string in ['A','T','C','G']):
            pwmdata.append(pwmdata_2[i])
            counter += 1
    
    ### Converting all the PWM values to float
    for x in range(4): #for every line with probabilities for the bases
        for y in range(1,len(pwmdata[1])):
            pwmdata[x][y] = float(pwmdata[x][y])
    
    ### Finding the sum of the columns, so if it's not 1, we can convert the numbers to probabilities.
    coltotals = []
    for n1 in range(1,len(pwmdata[1])): #for every column
        coltotal = 0
        for n2 in range(4): #for every line with probabilities for the bases
            coltotal += pwmdata[n2][n1]
        coltotals.append(coltotal) 
    if any(x > 1 for x in coltotals): #if any of the values in the list are greater than 1
        for n1 in range(1,len(pwmdata[1])): #for every column
            for n2 in range(4): #for every line with probabilities for the bases
                pwmdata[n2][n1] = pwmdata[n2][n1]/(coltotals[n1-1])
    
    ## Creating a new list, where we can convert the values to a probability matrix
    #print 'Converting values to a log probability matrix'
    pwm_logs = copy.deepcopy(pwmdata) 
    for x in range(4): #for every line with probabilities for the bases
        for y in range(1,len(pwm_logs[1])):
            if pwm_logs[x][y] == 0: pwm_logs[x][y] = 0 # continue
            else: pwm_logs[x][y] = math.log(pwm_logs[x][y]/0.25) #converting to log probability adjusted for base frequence
    
    pwm = pwm_logs #for if we're not trimming the PWM, otherwise, comment this out and uncomment everything below
    ### Trimming the PWM, if needed
    #print 'Trimming the PWM, if needed'
    ### Getting the info content
#     infocontent = []
#     for n1 in range(1,len(pwmdata[1])):  #for each column with values (first has the base)
#         #Getting the information content for each column: 2+sum(p*log2 p)
#         info = 2
#         for n2 in range(4): #n2 is the value for each column (the line), skipping the first (header)
#             info = info + (pwmdata[n2][n1] * math.log(pwmdata[n2][n1],2)) 
#         infocontent.append(info)
#     ### Doing the trimming
#     pwm = copy.deepcopy(pwm_logs)
#     for n1 in range( len(infocontent) ):  #n1 is the column for the pwm, which corresponds to n1-1 in infocontent
#         '''Checks to make sure this column AND the following have info contents of 
#             greater than 0.3 (need 2 in a row), for the Left side'''
#         if n1 < (len(infocontent) / 2 ): #working on the first half of the sequence (so trimming from the right)
#             if not ( infocontent[n1] > 0.3 and infocontent[n1+1] > 0.3 ):
#                 for n2 in range(4): # for each value in the column (each line)
#                     #pwm values in columns with information content less than 0.3, are replaced with blanks
#                     pwm[n2][n1+1] = '' #Need to use n1+1 because pwm file has headder column
#         '''Checks to make sure this column AND the preceeding have info contents of 
#             greater than 0.3 (need 2 in a row), for the Right side'''
#         if n1 >(len(infocontent) / 2 ): #working on the last half of the sequence (so trimming from the left)
#             if not ( infocontent[n1] > 0.3 and infocontent[n1-1] > 0.3 ):
#                 for n2 in range(4):
#                     pwm[n2][n1+1] = '' #Need to use n1+1 because pwm file has headder column
#     for n3 in range(4): #for each value in the column...
#         pwm[n3] = filter(None, pwm[n3]) #filter out any values that are blank
    
    return pwm 

def reverse_complement(seq):
    '''Takes any sequence, and gets the reverse complement.  Note, only changes A, T, C, or G. Anything else will be left as is.'''
    comp = string.maketrans("ATCG", "TAGC") #creates the translation definitions, required by line below
    rev_comp = (seq[::-1]).translate(comp) #First reverses the sequence, then translates it according to rules above
    return rev_comp

def load_genome(genome_file):  #created by Alina
    '''load the whole genome: read from file and then construct a dictionary with a 
    loooong sequence for each chromosome.'''
    if not args.genomefile:
        sys.exit("A genome file must be specified using the -g flag")
    if os.path.isfile(genome_file): f = open(genome_file,'r')
    else:
        sys.exit("The file", genome_file, "does not exist! Restart the program and specify a genome sequence file with the -g flag")
    try:
        memory = int(os.popen("free -m").readlines()[1].split()[1])
        if memory < 3000: print "Warning, This program needs at least 3GB of free memory for loading the genome file, and it looks like you have less than that. Trying to continue..."
    except: 
        print "This program needs at least 3GB of free memory for loading the genome file"
        pass
    print 'Loading the whole genome from', genome_file
    genome = {}
    
    curr_chr = 'chr1'
    curr_seq = ''
    for line in f:
        # removing the '\n' at the end
        line = line.strip()
        # verifying if it is a chromosome name line
        if line[0]=='>':
            genome[curr_chr] = curr_seq # append the previous chromosome to the genome
            #print curr_chr+' - '+str(len(curr_seq))
            # if the next chromosome is a variation (ex: 'chr17_ctg5_hap1'), then stop reading
            if len(line)>6:
                break
            curr_chr = line[1:] #the name of the new chromosome
            curr_seq = ''
            print 'Reading '+curr_chr+' ...'
        else:
            curr_seq +=line
    genome[curr_chr] = curr_seq #adding the last chromosome to the dictionary
    f.close()
#    for line in genome: print line
    return genome

def get_peak_lengths(npfile):
    '''Takes a narrow peaks file, gets the length of each peak, and writes this to a new column, in a new file'''
    peaks = read_data(npfile)
    outdata = []
    #outfile = os.path.splitext(npfile)[0]+'_peaklengths.pk'
    #f=open(outfile, 'w')
    for line in peaks:
        seq,beg,end = line[0],int(line[1]),int(line[2])
        total = end - beg
        if total > 1000: total = 1001
        newline = [seq,beg,end,total]
        #print >>f, ("\t".join(map(str,newline)))
        outdata.append(newline)
    #f.close()
    return outdata

def normalize_peak_lengths(peak1file,peak2file):
    '''Takes two narrow peak or bed files, and reduces the size of the larger one to the size of the smaller file
    by taking each sequence in the smaller file, finding an equal or larger sequence in the larger file, then
    reduces the new sequence to the length of the smaller (if not equal), and writes this to a new list, and
    also deletes the peak from the larger file so it won't be re-used. Returns the complete smaller set, and 
    the shrunk larger set. Data must be in the format with col1 = chromosome, col2 = peakstart, col3 = peakend'''
    
    ### Reading the np files and setting the lengths of the peaks as the 4th column, and getting rid of the 
    ###  rest of the fluff in the np files
    peak1data = get_peak_lengths(peak1file)
    peak2data = get_peak_lengths(peak2file)
    
    ### setting the larger file to set1, and the smaller to set2
    if len(peak1data) > len(peak2data):
        set1 = peak1data
        set2 = peak2data
        outfile = os.path.splitext(peak1file)[0]+'_normalized-to_'+os.path.splitext(peak2file)[0]+'.pk'
    else:
        set1 = peak2data # we need to randomize this
        set2 = peak1data
        outfile = os.path.splitext(peak2file)[0]+'_normalized-to_'+os.path.splitext(peak1file)[0]+'.pk'
    
    if os.path.isfile(outfile):
        print "Normalized file already exists.  Using this file."
        return outfile, read_data(outfile)
    
    random.shuffle(set1)
    print "larger file is", len(set1), "smaller is", len(set2)
    
    ### Creating a new set of peaks, from the larger set, that has the same number of peaks, with 
    ###  identical lengths to those in the smaller file
    newdata = []
    for peak in set2: #for every peak in the smaller set
        size = peak[3] #getting the length of the peak
        #print "\npeak size is", size
        for x in range(len(set1)): #find a peak in the randomized larger set 
            peak2 = set1[x]
            size2 = peak2[3] #getting the length of the larger peak
            #print "\ncomparing", size, "to", size2
            if size2 >= size: #only if the length of the new peak is more than the starting peak
                #print "peak 2 is larger!"
                start,stop = peak2[1],peak2[2] #getting the start and stop of the new peak
                newstart = random.randint(start,(stop-size))
                newstop = newstart + size
                newpeak = [peak2[0],newstart,newstop,(newstop-newstart)]
                newdata.append(newpeak)
                del set1[x]
                #print "old", start, stop, (stop-start), "- new", newstart, newstop, (newstop-newstart)
                #print peak,peak2,newpeak,"\n"
                break
    
    print "writing data to", outfile, "newsize is", len(newdata) 
    f=open(outfile,"w")
    for line in newdata: print >>f, ("\t".join(map(str,line)))
    f.close()

    return outfile, newdata #output file name, then the actual data in the file

def binding_prediction_chipall(modelfile,seqfile,genomefile):
    '''Takes a ChIP peak file, uses a libsvm regression model file and libsvm to predict the binding intensity
    for every kmer (length according to the model) in each peak, returns all scores.'''
    dataseqfile = os.path.splitext(seqfile)[0]+'_sequences.txt' #getting the file name
    
    ### Checking if this file already contains sequences, so we don't have to bother getting them
    seqdata = read_data(seqfile)
    linetest = seqdata[0] #getting the first line from the file for testing
    seqloc = 0
    
    if len(linetest) > 3:
        for col in linetest[3:]: #for each column, starting with the 4th (cant be in first three)
            if not re.search(r'[^ATGC]', col) and len(col) == int(linetest[2]) - int(linetest[1]): # if this column contains only A, T, G, or C, and is the right length
                    seqloc = linetest.index(col) #getting the column number
    if seqloc > 0: # if this exists, then we can use the sequence listed in this file
        print "\nInput file already contains sequences, using these."
        if seqloc != len(linetest)-1: #if sequence is not the last column, adjust the file to make it so
            for i in range(len(seqdata)):
                seq = seqdata[i][seqloc] #getting the sequence
                del seqdata[i][seqloc] #removing the sequence from the midle of the line
                seqdata[i].append(seq) #adding it to the end of the line
        ### Writing this to a new file name. Note, if the last column was the sequence, we're only renaming the file
        print "writing these to", dataseqfile
        f_seqs=open(dataseqfile, 'w', 1)
        for line in seqdata: 
            print >>f_seqs, ("\t".join(map(str,line)))
        f_seqs.close
        if len(read_data(dataseqfile)) == 0: sys.exit("The sequence file was not created properly: exiting!")
        else: print "sequence file was not empty, using this..."
    else: print "no sequence found in input file"
     
    ### Getting the genome
    if not os.path.isfile(dataseqfile): #If the file doesn't already exist (i.e. wasn't created above, or any time in the past
        print 'Generating the sequences for the data sets...'
        ###getting the sequence data files, then getting the sequences for each from the genome dictionary
        data = read_data(seqfile)
        ### Getting the sequences for the peaks
        genome = load_genome(genomefile)
        for n1 in range(len(data)): #for each line in the data set
            chrom,start,stop = data[n1][0],int(data[n1][1]),int(data[n1][2])
            sequence = genome[chrom][start:stop].upper()
            data[n1] = [chrom, start, stop, sequence] #adding the sequence to the line, to write to a file
        ### Writing the results to new files
        f1=open(dataseqfile, 'w')
        for line in data:print >>f1, "\t".join(map(str,line))
        f1.close()
        genome.clear() #to clear this out of system memory
    else: #If this file already exists, lets just use it instead
        print 'Loading the data sets with sequences generated previously (delete or rename them if you want generate new ones)'
        #print >>f, 'The file'
        #data = read_data(seqfile)
    
    if seqtype == 'fwd': 
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_FWD_SVR-scores.txt'
        print >>f_info, "Finding the scores for only sequences in the forward direction"
    elif seqtype == 'rev': 
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_REV_SVR-scores.txt'
        print >>f_info, "Finding the scores for only sequences in the reverse direction"
    elif seqtype == 'best':
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_SVR-scores.txt'
        print >>f_info, "Finding the best score of both the forward and reverse sequences"
    
    print >>f_info, 'Scores are being written to: ', datascorefile
    print 'Getting scores for this sequence using the SVR model...'
    seqdata = read_data(dataseqfile)
    scores_by_SVR(seqdata,modelfile,datascorefile)
    return datascorefile

def scores_by_SVR(seqdata,modelfile,outfile):
    '''Takes a ChIP file in a format with col1 = chromosome, col2 = seqstart, col3 = seqend, and the last 
    column has the sequence for that peak, gets the predicted intensity for all sub-sequences (length is 
    defined in the model file) in each peak, Now also checks if the outfile already exists, and starts 
    where it left off in case of incomplete run. "outfile" is the name of the desired output file. Make
    sure to edit the lines below when testing what the core sequence so that it matches the models used.
    Now optimiezed for speed'''
    
    if os.path.isfile(outfile):
        if len(read_data(outfile)) == len(seqdata):
            print "Looks like this file already has sequences. Delete this file and re-run to start over:", outfile
            return 
    
    ### Breaking the allseqs file into smaller chunks if it's too big (500,000 sequences per chunk)
    print "total number of sequences", len(seqdata)
    
    ### Splitting list into smaller lists, and saving as a list of lists
    chiplists = [seqdata[i:i + chunksize] for i in range(0, len(seqdata), chunksize)]
    
    chunks = len(chiplists)
    if chunks > 1: print "Splitting this into", chunks, "files, each with", chunksize, "main sequences"
    #else: print "This is small enough, using whole file"
    
    f=open(outfile,'w', 1)
    for n in range(len(chiplists)): #for each 
        #print "testing set", n+1, "of", len(chiplists)
        smallchip = chiplists[n]
    
        ### Making master list of sequences to get scored
        allseqs = [] #resulting colums = chromosome, peakstart, peakend, kmer-number, kmer-seq
        for line in smallchip:
            peak = line[-1] #getting the peak sequence
            for n2 in range(len(peak)-(length-1)): #for ever possible kmer (with length of the model) in the larger sequence
                kmer = peak[n2:n2+length]
                allseqs.append(line[:3]+[n2,kmer]) #adding this new sub-sequence to the list, along with relevant info

        ### Finding the kmers (and reverse complements) with good cores, writing to a new list
        seqlist = []
        goodcorekmers = []
        searchstrings = ['GCGG','CCGC','GCGC']
#         searchstrings = ['']
        for line in allseqs:
            kmer = line[-1]
            coreF = kmer[(length/2)-2:(length/2)+2] #finding the central core 4-mer
            kmerR = reverse_complement(kmer)
            coreR = kmerR[(length/2)-2:(length/2)+2] #finding the central core 4-mer for reverse
            if any(coreseq in coreF for coreseq in searchstrings): seqlist.append(kmer) #if the core is a good core, add to the list
            if any(coreseq in coreR for coreseq in searchstrings): seqlist.append(kmerR) #if the core is a good core, add to the list
         
        print "  Getting the SVR scores for", len(seqlist), "small sequences:  set", n+1, "of", len(chiplists)
        scores = apply_model_to_seqs(seqlist,modelfile) #columns are Score, Sequence, with header line
        #for line in scores: print "\t".join(map(str,line))
         
        ### Assigning SVR scores (or badscore) to allseqs list
        #print "Assigning scores to the right sequences"
        starttime2 = time.time()
        for i in range(len(allseqs)): # for each 36mer in each peak
            chrom, start, stop, kmerpos, seq1 = allseqs[i]
            #print allseqs[i][4], len(scores)
            ### Getting the score. Note, kmers that were scored have reverse compliments, getting best of those two scores, or badscore
            newscore = badscore #we'll change this if there is a score assigned to this sequence
            if len(scores) > 1:
                for j in range(1,(len(scores)),2): #for each sequence we have a SVR score for, counting by 2 
                    #print "lines", j, "and", j+1, "out of", len(scores)-1, scores[j], scores[j+1]
                    score_f, seq_f = scores[j]
                    score_r, seq_r = scores[j+1]
                    if seq1 == seq_f or seq1 == seq_r: 
                        #print "This sequence has a score!", score_f, score_r, seq_f
                        
                        #changethis (remove line below) when done testing
                        if seqtype == 'fwd': newscore = float(score_f)
                        elif seqtype == 'rev': newscore = float(score_r)
                        elif seqtype == 'best': newscore = max([float(score_f),float(score_r)])
                        else: sys.exit("seqtype parameter in the program is not valid; must be 'fwd', 'rev', or 'best' ")
                        
                        del scores[j:j+2]
                        break #once we find a score for this sequence, exit out of the "if seq1 ==  seq_f" loop
                    else: 
                        continue
                    break #once we've found a score
            if newscore < 0: newscore = 0
            if newscore > 1: newscore = 1
            allseqs[i].append(newscore)
            #print "score is", newscore, "for", allseqs[i]
            
        #for line in allseqs: print line
     
        ### Formatting results and writing to output file
#         print "Formatting results for part", n+1,"out of", len(chiplists), "and updating output file:", outfile
        lineout = []
        for n2 in range(len(allseqs)):
#             print n2+1, "out of", (len(allseqs))
            line = allseqs[n2]
            chrom, start, stop, kmerpos, seq1 = line[:5]
            if len(lineout) == 0: #starting the new line, with this first score
                lineout = [chrom, start, stop, line[-1]]
            else: #adding the score to the line
                lineout.append(line[-1]) #appending the score
            if n2 == len(allseqs) - 1: #if this is the last line, write the outline to the file, and quit 
                print >>f, "\t".join(map(str,lineout))
                break #if this is the last line, stop
            if line[:3] != allseqs[n2+1][:3]: #if the next line is for a different peak, write the results to the file, and start a new line
                print >>f, "\t".join(map(str,lineout))
                lineout = []
        
        endtime2 = time.time()
        runtime2 = int(endtime2 - starttime2)
    f.close()

def apply_model_to_seqs(seqlist,model):
    '''Getting predicted scores for set of sequences'''
    starttime = time.time()
    
    ### Generating the matrix file
    svrmatrix = []
    #print 'Getting matrix for data set, with', len(seqlist), 'sequences'
    for seq in seqlist:
        #score,seq = 0,line #assigning a default score of 0 to every sequence
        ###Creating the list of features with relevent info
        featureinfo = [['feature','position','featnum','featvalue'],['start','na',1,1]] #header, and first feature which never changes
        featnum = 2 #the first feature is already definde as 1:1, so we start with 2
        #print seq, score
        for k in kmers:
            ###Getting list of all possible combinates of bases, of length k, and all 4mers (for core)
            bases = [] #creating empty lists needed later
            for n in itertools.product('ACGT', repeat=k): bases.append(''.join(n)) #all possible kmers with lenth k
            for n1 in range(len(seq)-(k-1)): #For each position in the sequence
                kmer =seq[n1:n1+k] #getting the actual kmer at this position
                for n2 in range(len(bases)): #for every possible kmer of the size we want
                    feature = bases[n2] #getting the feature
                    if feature == kmer: featvalue = 1 #testing if the actual k-mer matches the feature, in which case the feature value is 1
                    else: featvalue = 0
                    featureinfo.append([feature, n1, featnum, featvalue]) #adding the info about the feature to the list
                    featnum += 1 #increasing the feature number by 1
        #for line in featureinfo: print line 
        features = [0] #starting a new list for building the matrix for SVR
        for x in range(1,len(featureinfo)): #for every feature, in this list (skipping first item because header)
            features.append(str(featureinfo[x][2])+':'+str(featureinfo[x][3])) # putting the feature values into the list
        svrmatrix.append(features) #adding the features for each sequence to the master list of features for this set of sequences

    
    ### Writing the matrix to a temporary file
    matrixfile = 'SVRmatrix_temp.txt'
    #print 'Writing matrix to', matrixfile
    f1 = open(matrixfile, 'w')
    for line in svrmatrix: print >>f1, ("\t".join(map(str,line)))
    f1.close()
    
    ### Getting the predicted scores using the SVR model
    outfile = matrixfile[:-4]+'-prediction.txt'
#     print "svm-predict", testfile, modelfile, outfile
    args = ["svm-predict", matrixfile, model, outfile]
    output, error = Popen(args, stdout = PIPE, stderr = PIPE).communicate() #running the command, and storing the results
    if len(error) > 0: #if running the command caused an error, print the error
        print "SVM-predict error:", error
    
    ### Getting and organizing the results
    results = [['Score','Sequence']]
    outdata = read_data(outfile)
    for x in range(len(seqlist)):
        #print seqlist[x], outdata[x]
        results.append([float(outdata[x][0]), seqlist[x]])
        #print outdata[x][0], seqs[x], seqs[x][16:20]
    endtime = time.time()
    runtime = int(endtime - starttime)
    #print "SVR runtime (hh:mm:ss) = ", time.strftime('%H:%M:%S', time.gmtime(runtime))
    return results
    
def SVR_ROC(modelfile,posscorefile,negscorefile):
    '''The goal here is to take a libsvm regression model for binding, a positive set of sequences (i.e. ChIP peaks), 
    and a negative set of sequences (i.e. DNAse not overlapping ChIP peaks), and builds a table for making an ROC plot'''
    
    posscoredata = read_data(posscorefile)
    negscoredata = read_data(negscorefile)
    
    ### Finding the max and min from the two data sets
    high, low = 1.0,0.0
    
    ### For each data set, getting the maximum score for each sequence, and only keeping this as 4th column
    for data in [posscoredata, negscoredata]:
        for i in range(len(data)):
            maxscore = max([float(x) for x in data[i][3:]])
            data[i] = data[i][:3] + [maxscore] #re-writing the line to only include chromosome, start position, stop position, and max score
     
    ### Getting the list of cutoffs to use
    cutoffs = []
    counter = low
    size = 1000 # how many points we want in the list
    for x in range(size):
        cutoffs.append(counter)
        counter += (((high-low)/size) )
    cutoffs.append(high) # we need to be inclusive of the max value for the cutoff
     
    ### Building the TPR/FPR table for every score cutoff
    ROCdata = []
    for cutoff in cutoffs:
        #cutoff = cutoffs[n1]
        TP, FN, FP, TN, = 0,0,0,0
        #print "looking at positive controls"
        for line in posscoredata: # Getting TP and FN
            #print line[3], cutoff
            if float(line[3]) >= cutoff: TP += 1 # If this positive control has a score above the cutoff, then this a true positive 
            else: FN += 1 # If this is below the cutoff, then this is a false negative
        #print "looking at negative controls"
        for line in negscoredata:
            #print line[3]
            if float(line[3]) < cutoff: TN += 1 # If this negative control has a score below the cutoff, then this a true negative 
            else: FP += 1 # If this is above the cutoff, then this is a false positive
        #print cutoff, TP, FN, FP, TN, '\n'
        TPR = float(TP)/(float(TP)+float(FN))
        FPR = float(FP)/(float(FP)+float(TN))
        ROCdata.append([cutoff, TPR, FPR])
    #for line in ROCdata: print line
    outfile = outprefix+'_SVR_ROC-curve-data.txt'
    f=open(outfile, 'w')
    for line in ROCdata: print >>f, "\t".join(map(str,line))
    f.close
    
    auc = AUC(ROCdata)
    print 'AUC is', auc
    print >>f_info, 'AUC is', auc
    

def binding_prediction_PWM_chipall(pwmfile,seqfile,genomefile):
    '''Takes a ChIP peak file, uses a PWM to predict the binding intensity
    for every kmer (length according to the model) in each peak, returns all scores.'''
    dataseqfile = os.path.splitext(seqfile)[0]+'_sequences.txt' #getting the file name
    
    ### Checking if this file already contains sequences, so we don't have to bother getting them
    seqdata = read_data(seqfile)
    linetest = seqdata[0] #getting the first line from the file for testing
    seqloc = 0
    
    if len(linetest) > 3:
        for col in linetest[3:]: #for each column, starting with the 4th (cant be in first three)
            if not re.search(r'[^ATGC]', col) and len(col) == int(linetest[2]) - int(linetest[1]): # if this column contains only A, T, G, or C, and is the right length
                    seqloc = linetest.index(col) #getting the column number
    if seqloc > 0: # if this exists, then we can use the sequence listed in this file
        print "\nInput file already contains sequences, using these."
        if seqloc != len(linetest)-1: #if sequence is not the last column, adjust the file to make it so
            for i in range(len(seqdata)):
                seq = seqdata[i][seqloc] #getting the sequence
                del seqdata[i][seqloc] #removing the sequence from the midle of the line
                seqdata[i].append(seq) #adding it to the end of the line
        ### Writing this to a new file name. Note, if the last column was the sequence, we're only renaming the file
        print "writing these to", dataseqfile
        f_seqs=open(dataseqfile, 'w', 1)
        for line in seqdata: 
            print >>f_seqs, ("\t".join(map(str,line)))
        f_seqs.close
        if len(read_data(dataseqfile)) == 0: sys.exit("The sequence file was not created properly: exiting!")
        else: print "sequence file was not empty, using this..."
    else: print "no sequence found in input file"
     
    ### Creating the sequence file
    if not os.path.isfile(dataseqfile): #If the file doesn't already exist (i.e. wasn't created above, or any time in the past
        print 'Generating the sequences for the data sets...'
        ###getting the sequence data files, then getting the sequences for each from the genome dictionary
        data = read_data(seqfile)
        ### Getting the sequences for the peaks
        genome = load_genome(genomefile)
        for n1 in range(len(data)): #for each line in the data set
            chrom,start,stop = data[n1][0],int(data[n1][1]),int(data[n1][2])
            sequence = genome[chrom][start:stop].upper()
            data[n1] = [chrom, start, stop, sequence] #adding the sequence to the line, to write to a file
        ### Writing the results to new files
        f1=open(dataseqfile, 'w')
        for line in data:print >>f1, "\t".join(map(str,line))
        f1.close()
        genome.clear() #to clear this out of system memory
    else: #If this file already exists, lets just use it instead
        print 'Loading the data sets with sequences generated previously (delete or rename them if you want generate new ones)'
        #print >>f, 'The file'
        #data = read_data(seqfile)
    
    if seqtype == 'fwd': 
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_FWD_PWM-scores.txt'
        print >>f_info, "Finding the scores for only sequences in the forward direction"
    elif seqtype == 'rev': 
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_REV_PWM-scores.txt'
        print >>f_info, "Finding the scores for only sequences in the reverse direction"
    elif seqtype == 'best':
        datascorefile = outprefix+'_'+os.path.splitext(modelfile)[0]+'_PWM-scores.txt'
        print >>f_info, "Finding the best score of both the forward and reverse sequences"
    
    print >>f_info, 'Scores are being written to: ', datascorefile
    print 'Getting scores for sequences using PWM...'
    peakscores_by_PWM(read_data(dataseqfile),pwmfile,datascorefile)
    return datascorefile

def peakscores_by_PWM(chipdata,pwmfile,outfile):
    '''Takes a ChIP file in a format with col1 = chromosome, col2 = seqstart, col3 = seqend, and the last 
    column has the sequence for that peak, gets the predicted intensity for all sub-sequences (length is 
    defined in the model file) in each peak, Now also checks if the outfile already exists, and starts 
    where it left off in case of incomplete run. "outfile" is the name of the desired output file. Make
    sure to edit the lines below when testing what the core sequence so that it matches the models used.'''
    
    ### Checking if outfile exists, and if incomplete, to continue where it left off
    linestart = 0
    if os.path.isfile(outfile): #if file exists
        print "Warning!  Output file already exists!"
        #print "desired length:", len(chipdata)+1, "actual length", len(read_data(outfile))
        if (len(read_data(outfile))) < (len(chipdata)): # if file is incomplete - note, outfile has 2 lines for every peak, thus dividing by 2
            linestart = len(read_data(outfile)) #the number of completed lines
            print "Restarting incomplete run at line", linestart
            f=open(outfile, 'a', 1) #append to file, to continue where it left off
        else: 
            print "File is already complete, delete it if you want to re-run.", outfile 
            return #exiting this module because the file is already done
    else: 
        print "Writing results to:", outfile
        f=open(outfile, 'w', 1) #if file doesn't exist, create a new one
    
    for n1 in range(linestart,len(chipdata)):
        #print "getting sequence scores for", n1+1, "out of", len(chipdata)
        line = chipdata[n1]
        chrom, start, stop = line[:3] #first three columns, chromosome, start position, stop position
        peak = line[-1] #last column (sequence)
        allscores = pwm_big_seq_score(pwmfile,peak)
        ### creating the output line, with all the scores
        lineout = [line[0],line[1],line[2]] #columns are Chromosome, Start, Stop (can add sequence here, if needed, as "line[3]")
        for line2 in allscores: 
            lineout.append(line2[0])
        print >>f, "\t".join(map(str,lineout)) #writing this to the file
    f.close()
        
def pwm_big_seq_score(pwmfile,sequence):
    '''Takes a sequence, and gets the PWM score for all kmers the size of the pwm. Note, gets the PWM score for the 
    reverse complement, and reports the larger of the two for that position'''
    import sys, re
    pwm = read_pwm(pwmfile)
    pwmsize = len(pwm[0])-1
#     print "sequence is", sequence
    #print "PWM size is", pwmsize, "so there will be", (len(sequence)-pwmsize), "scores. (should be", (150-pwmsize), "scores)"
    results = []
    for i in range(len(sequence)-pwmsize+1):
        kmer = sequence[i:(i+pwmsize)]
#         print "kmer is", kmer
        baseorder = [pwm[0][0],pwm[1][0],pwm[2][0],pwm[3][0]] #getting the order of the bases in the PWM
        if re.search(r'[^ATGC]', kmer): #if the kmer contains something other than A, T, G, or C...
            break
        else:
            kmerR = reverse_complement(kmer)
            score_f,score_r = 0.0,0.0
            for i in range(len(kmer)):
                base = kmer[i]
                pwmline = baseorder.index(base)
                basevalue = pwm[pwmline][i+1]
                score_f += basevalue
            for i in range(len(kmerR)): # testing the reverse complement
                base = kmerR[i]
                pwmline = baseorder.index(base)
                basevalue = pwm[pwmline][i+1]
                score_r += basevalue
            
#             score = max(scoreF,scoreR) #using the highest PWM score from the forward or reverse sequence
            #changethis (remove line below) when done testing
            if seqtype == 'fwd': newscore = float(score_f)
            elif seqtype == 'rev': newscore = float(score_r)
            elif seqtype == 'best': newscore = max([float(score_f),float(score_r)])
            else: sys.exit("seqtype parameter in the program is not valid; must be 'fwd', 'rev', or 'best' ")
            
#             print "fwd:", kmer, scoreF, "rev:", kmerR, scoreR, "max:", score
                #print basevalue, "score so far", score
#         print "...score is", score
        results.append([newscore,kmer])
    return results

def scores_by_PWM(chipdata,pwmfile,outfile):
    '''Takes a ChIP file in a format with col1 = chromosome, col2 = seqstart, col3 = seqend, and the last 
    column has the sequence for that peak, gets the predicted intensity for all sub-sequences (length is 
    defined in the model file) in each peak, Now also checks if the outfile already exists, and starts 
    where it left off in case of incomplete run. "outfile" is the name of the desired output file. Make
    sure to edit the lines below when testing what the core sequence so that it matches the models used.'''
    
    ### Checking if outfile exists, and if incomplete, to continue where it left off
    linestart = 0
    if os.path.isfile(outfile): #if file exists
        print "Warning!  Output file already exists!"
        if (len(read_data(outfile))) < (len(chipdata)): # if file is incomplete 
            linestart = len(read_data(outfile)) #the number of completed lines
            print "Restarting incomplete run at line", linestart
            f=open(outfile, 'a', 1) #append to file, to continue where it left off
        else: 
            print "File is already complete, delete it if you want to re-run.", outfile 
            return #exiting this module because the file is already done
    else: 
        print "Writing results to:", outfile
        f=open(outfile, 'w', 1) #if file doesn't exist, create a new one
    
    
    for n1 in range(linestart,len(chipdata)):
        #print "getting sequence scores for", n1+1, "out of", len(chipdata)
        line = chipdata[n1]
        chrom, start, stop = line[:3] #first three columns, chromosome, start position, stop position
        peak = line[-1] #last column (sequence)
        allscores = pwm_big_seq_score(pwmfile,peak)
        ### creating the output line, with all the scores
        lineout = [line[0],line[1],line[2]] #columns are Chromosome, Start, Stop (can add sequence here, if needed, as "line[3]")
        for line2 in allscores: 
            lineout.append(line2[0])
        print >>f, "\t".join(map(str,lineout)) #writing this to the file
    f.close()
        
def PWM_ROC(modelfile,posscorefile,negscorefile):
    '''The goal here is to take a libsvm regression model for binding, a positive set of sequences (i.e. ChIP peaks), 
    and a negative set of sequences (i.e. DNAse not overlapping ChIP peaks), and builds a table for making an ROC plot'''
    
    posscoredata = read_data(posscorefile)
    negscoredata = read_data(negscorefile)
    
    
    ### For each data set, getting the maximum score for each sequence, and only keeping this as 4th column
    print "\n\nfinding the max score"
    for data in [posscoredata, negscoredata]:
        for i in range(len(data)):
            try: maxscore = max([float(x) for x in data[i][3:]])
            except: 
                print "an error occured finding the max PWm score:", i
                print "scores are", [data[i][3:]]
            data[i] = data[i][:3] + [maxscore] #re-writing the line to only include chromosome, start position, stop position, and max score
    
    ### Finding the max and min from the two data sets
    high, low = -99,99
    for data in [posscoredata, negscoredata]:
        for line in data:
            score = float(line[3])
            if score > float(high): high = score
            if score < float(low) and score != 0: low = score
    print "For ROC, high score is", high, "and low score is", low
     
    ### Getting the list of cutoffs to use
    cutoffs = []
    counter = low
    size = 1000 # how many points we want in the list
    
    for x in range(size):
        cutoffs.append(counter)
        counter += (((high-low)/size) )
    cutoffs.append(high) # we need to be inclusive of the max value for the cutoff
     
    ### Building the TPR/FPR table for every score cutoff
    ROCdata = []
    for cutoff in cutoffs:
        #cutoff = cutoffs[n1]
        TP, FN, FP, TN, = 0,0,0,0
        #print "looking at positive controls"
        for line in posscoredata: # Getting TP and FN
            #print line[3], cutoff
            if float(line[3]) >= cutoff: TP += 1 # If this positive control has a score above the cutoff, then this a true positive 
            else: FN += 1 # If this is below the cutoff, then this is a false negative
        #print "looking at negative controls"
        for line in negscoredata:
            #print line[3]
            if float(line[3]) < cutoff: TN += 1 # If this negative control has a score below the cutoff, then this a true negative 
            else: FP += 1 # If this is above the cutoff, then this is a false positive
#         print cutoff, TP, FN, FP, TN, '\n'
        TPR = float(TP)/(float(TP)+float(FN))
        FPR = float(FP)/(float(FP)+float(TN))
        ROCdata.append([cutoff, TPR, FPR])
    #for line in ROCdata: print line
    outfile = outprefix+'_PWMscores_ROCdata.txt'
    f=open(outfile, 'w')
    for line in ROCdata: print >>f, "\t".join(map(str,line))
    f.close
    
    auc = AUC(ROCdata)
    print 'AUC is', auc    
    print >>f_info, 'AUC is', auc
    
def AUC(data):  
    '''Taking a set of x,y data points, and calculating the area under the curve from these points'''
    AUC = 0
    for i in range(len(data)-1): #needs to be -1 because (i+1) for the last point won't work.
        ###Approximating the area for this point
        x1 = float(data[i][2])
        y1 = float(data[i][1])
        x2 = float(data[i+1][2])
        y2 = float(data[i+1][1])
        area = abs(x2-x1) * ( (y1+y2)/2 )
        AUC += area
    #print 'AUC is', AUC
    return AUC




'''Running the program ==========================================================='''

# ### Finding the length of the sequence used for the model, and feature kmer (1mers, 2mers, 3mers, etc)
# modeldata = read_data(modelfile)
# featnum = float(len(modeldata[7][0].split(' '))-2)
# seqlengths = []
# for n in range(1,6): # for possible feature lengths (i.e. 1mers, 2mers, 3mers, etc)
#     length = (featnum/4**n)+(n-1)
#     if length <= 36 and length%2 == 0:
#         seqlengths.append([int(length),n])
# if len(seqlengths) == 1:
#     print "Length of sequence is", seqlengths[0][0], ": features are kmers of length", seqlengths[0][1]
#     length,k = seqlengths[0]
# elif len(seqlengths) == 0:
#     print "Sequence length and feature kmer length cannot be determined given a total of", featnum, "features. Check the libsvm model file for errors."
# else: 
#     print "\nMultiple possibilities for sequence length exist. Choose the appropriate one."
#     for j in range(len(seqlengths)):
#         print "["+str(j+1)+"]",  "Length of sequence", seqlengths[j][0], ": feature kmers of length", seqlengths[j][1]
#     while True:
#         try: 
#             answer = int(raw_input("\nChoose the best option from above:\n"))
#         except ValueError:
#             print "Please enter an appropriate number"
#             continue
#         if answer > len(seqlengths) or answer <= 0:
#             print "Please enter an appropriate number"
#             continue
#         else:
#             break
#     length,k = seqlengths[answer-1]
# print "Length is ", length, ": k is", k
    
if runtype == 'SVR':
    seqscorefile = binding_prediction_chipall(modelfile,seqfile,genomefile)
    if args.negseqfile:
        print "\nNormalizing the negative sequence file to the input sequence file..."
        negseqfilenorm = normalize_peak_lengths(seqfile,negseqfile)[0]
        print "\nRunning the model prediction on the negative sequence file..."
        outprefix = outprefix + '_negseq'
        negscorefile = binding_prediction_chipall(modelfile,negseqfilenorm,genomefile)
        print "\nGetting the ROC curve"
        SVR_ROC(modelfile,seqscorefile,negscorefile)

if runtype == 'PWM':
    seqscorefile = binding_prediction_PWM_chipall(modelfile,seqfile,genomefile)
    if args.negseqfile:
        print "\nNormalizing the negative sequence file to the input sequence file..."
        negseqfilenorm = normalize_peak_lengths(seqfile,negseqfile)[0]
        print "\nRunning the model prediction on the negative sequence file..."
        outprefix = outprefix + '_negseq'
        negscorefile = binding_prediction_PWM_chipall(modelfile,negseqfilenorm,genomefile)
        print "\nGetting the ROC curve"
        PWM_ROC(modelfile,seqscorefile,negscorefile)



f_info.close()
