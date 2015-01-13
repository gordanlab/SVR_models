#!/usr/bin/python

'''
Created on Nov 12, 2014

@author: josh
'''

import argparse, random, string, sys, os, time, itertools, numpy
from subprocess import *
from distutils.util import strtobool
from operator import itemgetter

#===============================================================================
#  Command Line Arguments

parser = argparse.ArgumentParser(description = 'E2F Binding Model')
parser.add_argument('-i', metavar = 'PBMFile', 
                    help = 'The results of a custom PBM experiment with sequences centered by binding site (i.e. using PWM)' , 
                    dest = 'pbmfile' ,
                    required=True)
parser.add_argument('-o', metavar = 'OutFilePrefix', 
                    help = 'Optional, the prefix that all output files will be based on (do not include file extension). See program notes for proper format of this file', 
                    dest = 'outprefix')
parser.add_argument('-g', "--gridsearch", 
                    help = "Flag for running a grid search, if optimal cost and epsilon values are not known" ,
                    action = "store_true")
parser.add_argument('--seqlength', metavar = 'SequenceLength', 
                    help = "Change the length of the PBM sequence (Default is 36, new sequence will remain centered according to original 36mer from PBM data)",
                    dest = "length")
parser.add_argument('--feature', metavar = 'FeatureType',
                    help = "Define the type of features, i.e. 2 for 2mers, 123 for 1, 2, and 3-mers, etc; default is 3 for 3mers",
                    dest = "rawkmer")
parser.add_argument('--extrafiles',
                    help = "Print extra files: including all matrix files, feature definitions (sequence and position), and model sequence files",
                    action = "store_true")
parser.add_argument('-c', metavar = 'SVR_cost',
                    help = 'The cost value input for LibSVM. If running a grid search, this should be a string of numbers in quotes, i.e. "0.05 0.1 0.5".',
                    dest = 'c')
parser.add_argument('-p', metavar = 'SVR_epsilon',
                    help = 'The epsilon value input for LibSVM. If running a grid search, this should be a string of numbers in quotes, i.e. "0.05 0.1 0.5".',
                    dest = 'p')
args = parser.parse_args()
pbmfile = args.pbmfile

### Getting the output prefix name
if args.outprefix: outprefix = args.outprefix
elif args.gridsearch: outprefix = os.path.splitext(pbmfile)[0] + '_SVR-gridsearch'
else: outprefix = os.path.splitext(pbmfile)[0] + '_SVR-model'



if args.extrafiles: 
    extrafiles = 'yes'
else: extrafiles = 'no'

infofile = outprefix + '_run-info.txt'
f_info = open(infofile, 'a', 1)


#===============================================================================

print "\n_____  Running SVR_model_maker _____"




''' Notes ======================================================================
PBM file format: Tab separated file has header line with column names, and 
columns are: Name, ID, Sequence, Orientation1, Orientation2, Best-orientation, 
and Replicate_intensity_difference.



'''

# pbmfile = 'E2F4-cust2-scores.txt'
# libsvm_generate_matrix_multimers_v5(pbmfile,36,3, 'good',12)
# testfile = 'E2F4-cust2-scores_regr-matrix-test_3mer-feat_36mer-seq_all_good-cores_0.5orientdiff_v5_linear_nocoresoutsidecenter12bp.txt'
# trainfile = 'E2F4-cust2-scores_regr-matrix-train_3mer-feat_36mer-seq_all_good-cores_0.5orientdiff_v5_linear_nocoresoutsidecenter12bp.txt'
# libsvm_run(0.1,0.15,trainfile,testfile)
# libsvm_feature_weights_V2(trainfile,3)


''' Called Parameters =========================================================='''
### These parameters are querried by the program

if args.gridsearch:
    if args.c: 
        c_list_in = args.c
        c_list = c_list_in.split()
        c_list = [float(x) for x in c_list]
    else:
        print '\nEnter the different cost values to test, separated by spaces (i.e. "0.01 0.1 1")'
        c_list_in = raw_input(':')
        c_list = c_list_in.split()
        c_list = [float(x) for x in c_list]
    if args.p: 
        p_list_in = args.p
        p_list = p_list_in.split()
        p_list = [float(x) for x in p_list]
    else:
        print '\nEnter the different epsilon values to test, separated by spaces'
        p_list_in = raw_input(':')
        p_list = p_list_in.split()
        p_list = [float(x) for x in p_list]
else:
    # Getting the libsvm cost and epsilon parameters
    if args.c: c = args.c
    else:
        print "\nLibSVM cost needs to be specified (if not known, try 0.01)"
        c = raw_input(':')
    if args.p: p = args.p
    else:
        print "\nLibSVM epsilon needs to be specified (if not known, try 0.01)"
        p = raw_input(':')
    
''' Other Parameters ================================================================='''
### Different parameters used during the process that can be changed from there defaults here

if args.length: length = int(args.length)
else: length = 36 # Assigning the default value
if args.rawkmer: rawkmer = args.rawkmer
else: rawkmer = "3" # Assigning the default value

# Defining what kind of features we want 
kmers = []
for item in str(rawkmer): kmers.append(int(item))
kinfo = ''
for x in kmers: kinfo = kinfo + str(x) + " + " 
kinfo = kinfo[:-3] + " mer features"

# In the case of E2Fs, we want to make sure that good sequences have a good 
# central core, without having a high-affinity 4-mer in the flanking sequences.
# Set this to [''] to use all sequences (i.e. every 4mer is good in the flanks
# and the core)
searchstrings = ['GCGG','CCGC','GCGC']  

# To avoid misalligned sequences in the PBM, we dont want to use sequences that 
# where the difference between the orientation of a sequence is greater than 
# some cutoff. 
orientcutoff = 0.02

# Defines the size of the center of the sequence; used to find the flanks for 
# looking for good sequences (see searchstrings, above).
coresize = 12

# Determining whether we want to use select sequenes based on good E2F cores
# ("good"), or use all sequences ("all")
coretype = "good"

# Number bins we split the sequences into for SVR, where one bin is used for
# testing the model, and the remaining bins are combined for the training sequences.
svrbins = 5


''' Defining the modules ======================================================='''

def yes_no_query(question):
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return strtobool(raw_input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')

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

def reverse_complement(seq):
    '''Takes any sequence, and gets the reverse complement.  Note, only changes A, T, C, or G. Anything else will be left as is.'''
    comp = string.maketrans("ATCG", "TAGC") #creates the translation definitions, required by line below
    rev_comp = (seq[::-1]).translate(comp) #First reverses the sequence, then translates it according to rules above
    return rev_comp

def random_subset(data,n):
    '''Takes any list, randomizes it, then divides the list into two lists, according the fraction "n".
    Note, listA is the fraction from the list we want, listB is the remaining items from the list'''
    random.shuffle(data) #randomizing the list
    k = int(n*len(data)) #Finding out how many items of the list to keep
    listA = data[:k] #New list, fraction "n" of complete list
    listB = data[k:] #New list, remaining part of the list
    return listA, listB

def list_bins(l,bins):
    '''Takes some list l and breaks it up into bins'''
    n = float(len(l))/bins
    return [l[int(n*i):int(n*(i+1))] for i in range(bins)]
    
def libsvm_generate_matrix(seqlist): 
    '''Generates the matrix file from a list of sequences and their scores'''
    svrmatrix = []
    for line in seqlist:
        score,seq = line
        ###Creating the list of features with relevent info
        featureinfo = [['feature','position','featnum','featvalue'],['start','na',1,1]] #header, and first feature which never changes
        featnum = 2 #the first feature is already definde as 1:1, so we start with 2
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
        features = [score] #starting a new list for building the matrix for SVR
        for x in range(1,len(featureinfo)): #for every feature, in this list (skipping first item because header)
            features.append(str(featureinfo[x][2])+':'+str(featureinfo[x][3])) # putting the feature values into the list
        svrmatrix.append(features) #adding the features for each sequence to the master list of features for this set of sequences
    return svrmatrix,featureinfo
    
def libsvm_run_gridsearch(p_list,c_list,pbmfile):
    ''' Using libsvm, runs a grid search varying cost and epsilon, then prints a table of results, followed by the best R squared'''
    test = 0
    results = []
    
    ### Create your own module for selecting sequences from the PBM file for your protein
    print "\nFinding good sequences from the pbmfile to use for SVR"
    seqlist = E2F_SVR_seq_selector(pbmfile)
    
    ### Creating the matrix file
    print "Generating matrix file for grid search..."
    svrmatrix = libsvm_generate_matrix(seqlist)[0]
    matrixfile = outprefix + '_gridsearch_matrix.txt'
    f = open(matrixfile, 'w')
    for line in svrmatrix: print >>f, ("\t".join(map(str,line)))
    f.close()
    
    ### Reading data in libsvm format
    print '\nDoing Libsvm grid search for', pbmfile, "with 5-fold cross validation"
    print >>f_info, '\nDoing Libsvm grid search for', pbmfile, "with 5-fold cross validation"
       
    for c in c_list: #testing different values of C (cost)
        row = []
        pvals = [' ']
        for p in p_list: #testing different values of epsilon, for each different value of C
            start_time = time.time()
            print "\nTesting epsilon =", p, "and cost =", c
            pvals.append(p)
            command = "svm-train -s 3 -v 5 -t 0 -c " + str(c) + " -p " + str(p) + " " + matrixfile #the command we want to run
            args = string.split(command) #we need to split these into individual items
            output, error = Popen(args,stdout = PIPE, stderr = PIPE).communicate() #running the command, and storing the results
            if len(error) > 0: #if running the command caused an error, print the error, then quit the loop
                print error
                break
            out = string.split(output, '\n') #Splitting the output by line
            SCC =  float(string.split(out[-2])[-1]) #the output is a list of strings (lines), the last word in the last line (with info) is the R squared
            print 'RSQ=', SCC
            print >>f_info, 'p=', p,' c=', c, ' RSQ=', SCC
            print >>f_info, 'Round completed in', (time.time() - start_time)/60, 'minutes'
            print 'Round completed in', (time.time() - start_time)/60, 'minutes'
            if SCC > test:
                test = SCC
                best = [SCC,c,p]
            row.append(SCC)
        results.append([c] + row)
    results.insert(0,pvals)
    print >>f_info, '\nTable of results - columns are epsilon, and rows are cost\n'
    for line in results: print >>f_info, ("\t".join(map(str,line)))
    print >>f_info, '\nBest R squared is ', best[0], ' - with c = ', best[1], ' and p = ', best[2]
    print '\nBest R squared is ', best[0], ' - with c = ', best[1], ' and p = ', best[2]

def E2F_SVR_seq_selector2(pbmfile):
    '''For E2F PBMs, takes a pbm file, and selects only those sequences we want to use for SVR modeling'''
    data = read_data(pbmfile)
    if coretype == 'good': coresearch = searchstrings #used for checking if present in the "core" (=good), leave this blank to include any core, i.e. ['']
    elif coretype == 'any': coresearch = ['']
    
    allseqlist = []
    for n1 in range(1,len(data)):  ###for every line, skipping the first(headers)
        seq = data[n1][2] # getting the sequence from the data
        ###  Only selecting those sequences we want using the following criteria.  A new seq_selector module can be made, using any desired criteria.
        if ( any(string in seq[16:20] for string in coresearch)  #if the core is in the coresearch list, continuing on next line
            and ( float(data[n1][6]) > -orientcutoff and float(data[n1][6]) < orientcutoff )  #if the difference in orientation is above a cutoff, continuing on next line
            and ( 'Bound' in data[n1][0] or 'Neg' in data[n1][0] or 'Flank' in data[n1][0] or 'PosCtrl2' in data[n1][0] )  #if the sequence is in one of these categories, continuing on next line
            and all(string2 not in seq[:18-(coresize/2)] for string2 in searchstrings) and all(string3 not in seq[(18)+(coresize/2):] for string3 in searchstrings) ): #if the flanks don't contain...
            allseqlist.append([float(data[n1][5]),seq])
            #print "Using this sequence"
    print "initially found", len(allseqlist), "sequences"
    ###Truncating the sequence if needed
    for i in range(len(allseqlist)):
        if len(allseqlist[i][1]) > length:
            allseqlist[i][1] = allseqlist[i][1][(36-length)/2:((36-length)/2)+length] #truncating the sequence by trimming the edges
    ### Getting the reverse complement for each sequence, but assigning the same score, and adding to the list of sequences
    seqlist2 = []
    for line in allseqlist:
        seqlist2.append([line[0],reverse_complement(line[1])])
    allseqlist = allseqlist + seqlist2
    print "Number of sequences with reverse complements:", len(allseqlist)
    return allseqlist
      
def E2F_SVR_seq_selector(pbmfile):
    seqsize = length
    data = read_data(pbmfile)
    if coretype == 'good': coresearch = searchstrings #used for checking if present in the "core" (=good), leave this blank to include any core, i.e. ['']
    elif coretype == 'any': coresearch = ['']
    
    print "\nFinding matching sequences from the PBM file"
    allseqlist = []
    for n1 in range(1,len(data)):  ###for every line, skipping the first(headers)
        seq = data[n1][2] # getting the sequence from the data
        ###  Only selecting those sequences we want...
        if ( any(string in seq[16:20] for string in coresearch)  #if the core is in the coresearch list, continuing on next line
            and ( float(data[n1][6]) > -orientcutoff and float(data[n1][6]) < orientcutoff )  #if the difference in orientation is above a cutoff, continuing on next line
            and ( 'Bound' in data[n1][0] or 'Neg' in data[n1][0] or 'Flank' in data[n1][0] or 'PosCtrl2' in data[n1][0] )  #if the sequence is in one of these categories, continuing on next line
            and all(string2 not in seq[:18-(coresize/2)] for string2 in searchstrings) and all(string3 not in seq[(18)+(coresize/2):] for string3 in searchstrings) ): #if the flanks don't contain...
            allseqlist.append([float(data[n1][5]),seq])
            #print "Using this sequence"
    ###Truncating the sequence if needed
    print "Initially found", len(allseqlist), "sequences"
    if len(allseqlist[0][1]) != seqsize:
        for i in range(len(allseqlist)):
            allseqlist[i][1] = allseqlist[i][1][(36-seqsize)/2:((36-seqsize)/2)+seqsize] #truncating the sequence by trimming the edges
    ### Getting the reverse complement for each sequence, but assigning the same score, and adding to the list of sequences
    seqlist2 = []
    for line in allseqlist:
        seqlist2.append([line[0],reverse_complement(line[1])])
    allseqlist = allseqlist + seqlist2
    print "with reverse complements:", len(allseqlist)    
    return allseqlist

def libsvm_run(c,p,pbmfile):
    ''' Using libsvm, for running the best set of values (best if obtained from a grid search), using the train and test matrix files'''
    
    ### Create your own module for selecting sequences from the PBM file for your protein
    print "\nFinding good sequences from the pbmfile to use for SVR"
    allseqlist = E2F_SVR_seq_selector(pbmfile)
    allseqfile = outprefix + '_allseqlist.txt'
    f_a = open(allseqfile, 'w')
    for line in allseqlist: print >>f_a, line

    print "Using", len(allseqlist), "sequences (includes reverse complements)"
    traincount = int(float(len(allseqlist)) * ((float(svrbins)-1)/float(svrbins)))
    print "Using about", traincount, "sequences (+/-1) out of", len(allseqlist), "for training the model" #plus or minus 1, for each of the binned sets
    print >>f_info, "Using about", traincount, "sequences (+/-1) out of", len(allseqlist), "for training the model" #plus or minus 1, for each of the binned sets
    
    random.shuffle(allseqlist) #randomizing the list
    seqbins = list_bins(allseqlist,svrbins) 
    rsqinfo = []
    #### Takes the set of bins and uses one for the test sequences, and the rest for training the model
    for x in range(len(seqbins)):
        print "\nRunning SVR on round", x+1, "of",len(seqbins) 
        testseq = seqbins[x]
        trainseq = list(itertools.chain(*(seqbins[:x]+seqbins[x+1:]))) #Note, the remaining bins need to be collapsed into a single list of lists    
        #print >>f_info, 'Number of sequences in the training data set is', len(trainseq), 'out of', len(allseqlist), 'sequences'
        print "Generating the matrix files..."
        trainmatrix, featureinfo = libsvm_generate_matrix(trainseq) 
        testmatrix = libsvm_generate_matrix(testseq)[0]
        
        ### Writting matrix related output files
        trainmatrixfile = outprefix+'_train_matrix'+'_'+str(x+1)+'.txt'
        f1 = open(trainmatrixfile, 'w',1)
        for line in trainmatrix: print >>f1, ("\t".join(map(str,line)))
        f1.close()
        testmatrixfile = outprefix+'_test_matrix'+'_'+str(x+1)+'.txt'
        f2 = open(testmatrixfile, 'w',1)
        for line in testmatrix: print >>f2, ("\t".join(map(str,line)))
        f2.close()
        
        ###Training the model
        print 'Training the model for run', x+1, '...'
        command = "svm-train -s 3 -t 0 -c " + str(c) + " -p " + str(p) + " " + trainmatrixfile #the command we want to run
        args = string.split(command) #we need to split these into individual items
        output, error = Popen(args, stdout = PIPE, stderr = PIPE).communicate() #running the command, and storing the results
        if len(error) > 0: #if running the command caused an error, print the error
            print error
        out = string.split(output, '\n') #Splitting the output by line
         
        ###testing the model
        print 'Testing the model for run', x+1, '...'
        modelfile = trainmatrixfile + ".model"
        outfile = testmatrixfile[0:-4]+'_SVR-prediction.txt'
        args = ["svm-predict", testmatrixfile, modelfile, outfile]
        output, error = Popen(args, stdout = PIPE, stderr = PIPE).communicate() #running the command, and storing the results
        if len(error) > 0: #if running the command caused an error, print the error
            print error
        out = string.split(output, '\n') #Splitting the output by line
        SCC =  float(string.split(out[-2])[-2]) #the output is a list of strings (lines), the last word in the last line (with info) is the R squared
        print 'Squared Correlation Coefficient for run', x+1, 'is', SCC
        print >>f_info, 'Squared Correlation Coefficient for run', x+1, 'is', SCC
        rsqinfo.append([SCC,x])
        
        ### Printing the extra files if we have the extrafiles flag
        if extrafiles == 'yes':
            if x == 0: #We only need this file once, not for every iteration
                featurefile = outprefix+'_example-feature-list.txt'
                f0 = open(featurefile, 'w', 1)
                for line in featureinfo: print >>f0, ("\t".join(map(str,line)))
                f0.close()
        if SCC >= max([ rsqinfo[i][0] for i in range(len(rsqinfo)) ]):
            ###organizing the results
            testdata = read_data(testmatrixfile)
            predictdata = read_data(outfile)
            if extrafiles == 'yes':
                trainseqfile = outprefix+'_train_sequences.txt'
                testseqfile = outprefix+'_test_sequences.txt'
    
    rsqinfo = sorted(rsqinfo, key=itemgetter(0), reverse=True) #sorting the list by the first column
    
    bestrun = rsqinfo[0][1]
    rsqlist = [ float(rsqinfo[i][0]) for i in range(len(rsqinfo)) ] #getting just the r squared scores (first column) from the list
    mean_rsq = numpy.mean(rsqlist) #getting the average r squared
    std_rsq = numpy.std(rsqlist) #getting the standard deviation of the r squareds
    print "\nBest R squared is", rsqinfo[0][0], "from run", bestrun+1
    print "Average R squared =", mean_rsq, "\nStandard deviation =", std_rsq,
    print >>f_info, "\nBest R squared is", rsqinfo[0][0], "from run", bestrun+1, "- Mean = ", mean_rsq, "- Standard Deviation = ", std_rsq

    if extrafiles == 'yes':
        f3 = open(trainseqfile, 'w', 1)
        for line in trainseq: print >>f3, ("\t".join(map(str,line)))
        f3.close()
        f4 = open(testseqfile, 'w', 1)
        for line in testseq: print >>f4, ("\t".join(map(str,line)))
        f4.close()
        
    ### Keeping only the files for the best run, and renaming them
    for x in range(len(seqbins)):
        trainmatrixfile = outprefix+'_train_matrix'+'_'+str(x+1)+'.txt'
        testmatrixfile = outprefix+'_test_matrix'+'_'+str(x+1)+'.txt'
        modelfile = trainmatrixfile + ".model"
        outfile = testmatrixfile[0:-4]+'_SVR-prediction.txt'
        if x != bestrun:
            os.remove(trainmatrixfile)
            os.remove(testmatrixfile)
            os.remove(modelfile)
            os.remove(outfile)
        else:
            os.rename(trainmatrixfile, trainmatrixfile[:-6]+'.txt')  
            os.rename(testmatrixfile, testmatrixfile[:-6]+'.txt')
            os.rename(modelfile, trainmatrixfile[:-6]+'.txt.model')
            os.rename(outfile, testmatrixfile[0:-6]+'_SVR-prediction.txt')  

    bestresults = [['Actual-Intensity','Predicted-Intensity']]
    for n1 in range(len(testdata)): #for every line, except the first, wich contains headders
        bestresults.append([testdata[n1][0],predictdata[n1][0]])
    resultsfile = outprefix+'_SVR-test_prediction-results.txt'
    f = open(resultsfile, 'w', 1)
    for line in bestresults: print >>f, ("\t".join(map(str,line)))
    f.close
     
    ### Writing info to the info file
    print >>f_info, "\nOutput files for best run are:"
    print >>f_info, ' ', outprefix+'_train_matrix.txt', "<-- The LibSVM matrix file for the sequences used to trian the model"
    print >>f_info, ' ', outprefix+'_test_matrix.txt', "<-- The LibSVM matrix file for the sequences used to test the model" 
    if extrafiles == 'yes':
        print >>f_info, ' ', featurefile, "<-- The LibSVM matrix for the sequences used to trian the model"
        print >>f_info, ' ', outprefix+'_train_sequences.txt', "<-- The actual sequences and the corresponding scores used for training the model" 
        print >>f_info, ' ', outprefix+'_test_sequences.txt', "<-- The actual sequences and the corresponding scores used for testing the model" 
        print >>f_info, ' ', featurefile, "<-- List of definitions for the features (the feature sequence and it's position in the complete sequence" 
    print >>f_info, ' ', outprefix+'_train_matrix.txt.model', '  <-- The model file that can be used by libsvm to predict binding affinities'
    print >>f_info, ' ', outprefix+'_test_matrix_SVR-prediction.txt', '  <-- The actual predicted intensity scores for the test set generated by LibSVM'
    print >>f_info, ' ', resultsfile, '  <-- Contains both the actual intensities, and predicted intensities for the test set'
                
''' Running the program ========================================================'''

### Testing to see if we need to normalize the data in the PBM file
pbmdata = read_data(pbmfile)
scores = [float(row[5]) for row in pbmdata[1:]]
maxval,minval = max(scores), min(scores)
if maxval > 1 or maxval < 0 or minval > 1 or minval < 0: #if the scores are not between 0 and 1, we need to normalize everything
    print >>f_info, "Normalizing the data in the PBM file"
    for i in range(1,len(pbmdata)):
        for j in [3,4,5]: # each of the columns we need to normalize
            score = float(pbmdata[i][j])
            pbmdata[i][j] = (score - minval)/(maxval-minval)
        pbmdata[i][6] = pbmdata[i][3]-pbmdata[i][4] # getting the new difference
    newpbmfile = os.path.splitext(pbmfile)[0] + '_normalized.txt'
    f = open(newpbmfile, 'w')
    for line in pbmdata: print >>f, "\t".join(map(str,line))
    f.close()
    pbmfile = newpbmfile #making sure we only use this new normalized file from now on

print >>f_info, '\n================================================================================\n' , \
'\nProgram started on', time.strftime("%m\%d\%Y %H:%M"), \
'\nGenerating the SVR model file, using libsvm', \
'\n\nParamaters used in this run:', \
'\n ', " ".join(map(str,searchstrings)),  '  <-- Good core 4-mers used for selecting sequences for building the model', \
'\n ', args.pbmfile, '  <-- Input PBM file', \
'\n ', outprefix, '  <-- Output file pfrefix', \
'\n ', length, '  <-- Sequence Length', \
'\n ', kinfo, '  <-- Feature type', \
'\n ', orientcutoff, '  <-- PBM orientation difference cutoff', \
'\n ', coresize, '  <-- Size of central sequence for excluding sequences with good cores in the flanks for building the model', \
'\n ', svrbins, '  <-- Number bins we split the sequences into for SVR, where for each bin, it is used for testing, with the rest used for training the model', \
'\n  Linear   <-- LibSVM support vector regression model type', \
'\n ', c, '  <-- LibSVM cost', \
'\n ', p, '  <-- LibSVM epsilon\n'


if args.gridsearch:
    print "Running a grid search. Use the results and re-run this program to refine another grid search, or run full LibSVM"
    print "\nCost values to be tested:", c_list, "\nEpsilon values to be tested:", p_list
    libsvm_run_gridsearch(p_list,c_list,pbmfile)
else:
    print "Running full libsvm"
    libsvm_run(c,p,pbmfile)

f_info.close()
    
