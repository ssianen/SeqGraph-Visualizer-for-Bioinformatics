#!/usr/bin/env python

#SeqGraph.py

import sys
import re
import graphviz

#Global Variables
p1_filename = sys.argv[0]
reference_genome = sys.argv[1]


def main():
    '''
    Input: Small FASTA file of sequenced reads, flags where the flags are the following options:
           -O for overlap graphs
           -D for deBruijn graphs
    Output: Overlap graph or deBruijn graph written to new files

    Examples:

    $ python SeqGraph.py testOData.fa -O 
    $ ls overlapG.png

    $ python SeqGraph.py testDBData.fa -D
    $ ls deBruijnG.png
    '''

    acceptedOptions = ["-O", "-D"]
    errorMessage = "Input not recognized. Please double check your input and try again."
  
    if (len(sys.argv) == 3):
        graphType_1 = sys.argv[2]

        if (graphType_1 in acceptedOptions):
            if (graphType_1 == "-O"):
                overlap_graph()
            elif (graphType_1 == "-D"):
                deBruijn_graph()

        else:
            print(errorMessage)
            
    else:
       print(errorMessage)

    return


###########################################################################
#For building overlap graphs

def overlap_graph():
   '''
   Generates an overlap graph from the sequence reads in the given FASTA file.
   '''
   
   completedOMsg = "Overlap graph successfully created."

   dictOfEdges = {} #{read1: [(read2, len(match2)), (read3, len(match3)), .... (readN, len(matchN))], ...}
   listOfReads = []

   listOfReads = process_readsO()

   for readToCheck in listOfReads:
      for otherRead in listOfReads:
         if readToCheck != otherRead:
            overlap_match_finder(readToCheck, otherRead, dictOfEdges)

   build_Ograph(dictOfEdges)
   print(completedOMsg)
   
   return


def process_readsO():
   '''
   Retrieve all of the reads from the FASTA file and append them to a list
   '''
   readsList = []

   openFile = open(reference_genome)

   for line in openFile:
      if (line.find("@") == -1 and line.find("+") == -1 and line.find("h") == -1): #extract just the read
         line = line.strip()
         readsList.append(line)

   openFile.close()
   return readsList


def overlap_match_finder(read1, read2, edgeDictionary):

   perfectMatch = ""

   indexR1 = -10
   indexR2 = 10

   #Initialize the suffix and prefix to check the first 10
   read1Suffix = read1[indexR1:] #minimum perfect match of 10
   read2Prefix = read2[0:indexR2]

   #Case 1: len(perfectMatch) == 10
   if read1Suffix == read2Prefix:
      perfectMatch = read1Suffix #maximal perfect match 

   else:
      #Case 2: len(perfectMatch) > 10
      read2Suffix = read2[indexR2:]
      for indexR2 in range(10, len(read2Suffix)):
         char2 = read2[indexR2]

         if read1Suffix != read2Prefix:
            indexR1 -= 1
            read1Suffix = read1[indexR1:]
            read2Prefix += char2 #Add the next char to the right of read2

         else:
            perfectMatch = read1Suffix
            break

   #If a match >= 10 exists, add it to the edge dictionary
   if (len(perfectMatch) >= 10):
      if read1 not in edgeDictionary:
         edgeDictionary[read1] = []
      edgeDictionary[read1].append( (read2, len(perfectMatch)) )

   return
   

def build_Ograph(edgeDictionary):
   
   #Instantiate digraph object
   overlapGraph = graphviz.Digraph()

   for readNode in edgeDictionary:
      tupleList = edgeDictionary[readNode]

      for readTuple in tupleList:
         otherRead = readTuple[0]
         edgeWeight = readTuple[1]
         overlapGraph.edge(readNode, otherRead, label = " " + str(edgeWeight))
   
   # Save to file and generate a .png
   overlapGraph.render('overlapG', format='png', cleanup=True)  
   
   return
   

###########################################################################
#For building deBruijn graphs

def deBruijn_graph():
    '''
    Generates a de Bruijn graph from the sequence reads in the given FASTA file
    '''

    completedDBMsg = "de Bruijn graph successfully created."

    dictOfEdges = {} #{read1: [(read2, len(match2)), (read3, len(match3)), .... (readN, len(matchN))], ...}
    listOfTuples = []
    listOf10Mers = []

    matchListStr = []

    listOf10Mers = process_readsDB()

    nineMer1 = ""
    nineMer2 = ""

    for index in range(0, len(listOf10Mers) - 1): 
        # Without this -1 on the list length, the output overlap graph looks correct.
        # However, without the - 1, the output overlap graph doesn't match the specification and adds 
        # an extra edge between GAAAGAAAG and AAAGAAAGA from testDBData.fa. 
        # I suspect that this fixes it because the incomplete 10mers that were 
        # found from the input file, belong to these nodes.
        
        tenMer = listOf10Mers[index]
        nineMerPrefix = tenMer[0:9]
        nineMerSuffix = tenMer[1:10]

        dictOfEdges = deBruijn_match_finder(nineMerPrefix, nineMerSuffix, dictOfEdges)

    build_dBgraph(dictOfEdges)
    print(completedDBMsg)
   
    return


def process_readsDB():
    '''
    Get the reads from the FASTA file and append them to a list of 10 mers
    '''

    readCharList = []
    lineAsList = []
    k10MerList = []
    k9MerList = []   #prefixes and suffixes

    tenMerToNineMer = {}

    k10MerStr = ""
    k9MerStr = ""

    openFile = open(reference_genome)

    #Extract just the read
    for line in openFile:
        if (line.find("@") == -1 and line.find("+") == -1 and line.find("h") == -1 and line.find(">") == -1):
            line = line.strip()
            lineAsList = list(line)
            readCharList.extend(lineAsList)

    openFile.close()
   
   #Chop the read string into 10-mers
    for startPointer in range(0, len(readCharList)):
        kmer = readCharList[ startPointer : startPointer + 10 ]
        k10MerStr = ''.join(kmer)
        
        if len(k10MerStr) == 10:
            k10MerList.append(k10MerStr)

    return k10MerList


def deBruijn_match_finder(nineMer1, nineMer2, edgeDictionary):

    suffixOfRead1 = ""
    prefixOfRead2 = ""
    match = ""
 
    suffixOfRead1 = nineMer1[1:9]
    prefixOfRead2 = nineMer2[0:8]

    if (suffixOfRead1 == prefixOfRead2):
        match = suffixOfRead1

        if len(match) == 8 and nineMer1 != nineMer2:

            if (nineMer1, nineMer2) not in edgeDictionary:
                edgeDictionary[(nineMer1, nineMer2)] = 0
            
            edgeDictionary[(nineMer1, nineMer2)] += 1

    return edgeDictionary


def build_dBgraph(edgeDictionary):
   
    #Instantiate digraph object
    deBruijnGraph = graphviz.Digraph()

    for edgePair in edgeDictionary:
        nineMer1 = edgePair[0]
        nineMer2 = edgePair[1]
        occurenceCount = edgeDictionary[edgePair]

        #Loop through the match tuple list value of each read
        for index in range(0, occurenceCount):
            deBruijnGraph.edge(nineMer1, nineMer2)
   
    # Save to file and generate a .png
    deBruijnGraph.render('deBruijnG', format='png', cleanup=True)  
   
    return
   
###########################################################################

if __name__ == "__main__":
   main()