
# coding: utf-8

# In[1]:


import math;
import itertools;
from collections import Counter


# Reading text files

# In[2]:


def writeOutput(outputPatterns):
    cLength = 1
    output_text = ''
    while True:
        curPattern = list()
        num = 0
        for pattern in outputPatterns:
            if getLength(pattern[0]) == cLength:
                curPattern.append(pattern)
                num = num + 1
        print('The number of length ' + str(cLength) + ' sequential patterns is ' + str(num))
        print(curPattern)
        
        cLength = cLength + 1
        if not curPattern:
            break;

    
        
#     print(output_text)
#     output_file = open(out_file, 'w')
#     output_file.write(output_text)
#     output_file.close()
        


# In[3]:


def getLength(pattern):
    while isinstance(pattern[0], list):
        pattern = list(itertools.chain(*pattern))
#     print(pattern)
    return len(pattern)


# Defining MSPS

# In[4]:


def msps(sequences, minSup, sdc):
        # calculating actual Supports
    actualSupport = dict()
    sequencesLength = len(sequences)
    # print(sequencesLength)
    for seq in sequences:
        for key, value in minSup.items():
            for s in seq:
                if key in s:
                    actualSupport[key] = actualSupport.get(key,0) + 1
                    break;

    for key in minSup.keys():
        if key not in actualSupport.keys():
            actualSupport[key] = 0
    # print(actualSupport)

    # finding frequent items
    frequentItems = dict()
    for key, value in minSup.items():
        if (float(value) <= actualSupport.get(key)/sequencesLength):
    #         print(key, float(value), actualSupport.get(key)/sequencesLength)
            frequentItems[key] = value;

    # Sort the items by MIS
    frequentItems = sorted(frequentItems, key= frequentItems.get)
    # print(frequentItems)

    # 
    for item in frequentItems:
        supItem = actualSupport.get(item)/sequencesLength

        SeqK = list()
        for sequence in sequences:
            dupSeq = list()
            f = 0
            for seq in sequence:
                dupS = list()
                for el in seq:
                    if(el != ','):
    #                     print(el + str(actualSupport.get(el)/sequencesLength))
    #                     print(item + str(mis))
    #                     print('\n')
                        if (abs(float(actualSupport.get(el)/sequencesLength) - float(supItem)) <= float(sdc)):

                            dupS.append(el)

                dupSeq.append(dupS)
    #             print(dupSeq)
                if item in seq:
                    f = 1
            if f == 1:
    #             print(dupSeq)
                SeqK.append(dupSeq)
        misCount = int(math.ceil(float(minSup.get(item)) * sequencesLength))
        rPrefixSpan(item, SeqK, misCount)

        sequences = removeItem(sequences, item)
        ####### del item from sequences and then continue
    #     rPrefixSpan(Seqk)
    #     print(SeqK)
    #     print('\n')
    #     print('\n')      


# In[5]:


def rPrefixSpan(item, SeqK, misCount):

    [frequentItemsSeq, itemsCount] = removeInfrequentItems(SeqK, misCount)
#     for seq in frequentItemsSeq:
#         print(seq)
#     print(itemsCount)

    len_1_pattern = list()

    for key in itemsCount.keys():
        len_1_pattern.append(key)

#     print(len_1_pattern)
#     print('Item : ' + item)
    for pattern in len_1_pattern:
        if containsBaseItem(pattern, item) == 'true':
            outputPatterns.append([pattern, itemsCount.get(pattern)])


    for prefix in len_1_pattern:
        prefixSpan(prefix, frequentItemsSeq, item, misCount)

    # //prefix = [[x]] || [[x,y]] || [[x][y]]
    # print(outputPatterns)
    # for abc in outputPatterns:
    #     print(abc)


# In[6]:


def removeItem(sequences, item):
    finalSeq = list()
    # print(sequences)
    # print(item)
    for sequence in sequences:
        dupSequence = list()
        for seq in sequence:
            dupSeq = list()
            for s in seq:
                if s != ',':
                    if s != item:
                        dupSeq.append(s)
            dupSequence.append(dupSeq)
        finalSeq.append(dupSequence)
#     print(finalSeq)
    return finalSeq;


# In[7]:


def prefixSpan(prefix, frequentItemsSeq, item, misCount):

    projectedDB = projection(prefix, frequentItemsSeq)

#     print(projectedDB)


    template1 = list() #Two templates {_, x} and {30, x}
    template2 = list() #<{30}{x}>

    for sequence in projectedDB:

        for seq in sequence:

            if '_' in seq:
                if seq[1:len(seq)]!=[]:
                    template1 += seq[1:len(seq)]
            else:
                lastItemSetCount = 0
    #                 while lastItemSetCount < len(prefix[-1]):
    #                     print(lastItemSetCount)
    #                     print('loops2')
    #                     print(prefix[-1][lastItemSetCount])
    #                     print(seq)
    #                     if prefix[-1][lastItemSetCount] in seq:
    #                         lastItemSetCount = 1
    #                     else:
    #                         break;
                for pre in prefix[-1]:
                    if pre in seq:
                        lastItemSetCount = lastItemSetCount + 1


                if lastItemSetCount == len(prefix[-1]):
                    if seq[seq.index(prefix[-1][-1])+1:] !=[]:
                        template1 += seq[seq.index(prefix[-1][-1])+1:]
                else:

                    template2 += seq

    allTemplate1 = list()                    
    for t1 in template1:
        if t1 not in allTemplate1:
            allTemplate1.append(t1)


    allTemplate2 = list()                    
    for t2 in template2:
        if t2 not in allTemplate2:
            allTemplate2.append(t2)

    temp1Dict = dict()
    temp2Dict = dict()




    for t1 in allTemplate1:
        currSeq1 = 0
        while currSeq1 < len(projectedDB):
            for ind, sequence1 in enumerate(projectedDB[currSeq1]):
    #             print(t1)
    #             print(projectedDB[currSeq1])
                if currSeq1 < len(projectedDB):
                    if '_' in sequence1 and t1 in sequence1:
    #                     print(sequence1)

                        temp1Dict[t1] = temp1Dict.get(t1, 0) + 1
                        currSeq1 = currSeq1 + 1
                        break;
                    else:
                        lastItemSetCount = 0
    #                     print(sequence1)
    #                         while lastItemSetCount < len(prefix[-1]):

    #                             if prefix[-1][lastItemSetCount] in sequence1:
    #                                 lastItemSetCount = 1
    #                             else:
    #                                 break;
                        for pre in prefix[-1]:
                            if pre in sequence1:
                                lastItemSetCount = lastItemSetCount + 1

                        if t1 in sequence1 and lastItemSetCount == len(prefix[-1]):

                            temp1Dict[t1] = temp1Dict.get(t1, 0) + 1
                            currSeq1 = currSeq1 + 1
                            break;
    #             print(currSeq1, ind)
                ind = ind+1;
                if ind == len(projectedDB[currSeq1]) :
                    currSeq1 = currSeq1 +1

    
    for t2 in allTemplate2:
        if t2 not in temp2Dict.keys():
            currSeq2 = 0
            while(currSeq2 < len(projectedDB)):
                for i,sequence2 in enumerate(projectedDB[currSeq2]):
        #             for sequence2 in sequences2:
        #                 print('seq')
        #             print(sequence2)
                    if currSeq2 < len(projectedDB):
#                         print(t2)
#                         print(sequence2)
#                         print(i)
#                         print(currSeq2)
                        if '_' not in sequence2:
                            if t2 in sequence2:
                                
#                                 print(temp1Dict.get(t2, 0))
                                temp2Dict[t2] = temp2Dict.get(t2, 0) + 1
                                currSeq2 = currSeq2 + 1
                                break;

                    i = i +1;
                    if i == len(projectedDB[currSeq2]) :
                        currSeq2 = currSeq2 +1
#     print('tempate 2 dict')
#     print(temp2Dict.items())
    frequentSeqPatterns = []

    for key, value in temp1Dict.items():
        if value >= misCount:

    #         print(str(len(prefix)))
            if len(prefix)>1:  
                appendItem = [prefix[:-1],[prefix[-1].append(key)]]
                appendItem = list(itertools.chain.from_iterable(*appendItem))

                frequentSeqPatterns.append((appendItem, value))
#                 print('\n' + '1111')

#                 print((appendItem, value))
            else:
               
                appendItem = [[prefix[-1],key]]
                appendItem = list(itertools.chain.from_iterable(*appendItem))

#                 print((appendItem, value))
                frequentSeqPatterns.append((appendItem, value))

    for key, value in temp2Dict.items():
        if value >= misCount:
            
#             l3 = list()
#             l3.append(prefix)
#             l3.append(key)
#             print(([list(prefix)] + [[key]], value))
            frequentSeqPatterns.append(([list(prefix)] + [[key]], value))
    # print(prefix)    
    # print(frequentSeqPatterns)
    # print('\n')



    for seqPattern in frequentSeqPatterns:
        l = len(seqPattern)
        if containsBaseItem(seqPattern[:l-1], item) == 'true':
            outputPatterns.append((seqPattern[:l-1], seqPattern[l-1]))
    #             print('outputPatterns')
    #             print((seqPattern[0], seqPattern[1]))
#         print('calling preficSpan on ' + str(seqPattern[:l-1]))
        prefixSpan(seqPattern[:l-1], projectedDB, item, misCount)  # Call prefix_span recursively with the pattern as prefix





# In[8]:


def projection(prefix, frequentItems):
#     print(frequentItems)
    projectedDB = []
#     print(prefix[-1][-1][-1])
    p = prefix[-1][-1][-1]
    currentSeq = 0
    
   
    while (currentSeq < len(frequentItems)):    
        for index, seq in enumerate(frequentItems[currentSeq]):
            
            if currentSeq < len(frequentItems):
                if len(p) == 1:
#                     print(seq)
                    if p in seq:
                        
#                         print(frequentItems[currentSeq][index:])
                        if getPSequence(p, frequentItems[currentSeq][index:]):
                            projectedDB.append(getPSequence(p, frequentItems[currentSeq][index:]))
                        currentSeq = currentSeq +1
                        break;
                else:
                    if '_' in seq:
                        
                        if p[-1] in seq:

                            if getPSequence(p[-1], frequentItems[currentSeq][index:]):
                                projectedDB.append(getPSequence(p[-1], frequentItems[currentSeq][index:]))
                            currentSeq = currentSeq +1
                            break;

                    else:
                        if p in seq:

                            if getPSequence(p, frequentItems[currentSeq][index:]):
                                projectedDB.append(getPSequence(p, frequentItems[currentSeq][index:]))
                            currentSeq = currentSeq +1
                            break;
            if index == len(frequentItems[currentSeq])-1:
                
                currentSeq = currentSeq +1
            
    return projectedDB

            


# In[9]:


def getPSequence(p, sequence):
    if p == sequence[0][-1]:
        
        return sequence[1:]
    else:
        
        sequence[0] = ['_'] + sequence[0][sequence[0].index(p)+1:]
        
        return sequence


# In[10]:


def removeInfrequentItems(SeqK, misCount):
    riiSup = dict()
    countSequences = [ list(set(itertools.chain(*sequence))) for sequence in SeqK ]
#     print(SeqK)
    
    for sequences in countSequences:
        for seq in sequences:
            
                riiSup[seq] = riiSup.get(seq,0) + 1
                
    infrequentItems = list()
#     print(misCount)
#     print(riiSup)
    for key, value in riiSup.items():
        if int(value) < int(misCount):
            infrequentItems.append(key)
    
#     print(infrequentItems)
    
    filteredSequence = [[[s for s in seq if riiSup.get(s) >= misCount]for seq in sequences] for sequences in SeqK]
    
#     print(misCount)
    
    finalFilteredSequence = [[seq for seq in sequences if seq != []]for sequences in filteredSequence]
#     print(finalFilteredSequence)

    for i in infrequentItems:
        if i in riiSup.keys():
            del riiSup[i]
            
    return [finalFilteredSequence, riiSup]


# In[13]:


def containsBaseItem(pattern, item):
   
    pattern = list(itertools.chain.from_iterable(*pattern))

    if item in pattern:
        return 'true'
    else:
        return 'false'


# In[14]:


dataFile = open("data.txt", "r").read().replace(' ','').strip()
dataFile = dataFile.replace('<{', '')
dataFile = dataFile.replace('}>', '')
seq = dataFile.strip().split('\n')

outputPatterns = []

sequences = list()
for sequence in seq:
    sequences.append(sequence.split('}{'))
# print(sequences)
# print('\n')


paramFile = open("parameter.txt", "r").read().replace(' ','').strip()

sdc = paramFile.split('\n')[len(paramFile.split('\n'))-1].split('=')[1]

paramFile = paramFile.split('\n')[0:len(paramFile.split())-1]
# print(paramFile)
# paramFile.replace('MIS')

minSup = dict()
for parameters in paramFile:
    parameters = parameters.split('=')
    key = parameters[0].split(')')[0].split('(')[1]
    val = parameters[1]
    minSup[key] = val;


# print(minSup)
# print('\n')
# print(sdc)
# print('\n')
msps(sequences, minSup, sdc)
writeOutput(outputPatterns)


