#!/usr/bin/python3

import os


base = 'benchmarks'
folders = [ "matrix-mpi","phase-retrieval-benchmarks",
            "lid-driven-cavity","FourierBenchmarks",
            "hpcg" ]

metaData = ['minutes','seconds','instructions']
data = dict()

for folder in folders:
    path = base + "/" + folder + "/"
    files = os.listdir(path)
   
    # GET TIMING DATA ###############################################
    if not 'time.txt' in files:
        print("ERROR: no time.txt in: " + folder)
        exit(1)

    rawTimeData = []
    with open(path + "time.txt","r") as ifh:
        # chomps \n's
        rawTimeData = [ x.split('\n')[0] for x in ifh.readlines() ]

    minutes = ""
    seconds = ""
    for time in rawTimeData:
        if 'real' in time:
            parts = time.split('\t')
            minutes = parts[1].split('m')[0]
            seconds = parts[1].split('m')[1].split('s')[0]

    # GET INSTRUCTION DATA ##########################################
    callgrindFiles = []
    for i in files:
        if i.split('.')[0] == 'callgrind':
            callgrindFiles.append(i)
    
    if len(callgrindFiles) == 0:
        print("ERROR: no callgrind files in: " + folder)
        exit(1)

    instructionCounts = []
    for i in callgrindFiles:
        with open( path + i ,"r") as ifh:
            instructionCounts.append(ifh.readlines()[-1].split(" ")[1].split('\n')[0])

    count = sum([ int(x) for x in instructionCounts])/len(instructionCounts)

    # Store data
    data[folder] = dict()
    data[folder]['minutes'] = minutes
    data[folder]['seconds'] = seconds
    data[folder]['instructions'] = count

for benchmark in data:
    for header in metaData:
        print("%s\t%s\t%s" % (benchmark,header,data[benchmark][header]))



