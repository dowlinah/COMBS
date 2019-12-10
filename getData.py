#!/usr/bin/python3

import os


base = 'benchmarks'
folders = [ 
        "2d-heat"
        ,"FourierBenchmarks"
        ,"advection-diffusion"
        ,"fidibench"
        ,"hpcg" 
        ,"ks-pde"
        ,"lid-driven-cavity"
        ,"matrix-mpi"
        ,"monte-carlo"
        ,"phase-retrieval-benchmarks"
        ,"radix_sort"
        ,"sombrero"
        ]

metaData = ['minutes','seconds','instructions','memory']
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
    for time in rawTimeData[-3:]:
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

    # GET INSTRUCTION DATA ##########################################
    if not 'mem.txt' in files:
        print("ERROR: no mem.txt in: " + folder)
        exit(1)

    memRaw = []
    with open(path+'mem.txt','r') as ifh:
        memRaw = ifh.readlines()

    # Store data
    data[folder] = dict()
    data[folder]['minutes'] = minutes
    data[folder]['seconds'] = seconds
    data[folder]['instructions'] = count
    data[folder]['memory'] = memRaw[0].split('\n')[0]

for benchmark in data:
    for header in metaData:
        print("%s\t%s\t%s" % (benchmark,header,data[benchmark][header]))

# output latex table
headers = ["Benchmark", "Time(s)", "Instruction Count", "Max Memory Usage(B)"]

with open('output.tex','w') as ofh:
    ofh.write("\\begin{tabular}{|c|c|c|c|}\n\\hline\n")
    ofh.write("\t%s \\\\ \hline \n" % (" & ".join(headers) ) )
    for benchmark in data:
        finalSeconds = (float(
            data[benchmark]['minutes'])*60
            )+float(data[benchmark]['seconds'])
        
        ofh.write(
            "\t%s & %s & %s & %s \\\\ \\hline \n" % (
                benchmark,
                finalSeconds,
                data[benchmark]['instructions'],
                data[benchmark]['memory']
            )
        )

    ofh.write("\\end{tabular}\n")

