#Rejection Sampling Monte Carlo
#David Makay on Rejection Sampling called :Information Theory Inference and Learning Algorithms

import numpy as np
import random
import argparse
import sys
import pickle


def generate_branching_time(time_inf_dict):
    total = 0
    keys = list(time_inf_dict.keys())
    p_dict = {}
    for index, k in enumerate(keys):
        if index == 0:
            total += k*time_inf_dict[k]
            p_dict[k] = k*time_inf_dict[k]
        else:
            total += (k-keys[index-1])*time_inf_dict[k]
            p_dict[k] = (k-keys[index-1])*time_inf_dict[k] + p_dict[keys[index-1]]
    for k in p_dict:
        p_dict[k] /= total

    r1 = random.random()
    for index, k in enumerate(list(p_dict.keys())):
        if r1 < p_dict[k]:
            r2 = random.uniform(list(p_dict.keys())[index-1],k)
            break
    return r2


def getArgs():
	parser = argparse.ArgumentParser(description='Perform Rejection Sampling')
	parser.add_argument('-sf', '--successful-favites-data-file', help='FAVITES output file', required=True)
	parser.add_argument('-ff', '--failed-favites-data-file-pickle', help='FAVITES output file', required=False)
	parser.add_argument('-b', '--beast-data-file', help='Processed BEAST output file', required=True)
	parser.add_argument('-d', '--pdf-date', help='Date for rejection sampling: nov17 or dec01', required=True)
	parser.add_argument('-o', '--output-file-name', help='name of output file', required=True)
	parser.add_argument('-s', '--skip-shift', help='ignore shifts in TMRCA post-Dec. 24, 2019', required=False, action='store_true')
	parser.add_argument('-a', '--include-unascertained', help='count ascertained and unascertained cases', required=False, action='store_true')
	parser.add_argument('-if', '--include_failed', help='include failed simulation', required=False, action='store_true')
	args = parser.parse_args()
	return args

arg = getArgs()

successfulFavitesPath = arg.successful_favites_data_file
failedFavitesPath = arg.failed_favites_data_file_pickle
beastPath = arg.beast_data_file
dec01 = 2019+((335-1)/365)
nov17 = 2019+((321-1)/365)

if arg.pdf_date == 'dec01': 
	pdf = dec01
elif arg.pdf_date == 'nov17': 
	pdf = nov17
else:
	print("Missing date for posterior density function: 'nov17' or 'dec01'")
	sys.exit(2)
outPath = arg.output_file_name

successfulFavitesDays = []
for line in open(successfulFavitesPath):
	if 'Sample' not in line:
		p = line.rstrip().split('\t')
		sample = p[0]
		coal_time = float(p[1])/365
		if arg.include_unascertained == True: 
			first_case = min(float(p[2]),float(p[3]))/365
		else: 
			first_case = float(p[2])/365
		successfulFavitesDays.append([sample, coal_time, first_case])

beastTMRCA = []
for line in open(beastPath):
	if 'TMRCA' not in line:
		if arg.skip_shift == False: 
			beastTMRCA.append([line.split('\t')[0],float(line.split('\t')[1]) + float(line.split('\t')[2])])
		else: 
			beastTMRCA.append([line.split('\t')[0],float(line.split('\t')[1])])

posterior = []
if arg.include_failed == False:
	for tmrca in beastTMRCA:
		IC = 2020.0
		#print (tmrca)
		while IC > pdf:
			favites = random.choice(successfulFavitesDays)
			sample = favites[0]
			timeToCoal = favites[1]
			timeToCase = favites[2]
			# date = TMRCA - (DaysToCoalescent-DaysToSymptoms)
			dateOfSymptoms = tmrca[1] - timeToCoal + timeToCase
			dateOfInfection = tmrca[1] - timeToCoal
			#print (tmrca, IC, d)
			outList = [tmrca[0],tmrca[1],sample,timeToCoal,timeToCase,dateOfInfection]
			if dateOfSymptoms < pdf:
				posterior.append(outList)
				IC = dateOfSymptoms
				#print(IC)

	outFile = open(outPath,'w')
	outFile.write('BEASTGeneration\tBEASTtmrca\tFavitesSim\tTimeToCoal\tTimeToCase\tDateOfIndexCase\tIndexAtCoalescent\n')
	for AB in posterior: 
		if AB[1] == AB[5]:
			outFile.write(AB[0]+'\t'+str(AB[1])+'\t'+str(AB[2])+'\t'+str(AB[3])+'\t'+str(AB[4])+'\t'+str(AB[5])+'\t1\n')
		else:
			outFile.write(AB[0]+'\t'+str(AB[1])+'\t'+str(AB[2])+'\t'+str(AB[3])+'\t'+str(AB[4])+'\t'+str(AB[5])+'\t0\n')
	outFile.close()

else:
	with open(failedFavitesPath, 'rb') as handle:
		failedFavitesDays = pickle.load(handle)

	# failedFavitesDays = json.load(open(failedFavitesPath))				
	for tmrca in beastTMRCA:
		IC = 2020.0
		#print (tmrca)
		while IC > pdf:
			# print(list(failedFavitesDays.keys()))
			successfulFavites = random.choice(successfulFavitesDays)
			sucFavites_sample = successfulFavites[0]
			sucFavites_timeToCoal = successfulFavites[1]
			sucFavites_case = successfulFavites[2]

			faiFavites_sample = random.choice(list(failedFavitesDays.keys()))
			if arg.include_unascertained == True:
				faiFavites_firstCase = failedFavitesDays[faiFavites_sample]['First Case']
			else:
				faiFavites_firstCase = failedFavitesDays[faiFavites_sample]['First Ascertained']
				while faiFavites_firstCase == -1:
					faiFavites_sample = random.choice(list(failedFavitesDays.keys()))
					faiFavites_firstCase = failedFavitesDays[faiFavites_sample]['First Ascertained']
			faiFavites_end = failedFavitesDays[faiFavites_sample]['End']
			faiFavites_timeInfDict = failedFavitesDays[faiFavites_sample]['TimeInfDict']
			branching_time = generate_branching_time(faiFavites_timeInfDict)

			# date = TMRCA - (DaysToCoalescent-DaysToSymptoms)
			timeToCase = min(faiFavites_firstCase, sucFavites_case)
			dateOfSymptoms = tmrca[1] - sucFavites_timeToCoal - branching_time + timeToCase
			dateOfFirstInfection = tmrca[1] - sucFavites_timeToCoal - branching_time
			dateOfSuccessfulEpidemic = tmrca[1] - sucFavites_timeToCoal
			dateOfFailedEpidemicEnd = tmrca[1] - sucFavites_timeToCoal - branching_time + faiFavites_end

			if dateOfFailedEpidemicEnd > tmrca[1]:
				timeAtFaiFavites = tmrca[1] - dateOfFirstInfection
				for index, k in enumerate(list(faiFavites_timeInfDict.keys())):
					if timeAtFaiFavites < k:
						faiFavitesInfected = faiFavites_timeInfDict[list(faiFavites_timeInfDict.keys())[index-1]]
						break
			else:
				faiFavitesInfected = 0

			#print (tmrca, IC, d)
			outList = [tmrca[0],tmrca[1],sucFavites_sample,faiFavites_sample,dateOfFirstInfection,dateOfSymptoms,dateOfSuccessfulEpidemic,dateOfFailedEpidemicEnd,branching_time,timeToCase,sucFavites_timeToCoal,faiFavites_end,faiFavitesInfected]
			if dateOfSymptoms < pdf:
				posterior.append(outList)
				# print(outList[-1])
				IC = dateOfSymptoms
				#print(IC)

	outFile = open(outPath,'w')
	# FAVITES sim = successful_failed if failed included; else just successful
	outFile.write('BEASTGeneration\tBEASTtmrca\tSuccessfulFavitesSim\tFailedFavitesSim\tDateOfFirstInfection\tDateOfSymptoms\tDateOfSuccessfulEpidemic\tDateOfFailedEpidemicEnd\tBranchingTime\tTimeToCase\tTimeToCoal\tTimeToFailedEnd\tFailedInfectionsAtTMRCA\n')
	for AB in posterior: 
		outFile.write(AB[0]+'\t'+str(AB[1])+'\t'+str(AB[2])+'\t'+str(AB[3])+'\t'+str(AB[4])+'\t'+str(AB[5])+'\t'+str(AB[6])+'\t'+str(AB[7])+'\t'+str(AB[8])+'\t'+str(AB[9])+'\t'+str(AB[10])+'\t'+str(AB[11])+'\t'+str(AB[12])+'\n')
	outFile.close()

"""
tIC = []
for AB in posterior: tIC.append(AB[2])

(np.mean(tIC)-2019)*365-1
(np.percentile(tIC,5)-2019)*365-1
(np.percentile(tIC,1)-2019)*365-1
"""