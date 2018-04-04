#!/usr/bin/env python
# -*-coding:Utf-8 -*


############# Imports 
import os 
import time
from Bio.Seq import Seq
import subprocess
import sys
import pickle
import os.path


### Log of the execution: (need to verify if a crash of LoRTE, does indeed give a usable log file--> opposition of an empty file)
path_main=str(os.getcwd())
LogName='/'+time.strftime('%d_%m_%Y_LogLoRTE_%H_%M_%S')

#Loading parametres:
if len(sys.argv)>1:#If user did enter something after 'python LoRTEx.py ...'
 print ("Specific path for parameter file entered")
 param_path=str(sys.argv[1])
 if os.path.isfile(param_path):#Check if it exists. if true create log file there.
  log=open(sys.argv[1][0:param_path.rfind('/')]+LogName,'w')
  log.close()
  log=open(sys.argv[1][0:param_path.rfind('/')]+LogName,'a')
  log.write('Program is starting'+"\n")
  print("File %s found, loading paramettres."%param_path)
  log.write("File %s found, loading paramettres."%param_path)
 else:
  print("Error %s doesn't exist/coudn't find it."%param_path)
  exit()
 
else:#User didn't enter any specific param file: Search for local one.
 print('No specific path entered for parameter file, searching for one in the curent location')
 param_path=str(os.getcwd()+"/LoRTE-Parameters")
 if os.path.isfile(param_path):#Check if it exists.
  log=open(path_main+LogName,'w')
  log.close()
  log=open(path_main+LogName,'a')
  log.write('Program is starting'+"\n")
  print("File %s found, loading paramettres."%param_path)
  log.write("File %s found, loading paramettres."%param_path)
 else:
  print("Error %s doesn't exist/coudn't find it."%param_path)
  exit()

with open(param_path) as param:
 param=param.read().split('\n')

param_spe=list()
for x in param:
 if len(x)>0:
  param_spe.append(x)#No empty lines


#Localisation of reference genome:
path_ref_genome=param_spe[1]
#Localisation of the Transposable Element identified on the reference genome:
path_TE_annotated=param_spe[3]
#Localisation of Long read file (Pac bio):
path_pacbio=param_spe[5]
#Consensus of transposable elements:
path_consensus_TE=param_spe[7]

##Name of output files: Need to be sure the file doesn't exist yet...
#Need to be sure of the input format /at the end or not
name_folder_results=param_spe[9]
if name_folder_results[-1]=="/":
	name_folder_results=name_folder_results[0:-1]
name_file_flank=name_folder_results+'/Out'
name_file_flank_step2=name_file_flank+'2'


#EValue for all aligments (Blastn/megablast):
e_value=param_spe[11]
e_value=str(e_value)
#Genome equivalent (reads total length/ reference genome length):
sequencing_dept=float(param_spe[13])
#Ordre of the elements in the blast input you entered:
ordre_output=param_spe[15]
ordre_output=ordre_output.split(",")
#Maximum length beetween two alligned flanking sequence to consider possible TE:
maximum_length_TE=param_spe[17]
maximum_length_TE=int(maximum_length_TE)
#Length of flanking sequence that will be exctracted:
length_flank_seq=param_spe[19]
length_flank_seq=int(length_flank_seq)
#Number of cores blast can use to run:
number_cores=param_spe[21]
number_cores=str(number_cores)


#Creating results folder.
os.system('mkdir -p '+name_folder_results) 


#Checking if all the files of the input are correct: No need to start if the user didn't enter the corect path to files...
if os.path.isfile(path_ref_genome) and os.path.isfile(path_TE_annotated) and os.path.isfile(path_pacbio) and os.path.isfile(path_consensus_TE) and os.path.isdir(name_folder_results):
	log.write('All the input files exist !')
else:
	log.write('Error in the files path (file not found/output file) ! Program is exiting, check files and retry.')
	print('Error in the files path ! Program is exiting, check files and retry.')
	exit()


### Launch of script for ram monitoring, and total time analysis. Usefull for unexpected crash (ram limit).Subprocess needs to be stopped at the end.
#ram_monitoring=subprocess.Popen("./RamMonitor.sh")









#
##
###
####All the functions used (def):




def sorterTE(blastOutput,TE_name,chromosome_name,Start,Stop,senss,collumn_separation,first_blast_line):
	"""
	Put all different type of blast output in the needed order. (Changes collumns organisation and ads +/- if not added.)
	Ouput has 'Unique ID'_'TE Name'\t'Contig of TE'\t'Start'\t'Stop'\t'Sens'
	"""
	#One last visual verification for the clumsy user...
	print('Collumn nb for  TE name: ',TE_name)
	print('Collumn nb for  genomeRefs contig: ',chromosome_name)
	print('Collumn nb for  Sstart: ',Start)
	print('Collumn nb for  Ssend: ',Stop)
	print('Collumn nb for  Sens/antisens: ',senss)
	print('Starting at the line Nb: ',first_blast_line)
	
	collumn_separation=str(collumn_separation)
	if collumn_separation.count('t')>0:
		collumn_separation='\t'
	else:
		print("Blast input is separated by: ",collumn_separation)#For testing purposes
	with open(blastOutput, 'r') as blast_lines:
		blast_lines=blast_lines.read().split("\n")

	ordered_output=list()

	i=int(first_blast_line)
	while i<((len(blast_lines))-1):
		seperated_line_blast=blast_lines[i].split(collumn_separation)
		ordered_line=list()
		ordered_line_final=list()
		a=(str(i)) 
		ordered_line.append(a)#Add the unique ID
		element_retirer_TE_name=[".fa",".contig","--","__"," "]#The usual suspects
		final_name=(seperated_line_blast[TE_name])
		for x in element_retirer_TE_name:
			if (final_name.count(x))>0:
				final_name=(seperated_line_blast[TE_name].replace(x,""))
		

		ordered_line.append(final_name)
		#Add a _ to differentiate 10xxx from 110xxx in futur search of the prog
		intermed="_" 
		intermed+="_".join(ordered_line)
		ordered_line_final.append(intermed)
		ordered_line_final.append(seperated_line_blast[chromosome_name]) #Ajout de la localisation: nom contig OU READ
	
		#Determining the direction:
		senstr=str(senss)
		#If length for senstr is above 1 that mean it isn't a collumn id '4','3' etc. We have to calculate the +/- with the different blast coordinates.
		if len(senstr)>2:
			if ((int(seperated_line_blast[int(Start)]))-(int(seperated_line_blast[int(Stop)])))<0:
				sens="+"
			else:
				sens="-"
			#Adding start and stop.
			ordered_line_final.append(seperated_line_blast[Start])
			ordered_line_final.append(seperated_line_blast[Stop])
			
		else: #+/- is defined we just have to retrieve it.
			sens=seperated_line_blast[int(senss)]
			
			if sens=='+': 
				ordered_line_final.append(seperated_line_blast[Start])
				ordered_line_final.append(seperated_line_blast[Stop])
			else:#If +/- is defined the blast output is given with rising number for start and stop: Need to invert it for - alignements.
				ordered_line_final.append(seperated_line_blast[Stop])
				ordered_line_final.append(seperated_line_blast[Start])
			
		ordered_line_final.append(sens)
		ordered_output.append("\t".join(ordered_line_final))
		i+=1	
	#Saving data	
	filename_Final=blastOutput+"-Sorted-LorTE"
	monFichierDeux=open(filename_Final, "w")
	print("The sorted TE list has been created")
	for item in ordered_output:
		monFichierDeux.write("%s\n" % item)	
	monFichierDeux.close()
	return(filename_Final)





def cleaningFlank(blast,flank_file,non_resolvable):
	"""
	Takes out all the flanking sequence that alligned more than once on the reference genome. 
	Marks them as unresolvable for the final results.
	"""
	with open(flank_file) as f1:
		f1=f1.read()
		f1n=f1.split("\n")

	with open(blast) as blast:
		blast=blast.read()
	unresolvable_list=str()
	i=0
	length=len(f1n)-1
	while i<(length):
		line=f1n[i]
		TE_name=line[1:len(line)]
		if blast.count(TE_name)>2: #2==> One in the description of blast output plus one for each alignment.
			start=f1.find(TE_name)-1 #includes the '>'
			end_ID=f1.find('\n',start+1)
			end_seq=f1.find('\n',end_ID+1)+1
			old=f1[start:end_seq]
			#Deletion of element from the flanking sequence.
			f1=f1.replace(old,"")
			TE_name_final=old[old.find("__")+1:old.rfind('__')]
			unresolvable_list+=(TE_name_final+"\n")
		i+=2
	#Saving datas
	unresolvable_file=open(non_resolvable,"a")
	unresolvable_file.write(unresolvable_list)
	unresolvable_file.close()

	flanking_clear=open(flank_file,'w')
	flanking_clear.write(f1)
	flanking_clear.close()

	return(non_resolvable)






def cleaningFlankMegablast(blast,flank_file,non_Resolvable):
	"""
	Flags and takes out the megablast aligment who appears multiple time on the references genome : 
	Their ID is stored in the 'unresolved' file
	They are removed from the megablast files (no extraction of their sequences etc)
	"""
	blast_name_file=blast
	with open(flank_file) as f1:
		f1n=f1.read().split("\n")

	with open(blast) as blast:
		blast=blast.read()
	blast_lines=blast.split('\n')
	blast_header=blast_lines.pop(0)

	unresolvable_list=list()
	i=0
	listeSupr=list()
	taille=len(f1n)-1 #file split('\n') is terminated by an empty entry
	while i<(taille):
		ligne=f1n[i]
		TE_name=ligne[1:len(ligne)]#Ignores the '>' of fasta ID:
		if blast.count(TE_name)>1: #If more than 1 alignment in the reference genome it is unresolvable
			for x in blast_lines:
				if TE_name in x:
					listeSupr.append(x)
			unresolvable_list.append(TE_name[1:(TE_name.rfind('__'))])
		i+=2


	#Deleting all the unresolvable
	blast_lines=list(set(blast_lines)-set(listeSupr))
	megablast_cleaned = open (blast_name_file,"w")
	megablast_cleaned.write(blast_header+'\n')
	for x in blast_lines:
		if len(x)>1:
			megablast_cleaned.write(x+'\n')
	megablast_cleaned.close()
	
	#Saving Unresolvable list. List generated externaly, we only append to it.
	unresolvable_list=list(set(unresolvable_list))
	unresolvable_list_file=open(non_Resolvable,"a")
	for x in unresolvable_list:
		if len(x)>1:
			unresolvable_list_file.write(x+'\n')
	unresolvable_list_file.close()
	return(blast_name_file)






def MonoLineFasta(filename):
	"""
	Put Fasta format un one line sequnce mod.
	Filename has to be a STR.
	Uses the ">" as reference for the end of the sequence n-1 and n.
	"""
	strfinal=str()
	with open(filename) as f:
		for line in f:
			if line.find(">")==0:
				line=line.replace(">","\n>")
			else:
				line=line.replace('\n','')
			strfinal+=line
	print('Saving Data')
	output_file=open(filename,"w")

	output_file.write(strfinal[1:len(strfinal)])#0 is a \n.
	#close the file
	output_file.close()
	return(filename)








def makedb(fichier):
	"""
	Creats a db file: nucl at the file location. Needed for alignments.
	"""
	command_file=(" -in "+fichier+" -dbtype nucl")
	os.system("makeblastdb"+command_file)
	return("Fichier DB crée")





def create_flanq_file(te_anoo, genomee,Output_filename):
	"""
	Creates 5' and 3' flanking sequences file.
	Takes TE annotated ends and extractes the matching flanking sequence.
	Reverse the sequence if needed (- aligment)
	length_flank_seq is the input entered in the param: length in nt.
	Has a 50 nt limit, or doesn't take the flanking sequence (if te is on the edge of a read it wont be in the file.
	"""
	print('Extraction of the flanking sequences')
	#File sep \t organised with 0TE-name 1-contig 2-Start 3-Stop 4-(+/-)Alignment
	with open(te_anoo,'r') as te_note:
		te_ano=te_note.read().split("\n")

	#Reference genome:
	with open(genomee,'r') as genome:
		genome=genome.read().split("\n")

	list_three_p=list()
	list_five_p=list()
	i=0
	lenght_TE_file=len(te_ano)-1
	while i<(lenght_TE_file):
		separated_line=te_ano[i].split("\t")#Separation with \t: 0-TE name 1-contig 2-Start 3-Stop 4-Sens(+/-)
		#Taking the sequence in the fasta
		I=0
		while I<(len(genome)-1):
			if genome[I].count(separated_line[1])>0: #If fasta name is equal to contig of TE:
				Bound_five_p=int(separated_line[2]) #Start 5'
				Bound_three_p=int(separated_line[3]) #Stop 3'
				if separated_line[4].count("+")>0:#  If + alignement
					#For 5'
					if Bound_five_p<length_flank_seq:
						Other_bound_five=0
					else:
						Other_bound_five=(Bound_five_p-length_flank_seq)
					#For 3'	
					if Bound_three_p+length_flank_seq>(len(genome[I+1])-1):
						Other_bound_three=(len(genome[I+1])-1)
					else:
						Other_bound_three=(Bound_three_p+length_flank_seq)

					#Exctracting sequnce +:
					#5'
					extracted_sequence_three_p=(genome[I+1][Bound_three_p:Other_bound_three])
					#3'
					extracted_sequence_five_p=(genome[I+1][Other_bound_five:Bound_five_p])
					
				#If - alignment: REVERSE complement !
				else:
					#5'
					if Bound_five_p+length_flank_seq>(len(genome[I+1])-1):
						Other_bound_five=(len(genome[I+1]))
					else:
						Other_bound_five=(Bound_five_p+length_flank_seq)
					#3'
					if Bound_three_p<length_flank_seq:
						Other_bound_three=0
					else:
						Other_bound_three=(Bound_three_p-length_flank_seq)
					
					#Extracting sequences (reverse complemented):
					#5'
					extracted_sequence_five_p=Seq(genome[I+1][Bound_five_p:Other_bound_five])
					extracted_sequence_five_p=extracted_sequence_five_p.reverse_complement()
					extracted_sequence_five_p=str(extracted_sequence_five_p)
					#3'
					extracted_sequence_three_p=Seq(genome[I+1][Other_bound_three:Bound_three_p])
					extracted_sequence_three_p=extracted_sequence_three_p.reverse_complement()
					extracted_sequence_three_p=str(extracted_sequence_three_p)
					
				

				#If their length is above 50nt only they are stored in the file !! Most used for step two when the sequence are exctracted from reads.
				if len(extracted_sequence_three_p)>50:
					header_three_p=[">",separated_line[0]+"_",separated_line[1],separated_line[4]]#>-NomTE -contig-sens
					str_header3prim="_".join(header_three_p)
					list_three_p.append(str_header3prim)
					list_three_p.append(extracted_sequence_three_p)
					
				if len(extracted_sequence_five_p)>50:
					header_five_p=[">",separated_line[0]+"_",separated_line[1],separated_line[4]]
					str_header5prim="_".join(header_five_p)
					list_five_p.append(str_header5prim)
					list_five_p.append(extracted_sequence_five_p)
				
				I+=1
			else:
				I+=1
		i+=1
		advancement=(i*100)/lenght_TE_file
		sys.stdout.write("\r%d%%" % advancement)
		sys.stdout.flush()
	
	#Creating 3' file
	file_name_3prim=(Output_filename+"-3prim")
	first_file_flank = open(file_name_3prim, "w")
	for item in list_three_p:
		first_file_flank.write("%s\n" % item)

	#Creating 5' file	
	file_name_5prim=(Output_filename+"-5prim")
	second_file_flank = open(file_name_5prim, "w")
	for wtf in list_five_p:
		second_file_flank.write("%s\n" % wtf)
	second_file_flank.close()
	first_file_flank.close()
	print("")
	return(file_name_3prim,file_name_5prim)










def FlankingAlignmentAnalyser(Blast_file,blastTwo_file,genomePacBio_file):
	"""
	Exctracts the sequence in between the flanking sequence alignment:
	If the two flanking aligned on the same contig/read: 5'NNNNNNN3'
	If only one of the two has aligned: 5'NNNNNN+2000  -2000NNNNN3' etc if - alignment.
	"""
	list_IDflank_same_contig=list()#Contains list of IDs that matched
	list_sequence_fiveANDthree=list()#Contains the sequence in between 5' and 5' alignement
	list_sequences_3prim_only=list()#  -2000NNNNN3' 3' only
	list_sequences_5prim_only=list()#  5'NNNNNN+2000   only
	with open(Blast_file,'r') as blastUN:
		BlastLINES_fileONE=blastUN.read().split("\n")

	with open(blastTwo_file,'r') as BlastTWO:
		BlastLINES_fileTWO=BlastTWO.read().split("\n")

	with open(genomePacBio_file) as fichierPB:
		lignePB=fichierPB.read().split('\n')
	#Creation the dict that will have the ID read in key and the sequence in value:
	PacBio_dictio=dict()
	for i,x in enumerate(lignePB):
		if len(x)>0: 
			if x.startswith(">"):
				s=x.find(" ")
				if s!=-1:#Blast output doesnt take the rest of Read ID if it contains a space bar in it. Everything is covered in this prog :) ! Or we sure hope so !
					PacBio_dictio[x[0:s]]=lignePB[i+1]
				else:
					PacBio_dictio[x]=lignePB[i+1]
	del(lignePB)#Doesn't change one bit on the memory usage but still. 
	Blast_file=Blast_file.replace("-megablast","")
	Blast_file=Blast_file.replace("-3prim","")
	
	
	print("Searching 5' 3' to pool together:")
	#Search stas from 3' file and search candidates in the 5' file.
	n=1 #Megablast output contains 1 lines of description.. 
	length_blastfileONE=len(BlastLINES_fileONE)
	while n<(length_blastfileONE-1):
		Blast_SeparatedONE=BlastLINES_fileONE[n].split("\t")
		I=1 #Same output file format.
		while I<(len(BlastLINES_fileTWO)-1):
			if ((BlastLINES_fileTWO[I].find(Blast_SeparatedONE[0]))>-1): #If the name is identical it also inclues the sens.
				if ((BlastLINES_fileTWO[I].find(Blast_SeparatedONE[1]))>-1): #If it's on the same contig/read
					Blast_SeparatedTWO=BlastLINES_fileTWO[I].split("\t")  
					first_number=max(int(Blast_SeparatedONE[8]),int(Blast_SeparatedONE[9]))#Le max de la colone3'
					if int(Blast_SeparatedONE[8])>int(Blast_SeparatedONE[9]):#Start>End ==> - alignment
						Blast_SensONE="-"
					else:
						Blast_SensONE="+"
					if int(Blast_SeparatedTWO[8])>int(Blast_SeparatedTWO[9]):
						Blast_SensTWO="-"
					else:
						Blast_SensTWO="+"
					second_number=min(int(Blast_SeparatedTWO[8]),int(Blast_SeparatedTWO[9]))
					if (abs(first_number-second_number))<maximum_length_TE and Blast_SensONE==Blast_SensTWO:
						list_IDflank_same_contig.append(BlastLINES_fileONE[n]) 
						list_IDflank_same_contig.append(BlastLINES_fileTWO[I])
						#Sequence:
						sequence_start_end=PacBio_dictio[">"+Blast_SeparatedONE[1]]
						first_number=min(int(Blast_SeparatedTWO[9]),int(Blast_SeparatedONE[8]))
						second_number=max(int(Blast_SeparatedONE[8]),int(Blast_SeparatedTWO[9]))
						a=str()
						a=(">"+Blast_SeparatedONE[1]+Blast_SeparatedONE[0]+"_/_"+str(Blast_SensTWO)+"_/_"+str(first_number)+"_/_"+str(second_number)) #New ID for the sequence
						list_sequence_fiveANDthree.append(a)
						#Rev complement if needed
						if Blast_SensTWO=="-":
							futurRev=Seq(sequence_start_end[first_number:second_number])
							futurRev=futurRev.reverse_complement()
							futurRev=str(futurRev)
							list_sequence_fiveANDthree.append(futurRev)
						else:#No need to rev compl
							list_sequence_fiveANDthree.append(sequence_start_end[first_number:second_number])

						I+=1 
					else:
						I+=1
				else:
					I+=1
			else:
				I+=1

		n+=1 #+1 first megablas and go for an other round !! 
		ii=(n*100)/length_blastfileONE
		sys.stdout.write("\r%d%%" % ii)
		sys.stdout.flush()
	#Saving data:
	lesflankreunit = open(Blast_file+"-Read-53both-IDs", "w")
	for wtf in list_IDflank_same_contig:
		lesflankreunit.write("%s\n" % wtf)
	lesflankreunit.close()
	
	lesseqreads = open(Blast_file+"-53-extraction", "w")
	for hum in list_sequence_fiveANDthree:
		lesseqreads.write("%s\n" % hum)
	lesseqreads.close()

	#Deleting all uneccessy data stored on ram:
	del(list_IDflank_same_contig)
	del(list_sequence_fiveANDthree)
	
	
	print("\nGetting all the 5' only")
	###Taking all the 5' only from the previous list created:
	with open(Blast_file+"-Read-53both-IDs") as lesdeux:
		lesdeux=lesdeux.read()
	specific_line=1#file has a header.
	while specific_line<len(BlastLINES_fileTWO)-1:
		if BlastLINES_fileTWO[specific_line] not in lesdeux:
			list_sequences_5prim_only.append(BlastLINES_fileTWO[specific_line])
			specific_line+=1
		else:
			specific_line+=1
	FiveUniqueFile = open(Blast_file+"-megablast-5prim-only", "w")
	for Fiveunique in list_sequences_5prim_only:
		FiveUniqueFile.write("%s\n" % Fiveunique)
	FiveUniqueFile.close()
	
	print("Getting all the 3' only")
	###Same here
	specific_line=1#file has a header.
	while specific_line<len(BlastLINES_fileONE)-1:
		if BlastLINES_fileONE[specific_line] not in lesdeux:
			list_sequences_3prim_only.append(BlastLINES_fileONE[specific_line])
			specific_line+=1
		else:
			specific_line+=1	
	ThreeUniqueFile = open(Blast_file+"-megablast-3prim-only", "w")
	for threeUnique in list_sequences_3prim_only:
		ThreeUniqueFile.write("%s\n" % threeUnique)
	ThreeUniqueFile.close()
	
	
	#Getting the sequences: 5' only 3' only

	print("Extracting 5' only sequences ")
	list_sequence_extraction_5prim_unique=list()
	for x in list_sequences_5prim_only:	
		xx=x.split("\t")
		megablastStart=int(xx[8])
		megablastStop=int(xx[9])
		if (megablastStart-megablastStop)<0: #Sens "+" 0s.start - Xs.end < 0 
			sensmegablast="+"
		else:
			sensmegablast="-"
		if sensmegablast=="+": #5' NNNNNNNNn
			second_number=megablastStop+2000	
			sequence_start_end=PacBio_dictio[">"+xx[1]]			
			if len(sequence_start_end)<second_number:
				second_number=len(sequence_start_end)
			a=str()
			a=(">"+xx[1]+xx[0]+"_/_"+str(sensmegablast)+"_/_"+xx[9])
			list_sequence_extraction_5prim_unique.append(a)
			list_sequence_extraction_5prim_unique.append(sequence_start_end[megablastStop:second_number])

		else: #sens megablast -    nNNNNN 5' -> need rev compl
			second_number=megablastStop-2000
			sequence_start_end=PacBio_dictio[">"+xx[1]]
			if second_number<0:
				second_number=0
			a=str()
			a=(">"+xx[1]+xx[0]+"_/_"+str(sensmegablast)+"_/_"+xx[9])
			list_sequence_extraction_5prim_unique.append(a)
			a=Seq(sequence_start_end[second_number:megablastStop])
			a=a.reverse_complement()
			a=str(a)
			list_sequence_extraction_5prim_unique.append(a)

					

					
	seqfiveunique = open(Blast_file+"-Sequence-5prim-only", "w")
	for seqs in list_sequence_extraction_5prim_unique:
		seqfiveunique.write("%s\n" % seqs)
	seqfiveunique.close()
	
	del(list_sequence_extraction_5prim_unique)
	###Time for 3'
	print("Extracting 3' only sequences")
	list_sequence_extraction_3prim_unique=list()
	for x in list_sequences_3prim_only:
		xx=x.split("\t")
		megablastStart=int(xx[8])
		megablastStop=int(xx[9])
		if (megablastStart-megablastStop)<0: #Sens + query start  -  query end
			sensmegablast="+"
		else:
			sensmegablast="-"
		if sensmegablast=="+": #nNNNNNNNN 3'
			second_number=megablastStart-2000
			sequence_start_end=PacBio_dictio[">"+xx[1]]
			if second_number<0:
				second_number=0
			a=str()
			a=(">"+xx[1]+xx[0]+"_/_"+str(sensmegablast)+"_/_"+xx[8])
			list_sequence_extraction_3prim_unique.append(a)
			list_sequence_extraction_3prim_unique.append(sequence_start_end[second_number:megablastStart])


		else:
			second_number=megablastStart+2000
			sequence_start_end=PacBio_dictio[">"+xx[1]]
			if second_number>(len(sequence_start_end)-1):
				second_number=len(sequence_start_end)
			a=str()
			a=(">"+xx[1]+xx[0]+"_/_"+str(sensmegablast)+"_/_"+xx[8])
			list_sequence_extraction_3prim_unique.append(a)
			a=Seq(sequence_start_end[megablastStart:second_number])
			a=a.reverse_complement()
			a=str(a)
			list_sequence_extraction_3prim_unique.append(a)


	seqthreeUnique = open(Blast_file+"-Sequence-3prim-only", "w")
	for seqs in list_sequence_extraction_3prim_unique:
		seqthreeUnique.write("%s\n" % seqs)
	seqthreeUnique.close()
	
	return(Blast_file+"-Read-53both-IDs",Blast_file+"-53-extraction",Blast_file+"-Sequence-3prim-only",Blast_file+"-Sequence-5prim-only")





def SequenceCleaner(input_file,input_type):
	"""
	Cleans multiple megablast aligment on a contig. (even though it should never happen in ideal condition. It happens...)
	Takes the further possible bound when needed (multiple alignments are found.)
	file_type is the three possible type of input file : 5' only 3' only and 5'-3'
	"""
	with open(input_file) as john:
		fasta=john.read().split("\n")
	length_fasta=len(fasta)
	list_alreadyDone=list()
	list_to_Save=list()
	#5' Only file
	if str(input_type)=="five":
		print("File 5' being processed...")
		x=0
		while x<(length_fasta-1):
			fastaXsep=fasta[x].split("_/_")
			fastaXid=str(fastaXsep[0]) #Contains TE name and contig/read ID
			fastaXsense=str(fastaXsep[1]) # - or +
			memory_position=x 
			found=0
			if any(fastaXid in s for s in list_alreadyDone):#Skip if ID has already been treated.
				x+=2
			else:
				n=x+2 #Going to next fasta ID
				if n>(length_fasta-1):#When we are at the last one.
					list_to_Save.append(fasta[x])
				else:
					while n<(length_fasta-1):
						fastaNsep=fasta[n].split("_/_")
						fastaNid=str(fastaNsep[0]) #Name of contig/read and TE
						fastaNsense=str(fastaNsep[1])
						if found==0:
							if fasta[n].find(fastaXid)>-1 and fastaNsense==fastaXsense: #Same ID and sense
								#Compares relative position on the contig/read.
								if fastaNsense=="+": #The closer to 0 if multiple alignment.
									if (int(fastaXsep[2])<int(fastaNsep[2])):
										memory_position=x
										MemoryRead_location=(int(fastaXsep[2]))
										found=1
									else: #If n is smaller than th firstr
										memory_position=n
										MemoryRead_location=(int(fastaNsep[2]))
										found=1			
								else: #5' megablast antisense: the closer to the end
									if (int(fastaXsep[2])>int(fastaNsep[2])):
										memory_position=x
										MemoryRead_location=(int(fastaXsep[2]))
										found=1
									else: #Second is bigger
										memory_position=n
										MemoryRead_location=(int(fastaNsep[2]))
										found=1
						else: #Already found one: comparing to MemoryRead_location
							if fasta[n].find(fastaXid)>-1 and fastaNsense==fastaXsense:
								if fastaNsense=="+":
									if (int(fastaNsep[2])<MemoryRead_location):
										memory_position=n
										MemoryRead_location=(int(fastaXsep[2]))
								else:
									if (int(fastaNsep[2])>MemoryRead_location):
										memory_position=n
										MemoryRead_location=(int(fastaNsep[2]))
										#No need for else: the memory doesn't change
						n+=2
						if n>=(length_fasta-1):
							list_alreadyDone.append(fastaXsep[0])
							list_to_Save.append(fasta[memory_position])
							list_to_Save.append(fasta[memory_position+1])
				x+=2
			ii=(x*100)/length_fasta
			sys.stdout.write("\r%d%%" % ii)
   			sys.stdout.flush()



	#3' Only file:
	elif str(input_type)=="three": 
		print("File 3' being processed...")
		x=0
		while x<(length_fasta-1):
			fastaXsep=fasta[x].split("_/_")
			fastaXid=str(fastaXsep[0])
			fastaXsense=str(fastaXsep[1])
			memory_position=x
			found=0
			if any(fastaXid in s for s in list_alreadyDone):
				x+=2
			else:
				n=x+2
				if n>(length_fasta-1):
					list_to_Save.append(fasta[x])
				else:
					while n<(length_fasta-1):
						fastaNsep=fasta[n].split("_/_")
						fastaNid=str(fastaNsep[0])
						fastaNsense=str(fastaNsep[1])
						if found==0:
							if fasta[n].find(fastaXid)>-1 and fastaNsense==fastaXsense:
								if fastaNsense=="+":#This time 3' we want the bigger int
									if (int(fastaXsep[2])>int(fastaNsep[2])):
										memory_position=x
										MemoryRead_location=(int(fastaXsep[2]))
										found=1
									else:
										memory_position=n
										MemoryRead_location=(int(fastaNsep[2]))
										found=1
										
								else:
									if (int(fastaXsep[2])<int(fastaNsep[2])):
										memory_position=x
										MemoryRead_location=(int(fastaXsep[2]))
										found=1
									else:
										memory_position=n
										MemoryRead_location=(int(fastaNsep[2]))
										found=1
						else:
							if fasta[n].find(fastaXid)>-1 and fastaNsense==fastaXsense:
								if fastaNsense=="+":
									if (int(fastaNsep[2])>MemoryRead_location):
										memory_position=n
										MemoryRead_location=(int(fastaXsep[2]))
										found+=1
								else:
									if (int(fastaNsep[2])<MemoryRead_location):
										memory_position=x
										MemoryRead_location=(int(fastaXsep[2]))
										found+=1
						n+=2 
						if n>=(length_fasta-1):
							list_alreadyDone.append(fastaXsep[0])
							list_to_Save.append(fasta[memory_position])
							list_to_Save.append(fasta[memory_position+1])
				x+=2
			ii=(x*100)/length_fasta
			sys.stdout.write("\r%d%%" % ii)
   			sys.stdout.flush()	
   			
   	#The most important: 5' and 3' flanking did align on contig/read.			
	elif input_type=="both":
		print("File 5'-3' being processed ...")
		list_alreadyDone=list()
		list_to_Save=list()
		x=0
		while x<(length_fasta-1):
			premierefois=0
			fastasep=fasta[x].split("_/_")
			fastaXid=str(fastasep[0])
			fastaXsense=str(fastasep[1])
			# old alert of i added... don't know why ..._/_ !!!!!!
			if any(fastaXid in s for s in list_alreadyDone):
				premierefois+=1
			else:
				memory_position=x
				tailleMemoire=len(fasta[x+1])
				n=x
				while n<(length_fasta-2):
					fastaNsep=fasta[n].split("_/_")
					fastaNid=str(fastaNsep[0])
					fastaNsense=str(fastaNsep[1])
					if fasta[n].find(fastaXid)>-1 and fastaNsense==fastaXsense:
						if len(fasta[n+1])>tailleMemoire:
							memory_position=n
							tailleMemoire=len(fasta[n+1])
					n+=2
					if n>=(length_fasta-2):
						list_alreadyDone.append(fastasep[0])
						list_to_Save.append(fasta[memory_position])
						list_to_Save.append(fasta[memory_position+1])
			x+=+2
			ii=(x*100)/length_fasta
			sys.stdout.write("\r%d%%" % ii)
   			sys.stdout.flush()

	else:
		print('Error')
		return('Error')
		
	input_filer = open(input_file+'-NoRepeats','w')
	for elem in list_to_Save:
		input_filer.write(elem+"\n")
	print(" Finished cleaning: ",input_file)

	return(input_file+"-NoRepeats")






def MEGABLAST(reads,query,eEevalue):
	"""
	For all flanking aligment on the reads or reference genome
	Blast+ alignment using megablast option. 
	Using E value entered by user.
	Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score::::qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
	"""
	resultsLocation=name_folder_results+'/'
	queryName=query[(query.rfind("/")+1):(len(query))]
	
	readscommand=" -db "+reads #Database (equivalent of subject)
	querycommand=" -query "+query 
	evaluePart=" -evalue "+eEevalue
	restOfCommand=" -num_threads "+number_cores+" -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -out "+resultsLocation+queryName+"-megablast"
	os.system("blastn -task blastn"+readscommand+querycommand+evaluePart+restOfCommand)
	
	return(query+"-megablast")




def BLAST(nameFileSequenceExtracted,consensus,eeValue):
	"""
	The one used for aligment of the extracted sequence 5'NNNNNNN3' 5'NNNNNNNn nNNNNNN3' on the consensus of transposable elements.
	Option -max_target_seqs 1 for best alignment only, we need a yes/no result only in the output.
	"""
	C1="blastn -task blastn -query "+nameFileSequenceExtracted	
	C2=" -evalue "+eeValue+" -db "+consensus+" -num_threads "+number_cores+" -outfmt \"7 sseqid sstart send length qseqid  evalue\" -max_target_seqs 1 -out "+nameFileSequenceExtracted+"-blastn"
	blast=C1+C2
	os.system(blast) 
	return(nameFileSequenceExtracted+"-blastn")







def BLAST2(nameFileSequenceExtracted,consensus,eeValue):
	"""
	Used to detect all remaining Transposable element in the PacBioReads.
	"""
	C1="blastn -task blastn -query "+nameFileSequenceExtracted	
	C2=" -evalue "+eeValue+" -db "+consensus+" -num_threads "+number_cores+" -outfmt \"7 sseqid qseqid sstart send length evalue\" -out "+nameFileSequenceExtracted+"-blastn"
	blast=C1+C2
	os.system(blast) 
	return(nameFileSequenceExtracted+"-blastn")






def SimpleBlast(reads,query):
	"""
	Used for aligment of the flanking sequence of negative results to see if they come from the same locus 
	Uses a specific e_value for this specific aligment of 1e-10.
	"""
	PathResults=name_folder_results+'/'
	nomquery=query[(query.rfind("/")+1):(len(query))]
	
	readscommand=" -db "+reads
	querycommand=" -query "+query
	SpecificEvalue=" -evalue 1e-10 " #The low e-value for aligning such short flanking sequence may need change one day.
	restCommand=" -num_threads "+number_cores+" -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -out "+PathResults+nomquery+"-blastN"
	os.system("blastn -task blastn "+readscommand+querycommand+SpecificEvalue+restCommand)
	return(PathResults+nomquery+"-blastN")
	





def CleanerBlastOutput0Hits(myFile):
	"""
	The original blast (input) file is modified and will only contain alignment hits.
	Lines with no alignments begin with #.
	We keep one line as header if user want to see the results. The results with no aligmen is added inb a seperate file "ZeroHits"
	And it is the file that the function returns. 
	"""
	with open(myFile,'r') as blastos: 
		blastoslines=blastos.read().split("\n")
	liste_TE_found=list()
	liste_TE_NOTfound=list()
	n=0
	i=0
	while i<(len(blastoslines)-1):
		if (len(blastoslines[i]))>0:	
			if blastoslines[i][0]!="#":
				if n==0:
					liste_TE_found.append(blastoslines[i-2]) #Adding header
					liste_TE_found.append(blastoslines[i])
					n+=1
					i+=1
				else:
					liste_TE_found.append(blastoslines[i])
					i+=1
			else:
				if blastoslines[i].find(" 0 hits found")>-1: # importance of the space before "0 hits found" to diferentiate from 10 ou xx0 hits found !
					liste_TE_NOTfound.append(blastoslines[i-2])
				i+=1
		else:
			i+=1
	#Saving datas
	positivAlignments = open(myFile, "w") #Saving clean positiv alignment on the same file.
	for hum in liste_TE_found:
		positivAlignments.write("%s\n" % hum)

	myFilezero=myFile+"-ZeroHit"
	zerohitsfile= open(myFile+"-ZeroHit","w")
	for thereIsNoTry in liste_TE_NOTfound:
		zerohitsfile.write("%s\n" % thereIsNoTry)

	positivAlignments.close()
	zerohitsfile.close()
	return(myFilezero)







def PacBionetoyeurBlastn(file5prim,file3prim,file53,readpb):
	"""
	Cleans Reads (pacBio or else) of all Transposable element already identified from previous steps: The TE considered as present in step one.
	For a read who has 3TE--> 3'only 53 and 5' only  NNNNN3'only-xxxxxxxx-5'NNNNNNNNNNNNN3'-XXXXXXXXXXX-5'onlyNNNNNNNNNN
	will become xxxxxxxxXXXXXXXXXXX
	Will use read name to find it in the three input blast files.
	"""
	
	with open(file5prim) as prim5:
		line5=prim5.read().split("\n")
	with open(file3prim) as prim3:
		line3=prim3.read().split("\n")
	with open(file53) as prim53:
		line53=prim53.read().split("\n")
	with open(readpb) as pb:
		lineRead=pb.read().split("\n")
	
	#La premiere boucle qui concerne le fichier pac bio... 
	a=0	
	listFinalRead=list()
	while a<(len(lineRead)-1):
		nameRead=lineRead[a] #">xxyy.../zr/zr RQ=0.157" --> "xxyy.../zr/zr" like in blast output
		if nameRead.find(" ")!=-1:
			nameRead=nameRead.replace(nameRead[(nameRead.find(" ")):(len(nameRead))],"")#here we delete it
		nameRead=nameRead.replace(">","")
		#Var re-initiated each read
		minimum_read=0
		maximum_read=len(lineRead[a+1])
		#Searching that ID is the blast files 53,5',3':
		#Start with 5'
		lineFound=list()
		b=1 #0 is header
		while b<(len(line5)-1):
			if line5[b].find(nameRead+"_")>-1:
				 lineFound.append(line5[b])
			b+=1
		#Min & max can be found depending on blast sens 5'+ maximum, 5'- minimum
		for i,x in enumerate(lineFound):
			lSplit=x.split("\t") #Initial split
			lSplit=lSplit[4].split("_/_")#NameID sens Start Stop
			if lSplit[1]=="+": #Sens + defines the maximum: we want the smallest
				if int(lSplit[2])<maximum_read:
					maximum_read=int(lSplit[2])
			else: #Sens - minimum: we want the biggest
				if int(lSplit[2])>minimum_read:
					minimum_read=int(lSplit[2])
				
		#3' file
		lineFound=list()
		b=1
		while b<(len(line3)-1):
			if line3[b].find(nameRead+"_")>-1:
				 lineFound.append(line3[b])
			b+=1
			
		for i,x in enumerate(lineFound):
			lSplit=x.split("\t")
			lSplit=lSplit[4].split("_/_")
			if lSplit[1]=="+": #3' +: min
				if int(lSplit[2])>minimum_read:
					minimum_read=int(lSplit[2])
			else: #3' -: max 
				if int(lSplit[2])<maximum_read:
					maximum_read3=int(lSplit[2])
		
		#Doing 53: will contain all the 53 found intervals 
		listCoordo=list()
		listCoordo.append([0,minimum_read])
		listCoordo.append([maximum_read,len(lineRead[a+1])])
		lineFound=list()
		b=1 
		while b<(len(line53)-1):
			if line53[b].find(nameRead+"_")>-1:
				 lineFound.append(line53[b])
			b+=1
			#Cleaning file to keep the intervals
		
		for x in (lineFound):
			lSplit=x.split("\t")
			lSplit=lSplit[4].split("_/_")#Same as before
			interval=list()
			#List containig both Start/stop
			interval.append(int(lSplit[2]))
			interval.append(int(lSplit[3]))
			#Adding data to the big list
			listCoordo.append(interval)

		n=0
		while n<len(listCoordo):
			startPos=listCoordo[n][0]
			endPos=listCoordo[n][1]
			x=0
			while x<len(listCoordo):
				otherStart=listCoordo[x][0]
				otherEnd=listCoordo[x][1]
				if ((startPos<otherStart and startPos<otherEnd) and (endPos<otherStart and endPos<otherEnd)) or ((startPos>otherStart and startPos>otherEnd) and (endPos>otherStart and endPos>otherEnd)):
					x+=1
				else:
					listCoordo[x][0]=min(startPos,otherStart,endPos,otherEnd)
					listCoordo[x][1]=max(endPos,otherEnd,startPos,otherStart)
					listCoordo[n][0]=min(startPos,otherStart,endPos,otherEnd)
					listCoordo[n][1]=max(endPos,otherEnd,startPos,otherStart)
					startPos=listCoordo[n][0]
					endPos=listCoordo[n][1]
					x+=1
			n+=1
		listNoDuplicata=list()
		for x in listCoordo:
			if x not in listNoDuplicata:
				listNoDuplicata.append(x)
		listNoDuplicata=sorted(listNoDuplicata)
		finalSequence=str()
		i=0
		while i<len(listNoDuplicata)-1:
			finalSequence+=lineRead[a+1][listNoDuplicata[i][1]:listNoDuplicata[i+1][0]]
			i+=1
			finalSequence=finalSequence.replace(" ","") #Checking for spaces
		if len(finalSequence)>50: #Arbitrairy value... should be set to a much higher value 500/1000 ?!!!!!!!		
			listFinalRead.append(">"+nameRead)
			listFinalRead.append(finalSequence)
		else:
			log.write('Read deleted:\n')
			log.write(nameRead+'\n')
		a+=2
	filereadfinal= open(readpb+"-purifiedTE1","w")
	for u in listFinalRead:
		filereadfinal.write(u+"\n")
	filereadfinal.close()
	return(readpb+"-purifiedTE1")










def CleanerDuplicateBlastn(cheminConsensus,cheminSortieBlastn):
	"""
	Merges if needed the TE identified on the Reads. 
	The TE identified as two xxxxx-5'NNNNN3'-xxxx-5'NNNN3'-xxx instead of xxxx-5'NNNNNNNNNNNNNNNN3'-xxxxx
	Merged them if the total length of the two seperate enteties is less than 20% more than the consensus length.
	ex For a TE of 1000Kb in the consensus: If two framents are detected on a read F1 and F2
	StartF1 to endF2 seperated by less than 1200 they will be fused
	"""
	#Used a dict: TE Id as Key, TE length as value.
	with open(cheminConsensus) as ho:
		listCons=ho.read().split('\n')
	dico_consensus={}
	i=0
	while i<(len(listCons)-1):
		a=listCons[i]
		if a.find(" ")>-1:
			a=a.replace(a[(a.find(" ")):(len(a))],"") #If there is as space i delete it.
		a=a.replace(">","")
		dico_consensus[a]=(len(listCons[i+1]))
		i+=2
	for x in dico_consensus.keys():
		if len(x)<2: 
			print('dico strange key !!',x)#If a key is very wierd...
	with open(cheminSortieBlastn) as thatfile: 
		listHit=thatfile.read().split('\n')
	header=listHit.pop(0)
	lengthListHit=len(listHit)
	i=0 # Organisation of folowing lines: subject id, query id, s. start, s. end, alignment length, evalue, q. start, q. end
	while i<lengthListHit-1:
		ReadTeName=listHit[i][0:listHit[i].find('\t', (listHit[i].find('\t')+1))]
		consensusLength=dico_consensus[listHit[i][listHit[i].find('\t')+1:(listHit[i].find('\t',listHit[i].find('\t')+1))]]
		n=1
		while n<lengthListHit-1:
			if n!=i:
				if listHit[n].find(ReadTeName)>-1:
					lineI=listHit[i].split('\t')
					lineN=listHit[n].split('\t')
					if int(lineI[2])<int(lineI[3]):
						sens="-"
					else:
						sens="+"
					maxim=max(int(lineI[2]),int(lineI[3]),int(lineN[2]),int(lineN[3]))
					minim=min(int(lineI[2]),int(lineI[3]),int(lineN[2]),int(lineN[3]))
					distance=(maxim-minim)
					if distance<(consensusLength+(consensusLength*0.20)):
						if sens=="+":
							## Fields: subject id 0, query id 1, s. start 2, s. end 3, alignment length 4, evalue 5, q. start 6, q. end 7
							newLine=lineI[0]+'\t'+lineI[1]+'\t'+str(minim)+'\t'+str(maxim)+'\t'+lineI[4]+'\t'+lineI[5]#+'\t'+sens+lineI[6]+'\t'+lineI[7]
						else: #sens - le 5' is bigger
							newLine=lineI[0]+'\t'+lineI[1]+'\t'+str(maxim)+'\t'+str(minim)+'\t'+lineI[4]+'\t'+lineI[5]#+'\t'+sens#+lineI[6]+'\t'+lineI[7]
						listHit[i]=newLine
						listHit[n]=newLine
						if len(newLine)==0:
							print(i)
			n+=1 #Small loop
		i+=1 #Big loop
		ii=(i*100)/(lengthListHit)
		sys.stdout.write("\r%d%%" % ii)
		sys.stdout.flush()
	listHit=list(set(listHit))
	outputfilename=open(cheminSortieBlastn+'-CleanPy','w')
	outputfilename.write(header+'\n')
	for x in listHit:
		if len(x)>1:
			outputfilename.write(x+"\n")
	outputfilename.close()
	print('File purged')
	return(cheminSortieBlastn+'-CleanPy')







def LengthExtractorZeroHits(fileSequenceExctraction,blastFile):
	"""
	Exctracts length of fragments that where not alligned. It is required for the final function.
	"""

	with open(blastFile) as blast:
		blast=blast.read().split('\n')

	with open(fileSequenceExctraction) as flank:
		flank=flank.read()#Not splited
	FinalList=list()
	#For each blast line going to searh it in the exctracted file.
	for x in blast:
		if len(x)>1:
			name=x.replace('# Query: ','')#replace the classic output of blast. !!!!!!!! Dependant of blast output
			StartName=flank.find(name)
			if StartName==-1:
				print("Error coudn't find ID in file : "+ name)
				#Written in the logs. If error is detected there will be a trace.
				log.write("Error coudn't find ID in file : "+ name+'\n')
			EndName=StartName+len(name)+1
			EndSequenceLine=flank.find('\n',EndName)
			FinalList.append(name+"\t"+str((EndSequenceLine-EndName))+"\t"+flank[EndName:EndSequenceLine])#Sequence name length of sequence and the sequence
	outputFileName=open(blastFile+'+Length',"w")
	for x in FinalList:
		outputFileName.write(x+'\n')
	outputFileName.close()
	return(blastFile+'+Length')




def LengthExtractorHitPositives(fileSequenceExctraction,blastFile):
	"""
	Exctracts length of positiv blastN matches. (
	Format is different from negativ match hence the new function)
	"""
	
	with open(blastFile) as blast:
		blast=blast.read().split('\n')
	
	with open(fileSequenceExctraction) as flank:
		flank=flank.read()
	FinalList=list()
	x=1#0 is header
	while x<len(blast):
		if len(blast[x])>3:#Avoid empty lines
			X=blast[x].split('\t')
			name=X[4]
			StartName=flank.find(name)
			if StartName==-1:
				print('Error ! Will continue prog, but results for x will be marked in log and should be ignored.',(str(x)))
				log.write('Error name of sequence not found'+str(x)+'\n')
			EndName=StartName+len(name)+1
			EndSequenceLine=flank.find('\n',EndName)
			FinalList.append(name+"\t"+str((EndSequenceLine-EndName))+"\t"+flank[EndName:EndSequenceLine])
		x+=1
	FinalList=list(set(FinalList))
	outputFileName=open(blastFile+'+Length',"w")
	for x in FinalList:
		outputFileName.write(x+'\n')
	outputFileName.close()
	return(blastFile+'+Length')





def CreateSumUp(The53,The5,The3,The053,The05,The03,listOfTransposons,outputName,unresolved,StepProgram):
	"""
	Creats Sum up files. For each TE of the given list it counts the number of occurence in each type of entry: 53,5' 3' only and 0 hits. (takes in count the unresolvable list)
	This sum up is needed for the next step : Mixes all the different files. We finaly identifie the absent of small length >50nt and the polymorph.
	"""
	StepProgram=int(StepProgram)
	threshold=sequencing_dept

	with open(unresolved) as nonresolved:
		nonresolved=nonresolved.read()
	outputName=str(outputName)
	
	with open(listOfTransposons) as listTeName:#List of TE. Step 1 TE list input user. Step 2 List of TE in Reads identified and processed
		listTeName=listTeName.read().split('\n')
	tailleF=len(listTeName)
	
	with open(The53) as The53:
		The53=The53.read()
		The53N=The53.split('\n')

	with open(The3) as The3:
		The3=The3.read()
		The3N=The3.split('\n')

	with open(The5) as The5:
		The5=The5.read()
		The5N=The5.split('\n')
	
	with open(The053) as The053:
		The053=The053.read()
		The053N=The053.split('\n')

	with open(The03) as The03:
		The03=The03.read()
		The03N=The03.split('\n')
		

	with open(The05) as The05:
		The05=The05.read()
		The05N=The05.split('\n')

	listResults=list()
	listResults.append("#Name of TE\tContig\tStart\tStop\tStrand\tPositiveHit-TE\tNegativeHit(Tot)\tNegativeHit-TE<50nt\tNegativeHit-TE>50nt\tResults")
	listPositive=list()
	listMegablast=list()
	listNegativeSmall=list()
	listNegativeLong=list()
	listPoly=list()
	listPoly.append("List of 5-3 sequences found.\nMay not contend the positiv results if they where 5' or 3' only\n")
	listStrange=list()
	i=0
	#Absent
	ab=0
	#Absent len>50nt
	am=0
	#Unsolvable
	un=0 
	#Present
	p=0
	#Non covored megablast alignement 0
	w=0
	#Polymorph
	po=0
	for xname in listTeName:
		bilanImposed='nope'
		listTemporaryPos=list()#Reset each time.
		listTemporaryNeg=list()
		if len(xname)>0:
			nom=xname[0:(xname.find('\t'))] #Name is [0:find('\t')]
			if nom in nonresolved:
				bilanImposed="Too many hits"
			#Simple counts in all the files.
			numberHitPositif=0
			numberHitPositif+=The53.count(nom)
			if The53.count(nom)>0:
				presence53=1
			else:
				presence53=0
			numberHitPositif+=The3.count(nom)
			if The3.count(nom)>0:
				presence3=1
			else:
				presence3=0		
			numberHitPositif+=The5.count(nom)
			if The5.count(nom)>0:
				presence5=1
			else:
				presence5=0
				
			numberHitNegative=0
			numberHitNegative+=The053.count(nom)
			if The053.count(nom)>0:
				presence053=1
			else:
				presence053=0
			numberHitNegative+=The03.count(nom)
			if The03.count(nom)>0:
				presence03=1
			else:
				presence03=0
			numberHitNegative+=The05.count(nom)
			if The05.count(nom)>0:
				presence05=1
			else:
				presence05=0
			#Primary conclusion
			if numberHitPositif<1 and numberHitNegative>0: 
				bilanprov="Absent"

			elif numberHitPositif>=1 : 
				if numberHitPositif>threshold*100 or numberHitNegative>threshold*100:#If we have a verry strange situation
					bilanImposed='Too many hits'
				bilanprov="TE Present"

			else:
				if bilanImposed!='nope':
					bilanprov='Unresolvable locus'
					un+=1
				else:			
					bilanprov='Uncovered locus'
					w+=1
			
			
			if bilanprov=='TE Present': #If present: check if not polymorphism !
				#Check 053 and possible length:
				Big=0
				Small=0
				listForLength=list()
				listSmallFound=list()
				#Checking if it isn't poly:
				if presence053>=1:
					for x in The053N:
						if nom in x:
							listForLength.append(x)
					#Searhing for length
					for gf in listForLength:
						gff=gf.split('\t')	
						if int(gff[1])>50:
							Big+=1
						else:
							Small+=1
							listSmallFound.append(gf)	
					#If small are detected (Present and small)!
					if Small>=1  and bilanImposed!="Too many hits":
						bilanprov="Possible polymorphism"	
						for fx in listSmallFound:
							fxx=fx.split('\t')
							listPoly.append((">"+fxx[0]+" "+fxx[1]+"\n"+fxx[2]))
					
				if presence53>0:
					for x in The53N:
						if nom in x:
							listTemporaryPos.append(x)
				if presence3==1:
					for x in The3N:
						if nom in x:
							listTemporaryPos.append(x)
				if presence5==1:
					for x in The5N:
						if nom in x:
							listTemporaryPos.append(x)
			
				for x in listTemporaryPos:
					x=x.split('\t')
					if bilanImposed=="Too many hits":
						listStrange.append((">"+x[0]+" "+x[1]+"\n"+x[2]))
					else:
						if bilanprov=="Possible polymorphism":
							listPoly.append((">"+x[0]+" "+x[1]+"\n"+x[2]))
							
						else:#The ultimate positiv sequences
							listPositive.append((">"+x[0]+" "+x[1]+"\n"+x[2]))

					# Read namex[O] Sequence lengthx[1] \n and the sequence: Gives us a nice Fasta format
				if bilanImposed=="Too many hits":
					listStrange.append('\n')
					un+=1
					bilanprov="Unresolvable locus"
				else:
					if bilanprov=="Possible polymorphism":
						listPoly.append('\n')
						po+=1
					else:
						listPositive.append('\n')
						p+=1
				Big=(numberHitNegative-Small) #A litle simple but it works.
				
			elif bilanprov=='Absent':
				#If absent I take all the 'The0xx' and then check length of 053.
				prov=list()
				Big=0
				Small=0
				if presence053==1:
					for x in The053N:
						if nom in x:
							listTemporaryNeg.append(x)
				#Length decision
					for x in listTemporaryNeg:
						x=x.split('\t')	
						if (int(x[1]))>50:
							Big+=1
						else:
							Small+=1	
					if Small>0 and Big>0:
						decision='wierd'
					elif Big==0 and Small>0:
						decision='Small'
					elif Big>0 and Small==0:
						decision='Big'
					else:
						decision='Big'
				else:
					decision='Big'
				
				if decision=='Small' or decision=='wierd':
					for x in listTemporaryNeg:
						x=x.split('\t')	
						prov.append((">"+x[0]+" "+x[1]+"\n"+x[2]))#Fasta format
					prov.append('\n')
				else:#Big: taking all the other one.
					if presence053==1:
						for x in The053:
							if nom in x:
								listTemporaryNeg.append(x)
					if presence05==1:
						for x in The05N:
							if nom in x:
								listTemporaryNeg.append(x)
					if presence03==1:
						for x in The03N:
							if nom in x:
								listTemporaryNeg.append(x)
					for x in listTemporaryNeg:
						x=x.split('\t')	
						prov.append((">"+x[0]+" "+x[1]+"\n"+x[2]))#Fasta
					prov.append('\n')				
				Big=(numberHitNegative-Small)
				#Extending the selected list:
				if decision=='Small'and bilanImposed=='nope':	
					listNegativeSmall.extend(prov)
					bilanprov="TE absent"
					ab+=1
				elif decision=='wierd' or bilanImposed!='nope':
					listStrange.extend(prov)
					bilanprov="Unresolvable locus"
					un+=1
				else:
					listNegativeLong.extend(prov)
					bilanprov="Ambiguous negative >50nt"
					am+=1
			#Not found in both:
			else:
				listMegablast.append(nom)
				Small=0
				Big=0
			#Format of tabular output:
			listResults.append(str(xname)+'\t'+str(numberHitPositif)+'\t'+str(numberHitNegative)+'\t'+str(Small)+'\t'+str(Big)+'\t'+bilanprov)
			i+=1
		ii=(i*100)/tailleF
		sys.stdout.write("\r%d%%" % ii)
		sys.stdout.flush()
		
	
	#Pos file:
	if StepProgram==1:
		outputFileNamePositive=outputName+"Step1-Positiv-TE"
		outputFileNameLong=outputName+"Step1-Negativ-TE-Big"
		outputFileNameShort=outputName+"Step1-Absent"
		outputFileNamePoly=outputName+"List-Present-Absent-Polymorph"
		outputFileNameStrange=outputName+"Step1-Strange-Sequence"
		outputFileNameAbsent=outputName+"Step1-NonCovered-Sequence"
		outputFileNameSumUp=outputName+"List-Present-Absent"
	else:
		outputFileNamePositive=outputName+"Step2-Positiv-TE"
		outputFileNameLong=outputName+"Step2-Negativ-TE-Big"
		outputFileNameShort=outputName+"Step2-Insertion"
		outputFileNamePoly=outputName+"ErrorIfNotEmpty"
		outputFileNameStrange=outputName+"Step2-Strange-Sequence"
		outputFileNameAbsent=outputName+"Step2-NonCovered-Sequence"
		outputFileNameSumUp=outputName+"List-new-insertions"
	
	#Positive
	filePosi=open(outputFileNamePositive,'w')
	for x in listPositive:
			filePosi.write(x+'\n')
	filePosi.close()
	
	
	#Negative long:	
	fileNegBig=open(outputFileNameLong,'w')
	for x in listNegativeLong:
		fileNegBig.write(x+'\n')
	fileNegBig.close()
	
	
	#Negative Small (the real interesting file!):	
	fileNegSmall=open(outputFileNameShort,'w')
	for x in listNegativeSmall:
		fileNegSmall.write(x+'\n')
	fileNegSmall.close()
	
	
	#Poly
	filepoly=open(outputFileNamePoly,'w')
	for x in listPoly:
		filepoly.write(x+'\n')
	filepoly.close()
	
	#Strange
	filestrange=open(outputFileNameStrange,'w')
	for x in listStrange:
		filestrange.write(x+"\n")
	filestrange.close()
	
	
	#No coverage megablast
	FileNoCoverage=open(outputFileNameAbsent,'w')
	for x in listMegablast:
		FileNoCoverage.write(x+'\n')
	FileNoCoverage.close()
	
	#General Stats:
	AB=round(((ab*100)/(i+0.00)),2)
	AM=round(((am*100)/(i+0.00)),2)
	P=round(((p*100)/(i+0.00)),2)
	W=round(((w*100)/(i+0.00)),2)
	PO=round(((po*100)/(i+0.00)),2)
	UN=round(((un*100)/(i+0.00)),2)
	
	fileSumUp=open(outputFileNameSumUp,'w')
	LineHeader=("Percentage:\tAbsent\t%s\tAmbiguous\t%s\tPresent\t%s\tNon sequenced\t%s\tPolymorph\t%s\tUnsolvable\t%s\n\n" %((str(AB)),(str(AM)),(str(P)),(str(W)),(str(PO)),(str(UN))))
	for c in listResults:
		fileSumUp.write(c+'\n')
	fileSumUp.close()
	
	return(outputFileNameShort,outputFileNameSumUp,outputFileNamePoly,LineHeader)






def SimpleCleanBlast(fileBlastInput):
	"""
	Keeps the positive alignment only and nothing else. Used for cleaning alignment flanking of step 2 on themself to check is new insertion detected are not coming from the same locus sequences multiple times.
	"""
	listHitPos=list()
	with open(fileBlastInput) as FileBlast:
		FileBlast=FileBlast.read().split('\n')
	for x in FileBlast:
		if len(x)>2:
			if x[0]!="#":
				listHitPos.append(x)
	FileBlastw=open(fileBlastInput+'-HitOnly','w')
	for x in listHitPos:
		FileBlastw.write(x+'\n')
	return(fileBlastInput+'-HitOnly')





def traitementCouvSurEuxMeme(filemegaBflank,fileSequences):
	"""
	Checks if the multiple insertion detected are the result of multiple sequencing of the same locus (wich is a good thing !)
	Regroups all those multiple reads and add the number of occurence
	""" 
	dictio=dict()#Keys are going to be the regrouped results. A key will have all the multiple reads it regroups as his values. (inclueding himself)
	listAlreadyDone=list()
	with open(filemegaBflank) as fi:
		fi=fi.read().split('\n')
	taillefi=str(len(fi))
	for i,x in enumerate(fi): #This time first line isn't alignment header. 

		sep=x.split('\t')
		if len(sep)>1:
			#If a lines aligns itself on itself... and has never aligned before.
			if (sep[0]==sep[1]) and (sep[0] not in listAlreadyDone):
					dictio[sep[0]]=[sep[0]]
					listAlreadyDone.append(sep[0])
				
			#Every other case: the sequence will be added to an existing Key as value.
			else:
				#Both have been already analysed need to check if they belong to the same Key.
				if (sep[0] in listAlreadyDone) and (sep[1] in listAlreadyDone):
					check1=0
					check2=0
					for x in dictio:
						if sep[0] in dictio[x]:
							a=x
							check1+=1
						if sep[1] in dictio[x]:
							b=x
							check2+=1
					if (check1>1 or check2>1):
						print('Error in the program')
						log.write("Error in the program ! ")
					#If a!=b need to pop B into A
					if b!=a:
						dictio[a]+=dictio.pop(b)
				#If one of the two is not in the list
				elif ((sep[0]  in listAlreadyDone)) and ((sep[1] not in listAlreadyDone)):
					listAlreadyDone.append(sep[1])
					for x in dictio:
						if sep[0] in dictio[x]:
							dictio[x]+=[sep[1]]
				elif ((sep[0] not in listAlreadyDone) and (sep[1] in listAlreadyDone)): 
					listAlreadyDone.append(sep[0])
					for x in dictio:
						if sep[1] in dictio[x]:
							dictio[x]+=[sep[0]]
				#Not in the list and two different sequencee:
				else:
					dictio[sep[0]]=[sep[0]]
					dictio[sep[0]]+=[sep[1]]
					listAlreadyDone.append(sep[0])
					listAlreadyDone.append(sep[1])
	
	#Saving the final datas:
	uniqueKeys=list()#Will have all the unique keys
	fileOccurence=open(fileSequences+'Occurence','w')
	for x in dictio:
		fileOccurence.write(x+"\t"+str(len(dictio[x]))+'\n')
		uniqueKeys.append(x)
	fileOccurence.close()
	pourfichier=list()
	for x in uniqueKeys:
		Kro=x[0:x.find('__')]
		sep=x.split('_/_')
		if len(sep)==4:
			start=sep[2]
			stop=sep[3]
			sens=sep[1]
			pourfichier.append((Kro+"\t"+start+"\t"+stop)+"\t"+sens+"\t"+x)
	
	fileRecapAlignment=open(fileSequences+'alignement',"w")
	for x in pourfichier:
		fileRecapAlignment.write(x+'\n')
	fileRecapAlignment.close()
	return(fileSequences+'alignement',fileSequences+'Occurence')








def SeparatorInput(fileSmallID):
	"""
	Separates informations in the lines with _/_ and puts \ t. Used for next step of prog.
	Format of ID lines is very specific: >ID_/_Infos_/_Info2 etc
	"""
	with open (fileSmallID) as fileSmallIDop:
		fileSmallIDop=fileSmallIDop.read().split('\n')
	pourfichier=list()
	for x in fileSmallIDop:
		if len(x)>1:
			if x[0]=='>':
				Kro=x[1:x.find('__')] #pop the >
				x=x[1:x.find(" ")] #pop the > and the end after space.
				sep=x.split('_/_')#Splits the usefull informations.
				if len(sep)==4:
					start=sep[2]
					stop=sep[3]
					sens=sep[1]
					pourfichier.append((Kro+"\t"+start+"\t"+stop)+"\t"+sens+"\t"+x)
	fff=open(fileSmallID+'-',"w")
	for x in pourfichier:
		fff.write(x+'\n')
	fff.close()
	return(fileSmallID+'-')





def rassembleurhitnegatifsequences(file5prim,file3prim,fileNegSmall):
	"""
	Regroups 5' and 3' flanking sequence together so we can align them and regroupe the reads who come from the same locus.
	"""
	with open(file5prim) as prim5:
		prim5=prim5.read()
	with open(file3prim) as prim3:
		prim3=prim3.read()
	with open(fileNegSmall) as negative:
		negative=negative.read().split('\n')
	finalList=list()
	for x in negative:
		if len(x)>0:
			if x[0]=='>':
				sequenceName=x[(x.find('__')):(x.find("_/_"))]#Format __SequenceName_/_OtherInfo
				startSequence=prim5.find('\n',(prim5.find(sequenceName)+len(sequenceName)))+1#+1 to avoid \n
				seq5=prim5[startSequence:(prim5.find('\n',startSequence))]
				#Same for file 3'
				startSequence=prim3.find('\n',(prim3.find(sequenceName)+len(sequenceName)))+1
				seq3=prim3[startSequence:(prim3.find('\n',startSequence))]
				finalList.append(x)
				seqBothFlank=seq5+seq3
				finalList.append(seqBothFlank)
	outputFileName=fileNegSmall+'-FlankingSequences'
	fff=open(outputFileName,'w')
	for x in finalList:
		fff.write(x+'\n')
	fff.close()
	
	return(outputFileName)









def ExtractorFinalSequences(fileNegative,genome,reads,blast):
	"""
	Extracts the final sequence for the results file: Sequences containing the transposon +flanking 5' & 3'.
	The other sequence is the one who points to a TE insertion or deletion (step 1 or 2) ie were the flanking sequences alligned.	
	"""
	finalList=list()
	with open(fileNegative) as fileN:
		fileN=fileN.read().split('\n')
	with open(genome) as ref:
		ref=ref.read()
	with open(reads) as reads:
		reads=reads.read()
	with open(blast) as blast:
		blast=blast.read()
	for x in fileN:
		if len(x)>1:
			lineN=x.split('\t')
			#ID >xxxxx\n -->In genome after the '\n' 
			startSequence=ref.find('\n',(ref.find(lineN[0])+1))+1
			endSequence=ref.find('\n',startSequence)
			sequence=ref[startSequence:endSequence]
			if (int(lineN[1])-length_flank_seq)<0:
				nbdebut=0
			else:
				nbdebut=(int(lineN[1])-length_flank_seq) #!!!!!!!!!!!!!!!!!!!!!!!!!Done need o be sure it still works :) 
			sequenceGenome=sequence[nbdebut:(int(lineN[2])+length_flank_seq)]
			
			if lineN[3]=="-":
				sequenceGenome=Seq(sequenceGenome)
				sequenceGenome=sequenceGenome.reverse_complement()
				sequenceGenome=str(sequenceGenome)
				
			lineSplit=lineN[4].split('__')#Getting the other sequence
			TEname=lineSplit[1] #Loosing _ at beginning need to re add one or'1hoho'='X1hoho'
			beginningB=blast.find('_'+TEname)+1
			endB=blast.find('\n',beginningB)
			lineBlast=blast[beginningB:endB]
			lineBlast=lineBlast.split('\t')
			readName=lineBlast[1]
			mini=min(int(lineBlast[2]),int(lineBlast[3]))
			maxi=max(int(lineBlast[2]),int(lineBlast[3]))
			#Getting DNA sequence
			beginningN=reads.find('\n',(reads.find(readName)+1))+1 #+1 to skip \n
			endN=reads.find('\n',beginningN)
			sequenceReadtTot=reads[beginningN:endN]
			if mini-length_flank_seq<0:
				mini=0
			else:
				mini=mini-length_flank_seq
			sequenceRead=sequenceReadtTot[mini:(maxi+length_flank_seq)]		
			if lineBlast[4]=="-": #Sens for RT if needed.
				sequenceRead=Seq(sequenceRead)
				sequenceRead=sequenceRead.reverse_complement()
				sequenceRead=str(sequenceRead)
				
			finalList.append(("_"+TEname+'\n'+">"+readName+'\t'+(str(mini))+'\t'+(str(maxi+length_flank_seq))+'\n'+sequenceRead+'\n'+'>'+lineN[0]+'\t'+(str(nbdebut))+'\t'+(str(int(lineN[2])+length_flank_seq))+'\n'+sequenceGenome))
	ff=open(fileNegative+"Sequence",'w')
	for x in finalList:
		ff.write(x+'\n\n')
	ff.close()
	return(fileNegative+"Sequence")






def CreateFinalSumup(sortieDiff,lineHeaderSumup,StepNumber,SumUpFile):
	"""
	Inserts in sequences files the number of ocurence and information stored in the sum up file.
	"""
	finalSTR=str()
	#Calculating the number of unique element in the list:
	with open(sortieDiff) as varstc:
		strtot=varstc.read()
	varstc=strtot.split('\n')
	listUniqueID=list()
	for x in varstc:	
		if len(x)>3:
			if x[0]=="_":#Each ID os TE begins with this
				if x not in listUniqueID:
					listUniqueID.append(x)
	
	numberEventsTOT=len(listUniqueID)
	if StepNumber==1:
		lineSumUp=lineHeaderSumup.split('\t')
		a=float()
		a=100.00-float(lineSumUp[8])
		b=float(lineSumUp[10])
		#Adding specific headers: Adding number of occurence in a classy way :')
		for x in listUniqueID:
			nbhitTEunique=strtot.count(x)#Getting the number of occurence
			newStrOccurence=x+'\t'+str(nbhitTEunique)
			strtot=strtot.replace(x,newStrOccurence,nbhitTEunique+20)#Replacing old x with new who has the occurence.(the +20 is to be sure all will be replaced !
		finalSTR+='Percentage of covered TE: '+str(a)+'\n'
		finalSTR+='Percentage of polymorph: '+str(b)+'\n'
		finalSTR+="Number of deletion: "+str(numberEventsTOT)+'\n'
		finalSTR+="\n"
		finalSTR+=strtot
	
	else:#Step two
		with open(filecompteurhit) as OccurenceCount:#File containning the Occurence of each ! Name is global !!!!!
			OccurenceCount=OccurenceCount.read()
		for x in listUniqueID:
			numberOcurrence=OccurenceCount[(OccurenceCount.find('\t',OccurenceCount.find(x))):OccurenceCount.find('\n',(OccurenceCount.find('\t',OccurenceCount.find(x))))]#From the first \t after the ID we want to  the \n after the first find: the number o occurence
			strtot=strtot.replace(x,(x+numberOcurrence))
		finalSTR+="Number of insertion: "+str(numberEventsTOT)+'\n'
		finalSTR+="\n"
		finalSTR+=strtot

	aVa=open(SumUpFile+"-Sequence",'w')
	aVa.write(finalSTR)
	aVa.close()
	return((SumUpFile+"-Sequence"))







def CleanAllDuplicateLines(fileDuplicateLines):
	with open(fileDuplicateLines) as DuplicateLines:
		DuplicateLines=DuplicateLines.read().split('\n')
	listAlreadyDone=list()
	for x in DuplicateLines:
		if x not in listAlreadyDone:
			if len(x)>1:
				if x[0]=="_":
					listAlreadyDone.append("")
			listAlreadyDone.append(x)
	strrr="\n".join(listAlreadyDone)
	outputFile=open(fileDuplicateLines,"w")
	outputFile.write(strrr)
	outputFile.close()
	return(fileDuplicateLines)







#

####

#######

#############

###################

#########################

#################################

########################################

##############################################

#####################################################

##########################################################

################################################################

#######################################################################

#				# main #					#



#Creating sorted TE list: def sorterTE(blastOutput,TE_name,chromosome_name,Start,Stop,senss,collumn_separation,first_blast_line):
separator=str(ordre_output[5])
path_TE_annotated=sorterTE(path_TE_annotated,int(ordre_output[0]),int(ordre_output[1]),int(ordre_output[2]),int(ordre_output[3]),(ordre_output[4]),separator,int(ordre_output[6])) #replaced '\t' by separator !!!!!!!!!!!!!!!!!!!!!!!!!


#Creating Signle line fasta file
print('Input files being converted in single line Fasta style.')
path_pacbio=MonoLineFasta(path_pacbio)
path_ref_genome=MonoLineFasta(path_ref_genome)
log.write('Input files transphormed in single line Fasta style.'+"\n")


#Creating flank 5',3' files
(myFile3prim,myFile5prim)=create_flanq_file(path_TE_annotated, path_ref_genome,name_file_flank)


print("Creation of the References Genome db file.")
time.sleep(1)
makedb(path_ref_genome)


#Creating file of UnresolvableStep1 file
print("Creating list of unsolvable element from reference list")
os.system("echo 'List of non resolvable'>"+name_folder_results+"/UnresolvableStep1")
NonResolvable=name_folder_results+'/UnresolvableStep1'
#Identifiyng TE unresolvable: Flanking sequence aligne more than once on th reference genome.
BlastFlanking5OnRef=MEGABLAST(path_ref_genome,myFile5prim,e_value)
cleaningFlank(BlastFlanking5OnRef,myFile5prim,NonResolvable)
BlastFlanking3OnRef=MEGABLAST(path_ref_genome,myFile3prim,e_value)
cleaningFlank(BlastFlanking3OnRef,myFile3prim,NonResolvable)


#Creating DB for Reads:
print("Creation of the Reads db file.")
time.sleep(1)
makedb(path_pacbio)


#Alignment of flanking sequence on Reads
myFile3primblast=MEGABLAST(path_pacbio,myFile3prim,e_value)
myFile5primblast=MEGABLAST(path_pacbio,myFile5prim,e_value)
print("Cleaning megalast alignment")
CleanerBlastOutput0Hits(myFile3primblast)
CleanerBlastOutput0Hits(myFile5primblast)
print("File 3' & 5' cleaned")
log.write("File 3&5  megablast created and cleaned"+"\n")


#Extracting sequence in the reads.
print("Pooling megablast aligments.")
log.write("Pooling megablast aligments."+"\n")
(flankpooled,sequencePooled,threePrimOnly,fivePrimOnly)=FlankingAlignmentAnalyser(myFile3primblast,myFile5primblast,path_pacbio)


#Limiting he number of sequence
sequencePooledClean=SequenceCleaner(sequencePooled,'both')
threePrimOnlyClean=SequenceCleaner(threePrimOnly,'three')
fivePrimOnlyClean=SequenceCleaner(fivePrimOnly,'five')
print("Files cleaned.")
log.write("Files cleaned"+"\n")


#Db for the TE file who will be used as subject.
print("Creating Db file for TE consensus.")
time.sleep(1)
makedb(path_consensus_TE) 
log.write("DB file done. Alignment of extracted sequences on the consensus"+"\n")


#BLAST(nameFileSequenceExtracted(query),consensus(subject),eeValue(e_value)) (takes all hits not like Blast2)
#Blast des sequences pooled 5-3 then 3' only and 5' only.
nomBlast53=BLAST(sequencePooledClean,path_consensus_TE,e_value)
nomBlast3=BLAST(threePrimOnlyClean,path_consensus_TE,e_value)
nomBlast5=BLAST(fivePrimOnlyClean,path_consensus_TE,e_value)


#Cleaning blast output and putting aligments that failed in a special list.
print("Cleanning blastn output")
zeroHitStepOne53=CleanerBlastOutput0Hits(nomBlast53)
zeroHitStepOne5=CleanerBlastOutput0Hits(nomBlast5)
zeroHitStepOne3=CleanerBlastOutput0Hits(nomBlast3)
print("Output cleaned")
log.write("Blastns cleaned"+"\n")


#Exctracts length of differents elemements(def LengthExtractorZeroHits(fichierSeqExtraite,fichierBlast):)
print("Creating length files for the sumup, can take a while.")
log.write("Creating file containing length of sequences 0 hits Blastn"+"\n")
StepOne053=LengthExtractorZeroHits(sequencePooledClean,zeroHitStepOne53)
StepOne05=LengthExtractorZeroHits(fivePrimOnlyClean,zeroHitStepOne5)
StepOne03=LengthExtractorZeroHits(threePrimOnlyClean,zeroHitStepOne3)
#LengthExtractorHitPositives(fichierSeqExtraite,fichierBlast):
log.write("Creating file containing length of sequences Blastn"+"\n")
StepOne53=LengthExtractorHitPositives(sequencePooledClean,nomBlast53)
StepOne5=LengthExtractorHitPositives(fivePrimOnlyClean,nomBlast5)
StepOne3=LengthExtractorHitPositives(threePrimOnlyClean,nomBlast3)


log.write("Generating sum up files \n")
#def bila nteur(le53,le5,le3,le053,le05,le03,listetrouver)
#Cycle1 Sum up:
#nameOutputSumupOne=name_file_flank+'Step1'
nameOutputSumupOne=name_folder_results+"/"
(smallNegative,sumUpStepOne,polymorph1,Line_Header)=CreateSumUp(StepOne53,StepOne5,StepOne3,StepOne053,StepOne05,StepOne03,path_TE_annotated,nameOutputSumupOne,NonResolvable,1)
log.write("Divers data for the sum up:\n %s \n\n" %(Line_Header))

smallNegativeSeparated=SeparatorInput(smallNegative)
smallNegativeSequences=ExtractorFinalSequences(smallNegativeSeparated,path_pacbio,path_ref_genome,path_TE_annotated)
smallNegativeOccurence=CreateFinalSumup(smallNegativeSequences,Line_Header,1,sumUpStepOne)


#Polymorphs:
polymorph1Separated=SeparatorInput(polymorph1)
polymorph1Final=ExtractorFinalSequences(polymorph1Separated,path_pacbio,path_ref_genome,path_TE_annotated)



###
####Starting second cycle:
###


#Cleanning read from previous detected transposons.
print('Cleaning PacBio reads')
log.write('Cleaning PacBio reads'+"\n")
path_reads_clean=PacBionetoyeurBlastn(nomBlast5,nomBlast3,nomBlast53,path_pacbio)#!!!!!!!!!!!! The one function I still can't convert...
log.write("Cleaning done"+"\n")


#Blast TE consensus: Identify all the remaining TE and possible new insertions.
path_consensus_TE=MonoLineFasta(path_consensus_TE)
makedb(path_reads_clean)
log.write("Db file on PacBio reads done"+"\n")
print("Blastn: Alignments of the TE consensus on the reads.")
fileBlastSecondStep=BLAST2(path_consensus_TE,path_reads_clean,e_value)#Max target 1==> only the consensus that match the most is saved.
log.write("Second blastn done."+"\n")
log.write("Cleaning Blast output "+"\n")
TEabsentReadCycle2=CleanerBlastOutput0Hits(fileBlastSecondStep)
fileBlastSecondStep=CleanerDuplicateBlastn(path_consensus_TE,fileBlastSecondStep)

#1- sorterTE(SortieBlast,nomTE,nomKro,Start,Stop,sens,separation_collones):
#Remplie celon les paramettre de blast
log.write("Sorting TE of second step"+"\n")
print('Creating proper TE list for the program.')
ListTeSecondStep=sorterTE(fileBlastSecondStep,1,0,2,3,'+++++++','\t',1)


#Creating flank file of step two.
log.write("Creating step 2 flanking sequence files."+"\n")
(Flank3primV2,Flank5primV2)=create_flanq_file(ListTeSecondStep,path_reads_clean,name_file_flank_step2)


#Alignment of flanking sequence on reference genome.
print('Alignment of the flanking sequence on the reference genome')
myFile3primblastV2=MEGABLAST(path_ref_genome,Flank3primV2,e_value)
myFile5primblastV2=MEGABLAST(path_ref_genome,Flank5primV2,e_value)
log.write("Megablast of flanking sequences on the reference genome done."+"\n")


#cleaning Flank files:
print('Cleaning megablast output.')
CleanerBlastOutput0Hits(myFile5primblastV2)
CleanerBlastOutput0Hits(myFile3primblastV2)
NonResolvable2=name_folder_results+'/UnresolvableStep2'
myFile5primblastV2=cleaningFlankMegablast(myFile5primblastV2,Flank5primV2,NonResolvable2)
myFile3primblastV2=cleaningFlankMegablast(myFile3primblastV2,Flank3primV2,NonResolvable2)


#Extracting sequence of interest:
print("Comparing megablast alignments")
log.write("Comparing megablast alignments"+"\n")
(flankpooledV2,sequencepooledV2,threePrimOnly2,fivePrimOnly2)=FlankingAlignmentAnalyser(myFile3primblastV2,myFile5primblastV2,path_ref_genome)
print("Sequence of interest exctracted")


#Cleaning sequences:
print("Cleaning extracted sequences")
sequencepooledV2pur=SequenceCleaner(sequencepooledV2,'both')
threePrimOnly2Clean=SequenceCleaner(threePrimOnly2,'three')
fivePrimOnly2Clean=SequenceCleaner(fivePrimOnly2,'five')



#Blast of the extracted sequence on the TE consensus:
print("Alignment of the sequence of interest on the consensus")
log.write("Alignment of the extracted sequence on the TE consensus, and cleaning 0 hits."+"\n")
blast253=BLAST(sequencepooledV2pur,path_consensus_TE,e_value)
zerohitC253=CleanerBlastOutput0Hits(blast253)
blast23=BLAST(threePrimOnly2Clean,path_consensus_TE,e_value)
zerohitC23=CleanerBlastOutput0Hits(blast23)
blast25=BLAST(fivePrimOnly2Clean,path_consensus_TE,e_value)
zerohitC25=CleanerBlastOutput0Hits(blast25)


#Creating the other files needed for the sum up:
print("Creating files needed for the sum up of step two.")
log.write("Creating needed file for the sum up: \n")
C2053=LengthExtractorZeroHits(sequencepooledV2pur,zerohitC253)
C203=LengthExtractorZeroHits(threePrimOnly2Clean,zerohitC23)
C205=LengthExtractorZeroHits(fivePrimOnly2Clean,zerohitC25)
C253=LengthExtractorHitPositives(sequencepooledV2pur,blast253)
C25=LengthExtractorHitPositives(threePrimOnly2Clean,blast23)
C23=LengthExtractorHitPositives(fivePrimOnly2Clean,blast25)


#Creating sum up of step two:
log.write("Creation of the sum up: \n")
nameOutputSumUpTwo=name_folder_results+"/"
(newInsertions,FileSumUpTwo,polymorph2,Line_Header2)=CreateSumUp(C253,C25,C23,C2053,C205,C203,ListTeSecondStep,nameOutputSumUpTwo,NonResolvable2,2)


#Checking and joinning all new insertion found: (If they come from the same locus sequenced multiple times)
print("Joining multiple sequenced locus as one ")
aa=rassembleurhitnegatifsequences(Flank5primV2,Flank3primV2,newInsertions)
makedb(aa)
print('Alignment of the flanking sequences on themself')
at=SimpleBlast(aa,aa)
at=SimpleCleanBlast(at)
print('Merging elements.')
(aaa,filecompteurhit)=traitementCouvSurEuxMeme(at,aa)


#ExtractorFinalSequences(filefinaleneg,genome,reads,blast):
print('Creating final results.')
zoooo=ExtractorFinalSequences(aaa,path_ref_genome,path_reads_clean,ListTeSecondStep)
finalite2=CreateFinalSumup(zoooo,Line_Header,2,FileSumUpTwo)


#Final clean of final fills: clears all duplicate lines for an easy to read file
CleanAllDuplicateLines(smallNegativeOccurence)
CleanAllDuplicateLines(polymorph1Final)
CleanAllDuplicateLines(finalite2)







#Creating the polymorph step 2 files, who depends on the finalite2 file and the reads purified:
results=finalite2
path_reads=path_reads_clean

zeList=list()
with open (aa) as flankJoined:#aa is the flanking sequence joined
 flankJoined=flankJoined.read()
with open(results) as insertion:
 insertion=insertion.read()
 insertions=insertion.split('\n')


for x in insertions:#Taking all the flanking sequences that made it in the final file:
 if len(x)>0:
  if x[0]=="_":
   zeList.append((">"+x))#Name
   nameNoIndicator=x[0:(x.find("\t"))]
   if nameNoIndicator==-1:
    print(x,"/has/failed/finding/",nameNoIndicator)
    log.write("Error in polymorph2 step, finding name in file:"+'\n'+">"+x+'\n')
   else:
    zeList.append(flankJoined[flankJoined.find('\n',flankJoined.find(nameNoIndicator))+1:flankJoined.find('\n',1+flankJoined.find('\n',flankJoined.find(nameNoIndicator))+1)])

saveFlank=open(aa+'cured',"w")
for x in zeList:
 saveFlank.write(x+'\n')

saveFlank.close()

#File name
fileFlanquing=aa+'cured'








blastForPoly=MEGABLAST(path_reads,fileFlanquing,'1e-40')
#print('File output')
#print(blastForPoly)
#print('Opening file output')








with open(blastForPoly) as blastFilePoly:
 blastFilePoly=blastFilePoly.read().split('\n')


#Taking all the positiv match with no TE in it =3/4 of length of both flanking sequences:
listOfGood=list()
for x in blastFilePoly:
 if len(x)>0:
  if x[0]!="#":
   y=x.split('\t')
   #print("Hahahahahha")
   #print(y[3])
   if int(y[3])>length_flank_seq*2*float(3/float(4)) and int(y[3])<length_flank_seq*4:#Alignment should be the length of both 5' and 3'. Not too short and not too long.
    listOfGood.append(y)
    #print(y[3])


with open (path_reads) as readsTot:
 readsTot=readsTot.read()


dicoFullPolymorph=dict()
for x in listOfGood:
 reverse=False 
 lineRead=str()
 TEname=x[0]
 readName=x[1]
 start=min(int(x[8]),int(x[9]))
 end=max(int(x[8]),int(x[9]))
 startSequence=readsTot.find('\n',(readsTot.find(readName)+1))+1
 endSequence=readsTot.find('\n',startSequence)
 lineRead=readsTot[startSequence:endSequence]
 sequence=lineRead[start:end]
 if int(x[8])>int(x[9]):
  sequence=Seq(sequence)
  sequence=sequence.reverse_complement()
  sequence=str(sequence)
  reverse=True
 if len(sequence)<10:
  print('Fail')
  break
 else:
  if reverse:
   header=">"+readName+"\t"+str(end)+"\t"+str(start)+"\tR"
  else:
    header=">"+readName+"\t"+str(start)+"\t"+str(end) 
  if TEname in dicoFullPolymorph:
   dicoFullPolymorph[TEname].append(header)
   dicoFullPolymorph[TEname].append(sequence)
  else:
   dicoFullPolymorph[TEname]=[header,sequence]


#Creating final list:
insertionMoveToPoly=list()
for x in dicoFullPolymorph:
 tooMove=(insertion[insertion.find(x):insertion.find("\n\n_",insertion.find(x))])
 insertionMoveToPoly.append(tooMove)
 for y in dicoFullPolymorph[x]:
  insertionMoveToPoly.append(y)
 insertion=insertion.replace(tooMove,"")

#Creating final files
lastStep=open(results+'-Polymorph',"w")
lastStep.write((str(len(dicoFullPolymorph)))+" insertions have been identified as having at least one polymorph variants in the reads.")
for x in insertionMoveToPoly:
 if len(x)>0:
  if x[0]=="_":
   lastStep.write("\n\n")
   lastStep.write(x+"\n")
  else:
   lastStep.write(x+"\n")

lastStep.close()

numberInsert=int((insertion[insertion.find(': ')+2:insertion.find('\n')]))
finalNumber=numberInsert-len(dicoFullPolymorph)
insertion=insertion.replace(str(numberInsert),str(finalNumber),1)
 
finalStep=open(results,'w')
if finalNumber==0 and numberInsert!=0:
 finalStep.write('All insertions found have been move to Polymorph file')
finalStep.write(insertion)

finalStep.close()









#Creation of output summary folder:
my_summary_file=name_folder_results+'/Summary'
os.system('mkdir -p '+my_summary_file) 


#Moving files of interest in summary folder:
os.system("mv "+smallNegativeOccurence+" "+my_summary_file)
os.system("mv "+polymorph1Final+" "+my_summary_file)
os.system("mv "+sumUpStepOne+" "+my_summary_file)


os.system("mv "+finalite2+" "+my_summary_file)
poly2file=results+'-Polymorph'
new_name=poly2file.replace("-Sequence-Polymorph","-Polymorph-Sequence")
os.system("mv "+poly2file+" "+new_name)
os.system("mv "+new_name+" "+my_summary_file)

#Last litle print for the user and log:
print("LoRTE finished ! Time to analyse the data !")
log.write("Step 2 finished, LoRTE."+"\n"+"Go tiger !")



#""" If you want to keep all the files
#Deleting all unecessary files:
#All the files 
def SupressionOfFiles(path,list_of_name):
	for x in list_of_name:
		file_deleted=str(path+x)
		os.remove(file_deleted)



#Deleting all the files related to the reference Genome: from here
SupressionOfFiles(path_ref_genome,[".nhr",".nin",".nsq"])
SupressionOfFiles(path_TE_annotated,[""])
SupressionOfFiles(path_pacbio,[".nhr",".nin",".nsq"])
SupressionOfFiles(path_reads_clean,["",".nhr",".nin",".nsq"])
SupressionOfFiles(path_consensus_TE,[".nhr",".nin",".nsq","-blastn","-blastn-ZeroHit","-blastn-CleanPy","-blastn-CleanPy-Sorted-LorTE"])

#Deleting all the files in the output results:
os.chdir(name_folder_results)
if os.getcwd()==(name_folder_results):
	os.system("rm Out-*")
	os.system("rm Out2*")
	os.system("rm Step1-*")
	os.system("rm Step2-*")
	os.system("rm UnresolvableStep*")
	os.system("rm List-Present-Absent-Polymorph*")
	os.system("rm ErrorIfNotEmpty")
        os.system("rm "+FileSumUpTwo)

#to here""" 


#Clossing all files and stopping subprocess:
log.close()
#subprocess.Popen.kill(ram_monitoring)

