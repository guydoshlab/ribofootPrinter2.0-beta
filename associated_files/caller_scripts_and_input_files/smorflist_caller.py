import smorflist
import sys

# This caller script opens the txt input file, parses it, and for each rocc file specified, calls smorflist.
# Pull in the name of the file with input params. Filename is only input for the caller.
# Wrapper function for finding dORFs.	
args1 = (sys.argv)
if len(args1)!=2:
	print("Error.")
else:
	f=open(args1[1])
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	tossedgenes=0

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	# Now call the wrapper function
	smorflist.main(variables["roccnames"],variables["lengththresh"],variables["shift"],variables["smallest"],variables["mismatches"],variables["UTR"],variables["outfile"])
