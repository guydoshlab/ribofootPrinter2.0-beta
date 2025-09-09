import posavg
import sys

# This caller script opens the txt input file, parses it, and for each rocc file specified, calls posavg.
# Pull in the name of the file with input params. Filename is only input for the caller.
args1 = (sys.argv)
if len(args1)!=2:
	print("Error.")
else:
	f=open(args1[1])
	
# Parse the file and extract the inputs.
inputparams=[line.strip() for line in f.readlines()] 
f.close()

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
posavg.main(variables["roccnames"],variables["motif"],variables["kind"],variables["frame"],variables["bkndwindowthresh"],variables["bkndwindow"],variables["ORFnorm"],variables["shift"],variables["UTRmode"],variables["subsetlist"],variables["outfile"])