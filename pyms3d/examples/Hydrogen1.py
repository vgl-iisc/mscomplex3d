import pyms3d


#data info
DataFile = "Hydrogen_128x128x128.raw"
Dim      = (128,128,128)

print pyms3d.get_hw_info()

# Create the mscomplex object 
msc = pyms3d.mscomplex()

# Compute the Mscomplex
msc.compute_bin(DataFile,Dim)

# Simplify the complex
msc.simplify_pers(thresh=0.05)

# print out the 2 saddles
print "2Saddles Ids = " + str(msc.cps(2)) + "\n"

# print out the maxima
print "Maxima Ids = " + str(msc.cps(3)) + "\n"

#print out connectivity
print "-"*36 + "\nSad. \tMax. \tMultiplicity\n"+"-"*36
for sad in msc.cps(2):
    for c in msc.asc(sad):
        print "%d\t%d\t%d"%(sad,c[0],c[1])
