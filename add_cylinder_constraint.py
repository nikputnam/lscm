#!/usr/bin/python3

import sys, getopt

def main(argv):
   rows = ''
   cols = ''
   infile = ''
   try:
      opts, args = getopt.getopt(argv,"hr:c:")
   except getopt.GetoptError:
      print ('test.py -r <rows> -c <cols> infile')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('test.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-r"):
         rows = int(arg)
      elif opt in ("-c"):
         cols = int(arg)
   infile = args[0]
   #print("{} x {}; {}".format( rows, cols, infile))
   
   with open(infile,'r') as f:
    counts = { "v": 0, "vt": 0, "vn":0, "f":0 }
    line = f.readline().strip().split()
    n = 0
    #print(n,line)

    while len(line)>0:
        n = counts.get(line[0],0)
        #print(n,line)
        i = n%rows
        j = int(n/rows)
        #if n:
            #if (line[0]=="vt"):
        #print(n,i,j,line)
        if (line[0] == "v"):
            if ((j==0) or j==(cols-1)):
                if j==0 and i==0:
                    print(" ".join(line+["fix","0","0"]))
                elif j==cols-1 and i==0:
                    print(" ".join(line+["fix","1","0"]))
                else:
                    print(" ".join(line+["zip",str(i)]))
            else:
                print(" ".join(line))

        else:
            print(" ".join(line))
        counts[line[0]] = 1+counts.get(line[0],0)
        line = f.readline().strip().split()


if __name__ == "__main__":
   main(sys.argv[1:])
