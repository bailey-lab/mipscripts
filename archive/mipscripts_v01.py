#!/usr/bin/python3

VERSION=0.01_20211004
###standard modules
import argparse #parse initial arguements from command line
import re    
import os
import os.path
import pickle
import sys
import subprocess
import textwrap
import glob
import copy
import csv
from collections import defaultdict, Counter, OrderedDict
#from collections import namedtuple

#key science modules
import numpy
import scipy
#import numpy.string_
#import csv 
## key bio python modules 

import scipy.interpolate



###################################################################
#2021.10.04 added more warnings 

#TODO allow option to merge over replicates#

################################################################### 
# INITIALIZE GLOBAL VARIABLES ##################################### 
################################################################### 

subcommands={} 

global_samples_tsv_header=[ 'sample_name', 'sample_set', 'replicate',   'library_prep',	'probe_set', 'capture_plate', 'quadrant',	'library_prep',	 'fw', 'rev', 'capture_plate_row', 'capture_plate_column',	'sample_plate', 'sample_plate_column', 'sample_plate_row','FW_plate','REV_plate', ]

global_samples_tsv_oldnames={
"Library Prep": "library_prep",
"384 Column": '384_column',
"Sample Set": "sample_set"
}

#below is all the different headers in *sample_tsv 
# find . -maxdepth 2 -name '*samples.tsv' -exec head -1  {} \; | tr ' ' 'x'  | xargs -n 1 | sort | uniq

#384xColumn
#capture_plate
#capture_plate_column
#capture_plate_row
#CapturexPlatexLocation
#CapturexPlatexName
#column
#diff
#fw
#FW_plate
#library_prep
#LibraryxPrep
#owner
#probe_set
#quadrant
#replicate
#rev
#REV_plate
#row
#sample_name
#sample_plate
#sample_set
#SamplexSet



###################################################################
###################################################################  

###################################################################
#  MAIN        ####################################################
###################################################################



def main( args ):
  """Main allows selection of the main subcommand (aka function).
  Each subcommand launches a separate function. The pydoc subcommand 
  launches pydoc on this overall program file. 
  :param args: the main command line arguments passed minus subcommand
  """
  #print globals().keys()
  #print "ARGUMENTS", args
  
  if len(args)  == 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
    verbosity= 'shortDesc'
    if args[0] in ["help" , "--help", "-help"]:
      verbosity = 'longDesc' 
    program_name=os.path.basename(__file__)
    print ("VERSION: ", VERSION)
    print ("USAGE:",program_name, "[-h] subcommand [suboptions]")
    print ("DESCRIPTION: various scripts complementing MIPTools pipelines")
    print ("SUBCOMMANDS:")
    #tw=TextWrap()
    for k in subcommands.keys():
      text=subcommands[k][verbosity]
      text= textwrap.dedent(text)
      if text:
        text =  "%s:   %s " %(k, text )
        print (textwrap.fill(text,77, initial_indent='', subsequent_indent='         ') )
    print ("HELP:")
    print ("pydoc      detailed documentation of program structure")
    print ("-h/-help   short / long  subcommand descriptions")
    print ("For specific options:",program_name,"[subcommand] --help")
  elif args[0] == '--pydoc':
    os.system( "pydoc " + os.path.abspath(__file__) )
  elif args[0] in subcommands.keys():
    globals()[args[0]](args[1:])
  else:
    print ("unknown subcommand (" + args[0] + ") use -h for list of subcommands!")
    sys.exit(-1)
  sys.exit(0) #normal exit

#------------------------------------------------------------------------------
###############################################################################
####  MERGE FASTQs and SAMPLESHEET for MIPSET  ################################
###############################################################################
#------------------------------------------------------------------------------

shortDescText="merge a sampleset from multiple samplesheets and sequencing runs"
longDescText=""" pull sampleset from multiple samplesheets and merge the R1 and R2 fastq records """
subcommands['merge_sampleset'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

#./mipscripts_v01.py  merge_sampleset --name jjj --samplesheet /work/bailey_share/raw_data/190312 --samplesheet /work/bailey_share/raw_data/190321_nextseq


def merge_sampleset(args):
  """ commandline routine to merge seqeuncing runs and mip captures for one or a few samplesets and probesets.  NOTE:  best used for a singular probe_set and sample_set.  
  """
  aparser=argparse.ArgumentParser( prog = os.path.basename(__file__)+" merge_sampleset", 
      description="allows for flexible merging: preferable to analzye one sampleset by one probeset",
      epilog="Note: for output recommend using unique name for merged sheet based on set/probe/and other parameters",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
      
  aparser.add_argument("--set", required=True, action='append', help='sampleset name to aggregate')
  aparser.add_argument("--probe", required=True, action='append', help='probe sets to include in merge')
  aparser.add_argument("--sheet", required=False, action="append", help=' sample sheet file paths (dir fastq must be in same dir as sheet)')
  aparser.add_argument("--mergeon", required=False, help='fields to merge on ', default ="sample_name-sample_set-replicate")
  aparser.add_argument("--newsheet", required=False, help='name of new merged samplesheet', default="mergedsheet.tsv")
  aparser.add_argument("--newfastqdir", required=False, 
  				help='name of new fastq directory',default="mergedfastq")
  aparser.add_argument("--exclude" , required=False, action ="append", help = 'exclude mergeon patterns matching the text')
  aparser.add_argument("--skipnewfastq", action='store_true', help='dry run skips fastq merge but makes tsv')
  aparser.add_argument("--collapse", action='store_true', help='collapse to unique values in columns')
  aparser.add_argument("--ignorereplicateredundancy", action='store_true', help='collapse multiple replicates')
  
  ##TODO: set default mergedsheet and fastq directory to set_probe
  
  args=aparser.parse_args(args=args)
  mergeonfields=args.mergeon.split('-')
  print (args)
 
  total_samples = []
  merged={}
  replicate={}
  header=global_samples_tsv_header[:]  #copy this as it will be modified
  headers_to_rename=global_samples_tsv_oldnames; #this will not modify
  fieldstocollapse=['probe_set', 'replicate_old',  'library_prep', 'capture_plate', 'quadrant',	  	 'fw', 'rev', 'capture_plate_row', 'capture_plate_column',	'sample_plate', 'sample_plate_column', 'sample_plate_row','FW_plate','REV_plate','384_column', 'fastq']
  
  print ("CREATING and CLEANING NEW FASTQ DIR:" , args.newfastqdir )
  if len (args.newfastqdir) > 0:
    os.system("mkdir " + args.newfastqdir )
    os.system("rm  " + args.newfastqdir + "/*" )
  else:
    print ("ERROR: bad directory name ("+args.newfastqdir + ")")
    exit()

  print ("PROCESSING SAMPLE SHEETS...")
  for thefile in args.sheet:
    print ("LOADING", thefile ,  "...")
    fastqdir=os.path.dirname(thefile) + "/fastq"
    fastqs =os.listdir(fastqdir)
    fastqs = [ f for f in fastqs if any(s in f for s in args.set) ]
    #print (fastqs[0:5])
    print ( "...   ",  len(fastqs)    , " fastqs in associated directory" )
    with open( thefile, newline = '') as items_file:                                                       
      items_reader = csv.DictReader( items_file, delimiter='\t' )
      items_total=0
      items_kept=0
      items_excluded=0
      fastqs_kept=0;
      for item in items_reader:
        #rename old names to new names#
        for oldname in headers_to_rename:
      	  if oldname in item:
            item[ headers_to_rename[oldname] ]=item[oldname]  
            item.pop(oldname)  	 
        if items_total==0:
          #process the headerline to see what headers are present
          #determine if we need to add another line to current header. 
          #print (item)
          for name in item:
            if not name in header:
              header.append(name) 
        items_total+=1
        if not any (s == item['sample_set']  for s in args.set ):
          continue
        if not any ( p in item['probe_set'] for p in args.probe ) :
          continue
        #keep item with proper probeset and sampleset
        items_kept+=1
        startfastqname=item['sample_name']+'-'+item['sample_set']+'-'+item['replicate']+'_'
        #item['fastqdir']=fastqdir
        item['replicate_old']=item['replicate']
        item['fastq']=[ fastqdir+"/"+i for i in fastqs if i.startswith(startfastqname) ]
        item['fastq'].sort()
        fastqs_kept+=len(item['fastq'])
        item['mergeon']=[]
        for f in mergeonfields:
          item['mergeon'].append(item[f])
        item['mergeon']="-".join(item['mergeon'] )
        if   args.exclude and any(ele in item['mergeon'] for ele in args.exclude) :
          items_excluded +=1
          continue
        #print (startfastqname)
        #find the fastqs
        total_samples.append(item)
        merged[item['mergeon']]=0
        replicate[item["sample_name"]+"-"+ item["sample_set"]]=1
      print ( "    ", items_excluded , " excluded samples")
      print ( "   ", items_kept, " of ", items_total, " total with ", fastqs_kept, "fastqs kept")
   # print ("HEADER",header)
    #print (EXtotal_samples[0:1])
  print ("Total samples ", len(total_samples), " merging into ", len(merged), " samples")
  total_samples.sort(key= lambda x: x['mergeon'])
    
  #merge the list and copy the files ####################
  merged_samples=[]

  count=0;
  for mkey in sorted(merged):
     mergeset=[i for i in total_samples if i['mergeon']==mkey]
     count+=1
     print (count, mkey, "records:",len(mergeset))
     mergeditem=None
     name= None  
     ### create a merge record from first one ###
     for i in mergeset:
       if mergeditem==None:
         mergeditem=copy.deepcopy(i)
         name= mergeditem["sample_name"]+ "-" + mergeditem["sample_set"]
         mergeditem["replicate"]= replicate[name]
         replicate[name] +=1
         for f in fieldstocollapse:
           if f in mergeditem:
             mergeditem[f]=[mergeditem[f]]
       else :
         for f in fieldstocollapse:
           if f in mergeditem:
             mergeditem[f].append (i[f]) 
     ### merge the fastqs ##########
     writeoperator=' > '  ##initial write
     for pair in mergeditem['fastq']:
       if len(pair) ==0:
         continue
       elif (len(pair) % 2) == 0:  #2,4,6,8
         #expectation is each replicate record has pair of fastqs (R1 & R2)
         #potential to have more than one replicate in a sequencing run but that is wierd
         #print (pair)
         if len (pair) > 2 and args.ignorereplicateredundancy==False:
           print ("ERROR? More than just a single pair for replicate!",  pair
                   , "\nThis could be due to not properly numbering replicates!"
                   , "\nOr multiple demultiplexing runs into same directory?"
                   , "\nOr because you meant to..."
                   , "\nTo override use --ignorereplicateredundancy!")
           badcount=len([m for m in total_samples if len(m['fastq'])>2])
           print ("NOTE: ", badcount, "have more than  2 fastqs!")
           exit(1)
         if not args.skipnewfastq: #this skips lengthy process of writin

           os.system( "cat  "+ pair[0]  +writeoperator + args.newfastqdir + "/"+ name + "-" 
              + str(mergeditem['replicate']) + "_R1_001.fastq.gz" )
           os.system ( "cat  "+ pair[1]  +writeoperator + args.newfastqdir + "/"+  name + "-" 
              + str(mergeditem['replicate']) + "_R2_001.fastq.gz" )
           writeoperator=' >> ' #now appending after initial files
         #else:
         #  for i in range (0, len(pair),2):
         #     print (pair[i],"\n   ", pair[i+1])
       else:
         print ("ODD PAIR ERROR ",  pair)
         print (" could a R1 or R2 fastq file been deleted???")
         exit (1)
       writeoperator=' >> ' #now appending
     ### make fields nonredudnant (always make probeset)
     uniquefields= ['probe_set']
     mergeditem['fastq']=[fq for sublist in mergeditem['fastq'] for fq in sublist ]
     #print (mergeditem)
     if args.collapse :
       uniquefields=fieldstocollapse
     for f in uniquefields :
       if f in mergeditem:
          #print (f, mergeditem[f])
          mergeditem[f]= list(set(mergeditem[f]))
     for f in fieldstocollapse :
       if f in mergeditem:
         mergeditem[f]= ",".join(mergeditem[f])
     merged_samples.append(mergeditem)
     #print ("NEW COPY ", mergeditem)
  #write to a file ################################
  #write header and additional data from merge
  addedoutput=['fastq','replicate_old', 'mergeon']
  print ("WRITING TSV: ", args.newsheet)
  with open ( args.newsheet, "wt") as out_file: 
     tsv_writer=csv.writer(out_file,delimiter="\t")
     tsv_writer.writerow(header + addedoutput)
     for item in merged_samples:
       row=[]
       for h in header:
         if  h in item:
           row.append(item[h])
         else:
           row.append('')
       for x in addedoutput:
         if x in item:
           row.append(item[x])
         else:
           row.append('')
       tsv_writer.writerow(row)
         

#------------------------------------------------------------------------------
#main is left to the bottom so that global variables can be populated
if __name__ == "__main__":
  print ("SYS VERSION",)
  print (sys.version)
  print ("WORKING DIRECTORY", os.getcwd() )
  print ("COMMAND:" ) 
  print ( " " . join(sys.argv) )
  if len (sys.argv)==1:
    sys.argv.append("--help")  #if no command then it is a cry for help
  main(sys.argv[1:])





