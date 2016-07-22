#!/usr/bin/env python

'''
#
# Create a GAF converter that merges  GPAD/GPI files into a GAF file -
# This way loads/reports scripts continue to use the GAF format as MGI standard format for GO
# and are independent to further file type changes from the external source.
# Our GO loads will be able to download GPAD/GPI files then use the converter
# to generate the GAF file
#
#GPAD2GAF converter
#         Uses the GPAD/GPI files to generate a GAF file
#
#
# MGI GO loads download the following GAF files:
#  a. gene_association.mgi GAF
#  b. gene_association.mgi_nonmouse file
#  c. gene_mp_association.mgi file
#
#
#Assumption:
# 1. The order of the first 17 fields of a GAF file does not change
# 2. The order of the first 11 fields of a GPAD file does not change
# 3. The order of the first 9 fields of  a GPI file does not change
#
#Author: Lucie Hutchins, Scientific Software Engineer - MGI
#
# History:
#
# lnh	TR11630
#	- created - 01/2015
#
'''
 
import getopt, sys 
import config,goa_parser 
import os
from datetime import datetime

#

#
#gpad2gaf usage
#
def gpad2gaf_usage():
    print """\
    \n********************************\ngpad2gaf converts a GPAD and the associated GPI file into a GAF file
    \nAssumptions: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa for specs
        1. The order of the first 17/15 fields of a GAF2.0/GAF1.0 file follows goa specs
        2. The order of the first 11 fields of a GPAD file follows goa specs
        3. The order of the first 9 fields of  a GPI file follows goa specs
    \nUsage: 
    Example1: gpad2gaf.py --help  => to display this help page
    Example2: gpad2gaf.py   --gpad=gpad_file --gpi=gpi_file [--gaf=gaf_file] [--version=gaf_version]
    Example: gpad2gaf.py   --gpad=path2/gene_association.mgi.gpad 
             --gpi=path2/gene_association.mgi.gpi --gaf=path2/gene_association.mgi.gaf --version=2.0
    Where:
       --gaf      => <optional> specifies the name of the resulting gaf file (default gpad_file.gaf)
       --version  => <optional> specifies the version of the resulting gaf file (default 2.0)
       --gpad     => <required> specifies the path/name of the input gpad_file
       --gpi      => <required> specifies the path/name of the input gpi_file
    \nNote: If you do not provide the name of the gaf file to generate, the program will create 
       a gaf file in the same directory the input gpad file resides with the extension *.gaf
    \n********************************
    """ 
#
# Main
#
def main():
    log_dir=config.GAF_CONVERTER_LOG_BASE
    if not os.path.isdir('%s' % (log_dir)):
          os.system('mkdir %s' % (log_dir))
    log_file=log_dir+"/"+os.path.basename(sys.argv[0])+".log"
    log=open(log_file,"w")
    i = datetime.now()
    log.write("Program Starts: "+i.strftime('%Y/%m/%d %I:%M:%S %P'))
    log.write("\n")
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:p:i:v:", ["help", "gaf=","gpad=","gpi=","version="])
    except getopt.GetoptError, err:
        # print help information and exit:
        log.write(str(err)) # will print something like "option -a not recognized"
        print "ERROR:\n"+str(err) 
        gpad2gaf_usage()
        sys.exit(2)
    #set program arguments
    gaf_file=""
    gpad_file=""
    gpi_file=""
    gaf_version="2.0"
    for o, a in opts:
        if o in ("-h", "--help"):
            gpad2gaf_usage()
            sys.exit()
        elif o in ("-g", "--gaf"):gaf_file = a
        elif o in ("-p", "--gpad"):gpad_file = a
        elif o in ("-i", "--gpi"):gpi_file = a
        elif o in ("-v", "--version"):gaf_version = a
        else:
            assert False, "unhandled option"
    #Check if the gpad file exists
    if not os.path.isfile(gpad_file):
        log.write("**********\n\nError: The gpad file: "+gpad_file+\
            "  does not exist\nRun gpad2gaf.py --help to  See program usage")
        print "**********\n\nError: The gpad file: "+gpad_file+\
            " does not exist\nRun gpad2gaf.py --help to  See program usage"
        gpad2gaf_usage()
        sys.exit() 
    if not os.path.isfile(gpi_file):
        log.write("**********\n\nError: The gpi file: "+gpi_file+\
            " does not exist\nRun gpad2gaf.py --help to See program usage")
        print "**********\n\nError: The gpi file: "+gpi_file+\
            " does not exist\nRun gpad2gaf.py --help to See program usage"
        gpad2gaf_usage()
        sys.exit()
    if gaf_file == "":gaf_file=gpad_file+".gaf"
    #Process the gpad
    log.write("\nProcessing GPAD file :"+gpad_file+" and GPI file :"+gpi_file);
    log.write("\nInitiating the converter\n")
    print "\nInitiating the converter\n"
    #
    # Create the expected directory structure and
    # downloads dependencies (gaf-eco-map, mrk_list2.rpt) if needed
    #
    goa_parser.converter_init(log)
    #
    #Check if the containers were initiated properly
    #
    if goa_parser.empty_goref_containers()>0:
        message="**********\n\nSome GO Ontology dependencies were not downladed properly. Check "
        log.write(message+LOCAL_GO_REF_FILE+" and "+LOCAL_GO_REPORT_FILE+" local files\n")
        sys.exit()
    if goa_parser.empty_eco_containers()>0:
        message="**********\n\nSome ECO Evidence dependencies were not downladed properly. Check "
        log.write(message+LOCAL_MGI_GO_ECO_REPORT_FILE+" and "+LOCAL_MAP_FILE+" local files\n")
        sys.exit()

    log_file=log_dir+"/"+os.path.basename(gpad_file)+".log"
    gpad_log=open(log_file,"w")
    i = datetime.now()
    log.write("\nProgram Starts: "+i.strftime('%Y/%m/%d %I:%M:%S %P')+"\n")
    gpad_log.write("\nProgram Starts: "+i.strftime('%Y/%m/%d %I:%M:%S %P')+"\n")
    goa_parser.generateGaf(gaf_file,gpad_file,gpi_file,gpad_log,gaf_version)
    log.close()
    gpad_log.close()

if __name__ == "__main__":
    main()
