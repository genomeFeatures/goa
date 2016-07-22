#!/usr/bin/env python

'''
#
#GAF2GPAD converter: Create a GAF converter that converts a GAF file GPAD/GPI files -
#GAF2GPAD converter
#         Provide GPAD/GPI files on our FTP site in addition to the gene_association.mgi GAF, and the
#         gene_association.mgi_nonmouse file, for these two files;
#         These would be generated every time the respective gafs are generated.
#
# MGI  generates the following GAF files:
#  a. gene_association.mgi GAF
#  b. gene_association.mgi_nonmouse file
#  c. gene_mp_association.mgi file
#
#
#Assumption:
# 1. The order of the first 17 fields of a GAF file does not change
# 2. The order of the first 12 fields of a GPAD file does not change
# 3. The order of the first 9 fields of  a GPI file does not change
#
#
#Author: Lucie Hutchins, Scientific Software Engineer - MGI
#
# History:
#
# lnh	TR11630
#	- created - 01/2015
#
# Usage: 
# a) gaf2gpad.py --help => to display the program help page
# b) gaf2gpad.py --gaf=gaf_file [--gpad=gpad_file] [--gpi=gpi_file] [--mgi]

'''
 
import getopt, sys 
import goa_parser 
import os
from datetime import datetime

#
#gaf2gpad  displays program usage
#
def gaf2gpad_usage():
    print """\
    \n********************************\ngaf2gpad converts a given GAF file to GPAD and GPI files
    \nAssumptions: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa for specs
        1. The order of the first 17/15 fields of a GAF2.0/GAF1.0 file follows goa specs
        2. The order of the first 12 fields of a GPAD file follows goa specs
        3. The order of the first 9 fields of  a GPI file follows goa specs
    \nUsage: 
    Example1: gaf2gpad.py --help  => to display this help page
    Example2: gaf2gpad.py --gaf=gaf_file [--gpad=gpad_file] [--gpi=gpi_file] [--mgi]
    Where:
        --gaf  => <required> specifies the name of the gaf file , gaf_file is the full path to the gaf file
        --gpad => <optional> specifies the name of the output gpad_file(default gaf_file.gpad)
        --gpi  => <optional> specifies the name of the output gpi_file(default gaf_file.gpi)
        --mgi  => <optional> if set, the MGI MRK_List2.rpt report in addition
            to the gaf file will be used to generate the gpi file that includes only MGI set
    \nNote: If you do not provide the name of gpad and gpi result files, the program will create
        these two files in the same directory the input gaf file resides
        with the extension *.gpad and *.gpi respectively
    \n********************************
    """

#
# Main program
#
def main():
    log_dir=goa_parser.GAF_CONVERTER_LOG_BASE
    data_dir=goa_parser.GAF_CONVERTER_DATA_BASE
    if not os.path.isdir('%s' % (log_dir)):
        os.system('mkdir %s' % (log_dir))

    if not os.path.isdir('%s' % (data_dir)):
        os.system('mkdir %s' % (data_dir))
    log_file=log_dir+"/"+os.path.basename(sys.argv[0])+".log"
    log=open(log_file,"w")
    today = datetime.now()
    log.write("Program Starts: "+today.strftime('%Y/%m/%d %I:%M:%S %P'))
    print "Program Starts: "+today.strftime('%Y/%m/%d %I:%M:%S %P')
    log.write("\n")
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:p:i:m", ["help", "gaf=","gpad=","gpi=","mgi"])
    except getopt.GetoptError, err:
        # print help information and exit:
        log.write(str(err)) # will print something like "option -a not recognized"
        print "ERROR:\n"+str(err)  # will print something like "option -a not recognized"
        gaf2gpad_usage()
        sys.exit(2)
    #set program arguments
    gaf_file=""
    gpad_file=""
    gpi_file=""
    filter_mgi = False
    for o, a in opts:
        if o in ("-m","--mgi"):filter_mgi = True
        elif o in ("-h", "--help"):
            gaf2gpad_usage()
            sys.exit()
        elif o in ("-g", "--gaf"):gaf_file = a
        elif o in ("-p", "--gpad"):gpad_file = a
        elif o in ("-i", "--gpi"):gpi_file = a
        else:
            assert False, "unhandled option"
   
    #Check if gaf file exists
    if not os.path.isfile(gaf_file):
        log.write("**********\n\nError: The gaf file: "+gaf_file+" does not exist - See program usage")
        print "**********\n\nError: The gaf file: "+gaf_file+" does not exist - See program usage"
        gaf2gpad_usage()
        sys.exit() 
    if gpad_file == "":gpad_file=gaf_file+".gpad"
    if gpi_file=="":gpi_file=gaf_file+".gpi"
    #
    #Setup the converter log and gaf file process log
    # 
    gaf_log_file=log_dir+"/"+os.path.basename(gaf_file)+".log"
    gaf_log=open(gaf_log_file,"w")
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
        message="**********\n\nSome GO Ontology dependencies were not downloaded properly. Check the logs: "
        log.write(message+LOCAL_GO_REF_FILE+" and "+LOCAL_GO_REPORT_FILE+" local files\n")
        sys.exit()
    if goa_parser.empty_eco_containers()>0:
        message="**********\n\nSome ECO Evidence dependencies were not downloaded properly. Check the logs:"
        log.write(message+LOCAL_MGI_GO_ECO_REPORT_FILE+" and "+LOCAL_MAP_FILE+" local files\n")
        sys.exit() 
    if filter_mgi:gpi_file+=".mgi"
    #
    #Process the gaf file and generate coresponding gpad and gpi files
    #
    today = datetime.now()
    gaf_log.write("\nProgram Starts: "+today.strftime('%Y/%m/%d %I:%M:%S %P')+"\n")
    log.write("\nProcessing GAF file :"+gaf_file)
    print "\nProcessing GAF file :"+gaf_file
    gaf_log.write("\nProcessing GAF file :"+gaf_file)
    goa_parser.generateGpiGpad(gaf_file,gpad_file,gpi_file,gaf_log,filter_mgi)
    today = datetime.now()
    log.write("\nProgram Ends: "+today.strftime('%Y/%m/%d %I:%M:%S %P')+"\n")
    gaf_log.write("\nProgram Ends: "+today.strftime('%Y/%m/%d %I:%M:%S %P')+"\n") 
    print "Program Complete\n"
    log.close()
    gaf_log.close()

if __name__ == "__main__":
    main()
