#!/usr/bin/env python

'''
#
# goa_parser uses goa_specs defined objects to
# handle different GO annotation files.
# The parser also uses global variables defined
# in config.py file
#
#
#Assumption:
# 1. The order of the first 17/15 fields of a GAF file does not change
# 2. The order of the first 11 fields of a GPAD file does not change
# 3. The order of the first 9 fields of  a GPI file does not change
#
# History:
#
# lnh	TR11630
#	- created - 01/2015
#
# 
'''
 
import getopt, sys 
import goa_specs,config 
import os,csv
from datetime import datetime

#
#Set up project contact
#
PNAME=config.PROJECT_NAME
PURL=config.DB_URL
CEMAIL= config.CONTACT_EMAIL

#
#Set up dependencies
#
REMOTE_ECO_MAP_FILE=config.REMOTE_MAP_FILE
REMOTE_GO_REF_FILE=config.REMOTE_GO_REF_FILE
REMOTE_GO_XREF_FILE=config.REMOTE_GO_XREF_FILE
REMOTE_MRK_REPORT_FILE=config.REMOTE_MRK_REPORT_FILE
REMOTE_GO_REPORT_FILE=config.REMOTE_GO_REPORT_FILE
REMOTE_MGI_GO_REPORT_FILE=config.REMOTE_MGI_GO_REPORT_FILE

LOCAL_ECO_MAP_FILE=config.LOCAL_MAP_FILE
LOCAL_GO_REF_FILE=config.LOCAL_GO_REF_FILE
LOCAL_GO_XREF_FILE=config.LOCAL_GO_XREF_FILE
LOCAL_MRK_REPORT_FILE=config.LOCAL_MRK_REPORT_FILE
LOCAL_GO_REPORT_FILE=config.LOCAL_GO_REPORT_FILE
LOCAL_MGI_GO_ECO_REPORT_FILE=config.LOCAL_MGI_GO_ECO_REPORT_FILE

GAF_CONVERTER_LOG_BASE=config.GAF_CONVERTER_LOG_BASE
GAF_CONVERTER_DATA_BASE=config.GAF_CONVERTER_DATA_BASE
GAF_CONVERTER_DATA_TEMP=config.GAF_CONVERTER_DATA_TEMP

#
GPAD_FIELDS=goa_specs.GPAD_FIELDS
GPI_FIELDS=goa_specs.GPI_FIELDS
MRK_FIELDS=goa_specs.MRK_FIELDS
GPAD_FIELDS_LABEL=goa_specs.GPAD_FIELDS_LABEL
GPI_FIELDS_LABEL=goa_specs.GPI_FIELDS_LABEL
GAF_FIELDS_LABEL=goa_specs.GAF_FIELDS_LABEL
GAF_EVIDENCE_CODES=goa_specs.evidence_codes

update_index=config.UPDATE_INDEX
#
#
# set standard field indexes
#
gpad_annot_properties_index=goa_specs.gpad_annot_properties_index
gpad_annot_extension_index=goa_specs.gpad_annot_extension_index
gpad_db_index=goa_specs.gpad_db_index
gpad_db_object_index=goa_specs.gpad_db_object_index
gpad_evidence_index=goa_specs.gpad_evidence_index
gpad_goref_index=goa_specs.gpad_goref_index
gpad_relationship_index=goa_specs.gpad_relationship_index
gpad_taxon_index=goa_specs.gpad_taxon_index
gpad_goid_index=goa_specs.gpad_goid_index


gpi_taxon_index=goa_specs.gpi_taxon_index
gpi_db_index=goa_specs.gpi_db_index
gpi_db_object_index=goa_specs.gpi_db_object_index
gpi_parent_object_id_index=goa_specs.gpi_parent_index
gpi_object_symbol_index= goa_specs.gpi_object_symbol_index
gpi_object_name_index=goa_specs.gpi_object_name_index
gpi_object_syn_index=goa_specs.gpi_object_syn_index
gpi_object_type_index=goa_specs.gpi_object_type_index

mrk_ftype_index=goa_specs.mrk_ftype_index
mrk_object_index=goa_specs.mrk_object_index
mrk_object_name_index=goa_specs.mrk_object_name_index
mrk_object_syn_index=goa_specs.mrk_object_syn_index
mrk_object_id_index=goa_specs.mrk_object_id_index


gaf=goa_specs.gaf()
eco=goa_specs.Eco()
goref=goa_specs.Goref()
goref_collection= goref.COLLECTION
go_id2aspect=goref.GO_ONTOLOGY_MAP
  

#
#
#Displays customized gpad/gpi header            
#
def displayGFile_header(gfh,gaf_header,title,g_fields,g_labels,gf_type,gaf_version):
    gfh.write(title)
    gfh.write("!\n!GAF File Info:\n")
    gfh.write("!\n".join(gaf_header))
    gfh.write("!\n!Project name: %s\n"%(PNAME))
    gfh.write("!Website: %s\n"%(PURL))
    gfh.write("!Contact: %s\n"%(CEMAIL))
    gpad_col="!GPAD - GAF Columns mapping:\n!\n!\tGPAD_column_name\tGAF column #"
    gpi_col="!GPI - GAF Columns mapping:\n!\n!\tGPI_column_name\tGAF column #"
    gaf_col="!GPAD/GPI - GAF Columns mapping:\n!\n!\tGAF Column\tGPAD column\tGPI column"
    i = datetime.now()
    gfh.write("!\n!Date: "+i.strftime('%Y/%m/%d %I:%M:%S %P'))
    gfh.write("\n!\n")
    if gf_type == "gpad":
        gfh.write(gpad_col)
        rows=gaf.getgpad2gaf()
        gfh.write(rows)
    elif gf_type == "gpi":
        gfh.write(gpi_col)
        rows=gaf.getgpi2gaf()
        gfh.write(rows)
    else:
        gfh.write(gaf_col)
        rows= gaf.getgaf2gpadgpi()
        gfh.write(rows) 
        
    gfh.write("\n!\n!GAF - ECO evidence map - dependency :%s"%(REMOTE_ECO_MAP_FILE))
    gfh.write("\n!GO Reference Collection - dependency :%s"%(REMOTE_GO_REF_FILE))
    gfh.write("\n!MGI Marker report - dependency :%s"%(REMOTE_MRK_REPORT_FILE))
    gfh.write("\n!MGI GO report - dependency :%s"%(REMOTE_GO_REPORT_FILE))
    fields=[]
    for i in range(len(g_fields)): fields.append(g_labels[g_fields[i]])
    gfh.write("\n!\n!"+"\t".join(fields)+"\n")

#
#   
#Get GAF version
#
def getGaf_version(gfh):
    gaf_version=""
    gfh.seek(0,0)
    for line in gfh:
        line=line.strip()
        if line.startswith("!"):
           if "gaf-version:" in line: gaf_version=line.replace("gaf-version:","")
        else:
           gfh.seek(0,0) #Again set the pointer to the beginning
           break
    return gaf_version


#
#Get the gaf/gpad/gpi file header
#
def getGFile_header(gfh,header):
    #Set the pointer to the beginning of the file
    gfh.seek(0,0)
    for line in gfh:
        line=line.strip()
        if line.startswith("!"):
           header.append(line)
        else:
           gfh.seek(0,0) #Reset the pointer to the beginning
           break    

#
# Displays corresponding gpad/gpi row from a gaf row
#
def writeGP_row(gaf_row,gf_row_displayed,gfh,g_fields,is_gpad,is_gene_variant,parent_gp_id,evidence_code):
    #set taxon
    taxa=gaf_row[gaf.taxon_index].split("|")
    inter_taxon="";
    taxon=taxa[0]
    if len(taxa)> 1:
        end=len(taxa)
        inter_taxon="".join(taxa[1:end])
    if is_gpad==1: taxon=inter_taxon
    gf_row=[]
    for i in range(len(g_fields)):
        if g_fields[i] in gaf.fields:
           gaf_index=gaf.fields.index(g_fields[i])
           if gaf_index==gaf.taxon_index:
              gf_row.append(taxon)
           else:
              gf_row.append(gaf_row[gaf_index])
        else: #both gpad.Annotation_Properties and gpi.Parent_Object_ID are not gaf fields
           if is_gpad:
              #set GPAD.Annotation_Properties to gaf.evidence code
              if i == gpad_annot_properties_index:
                 gf_row.append(evidence_code)
              else:
                 gf_row.append(" ")
           else:
              if i==gpi_parent_object_id_index:
                 if is_gene_variant:gf_row.append(parent_gp_id)
                 else:gf_row.append(" ")
              else:gf_row.append(" ")

    line="\t".join(gf_row) 
    if not line in gf_row_displayed:
       gf_row_displayed[line]=1
       gfh.write(line+"\n")
#
#Load GPI file into memory - assuming GPI file not too big
#
def loadGpi(gpi_file,gpi_map,log,gpad_map_keys):
    log.write("\nGPI file data log:\n")
    reader = csv.reader(open(gpi_file, 'rb'), dialect='excel-tab')
    row_count=0
    fieldCountMis=0
    missFields=0
    for line in reader:
        if not line[0].startswith("!"):
            gpi_row=[]
            row_count+=1
            if row_count%10000==0: print "%d lines processed"%(row_count)
            for field in line:
                gpi_row.append(field)
            #
            #Data validation step
            #if number of fields does not match, store line in log
            if len(gpi_row) != len(GPI_FIELDS):
                log.write("\t".join(gpi_row)+"\tFields count mismatch:%d - %d" %(len(gpi_row),len(GPI_FIELDS))+"\n")
                fieldCountMis+=1
                continue
            #if missing required fields, store line in log
            field_missing_index=gaf.gpi_has_missing_fields(gpi_row)
            if field_missing_index :
                log.write("\t".join(gpi_row)+"\tFirst required GPI field missing is at index:%d"%(field_missing_index)+"\n")
                missFields+=1
                continue
            db_object_key=gpi_row[gpi_db_index]+":="+gpi_row[gpi_db_object_index]
            parent_id=gpi_row[gpi_parent_object_id_index]
            key= db_object_key+":="+parent_id
            gpi_map[key]=gpi_row
            if db_object_key in gpad_map_keys:
                if parent_id in gpad_map_keys[db_object_key]: gpad_map_keys[db_object_key][parent_id]+=1
                else: gpad_map_keys[db_object_key][parent_id]=1
            else:
                gpad_map_keys[db_object_key]={}
                gpad_map_keys[db_object_key][parent_id]=1
           
           
    log.write("Total number of lines from the gpi file:%s is :%d; total Indexed:%d"%(gpi_file,row_count,len(gpi_map))+"\n")

#
# Generates a GAF file from the specified  gpad and gpi files 
# The default GAF version is 2.0 but the user 
# can specify the GAF version
#
def generateGaf(gaf_file,gpad_file,gpi_file,log,gaf_version):
    gafh=open(gaf_file,"w")
    gpi_map={}         #Loads gpi file into a dictionary indexed by DB:=Object_ID:=parent_id
    gpad_map_keys={}   #Indexes Object_IDs by parent_id - to detect cases where an object is assigned 
                       # to more than one parent
    gaf_header=[] 
    loadGpi(gpi_file,gpi_map,log,gpad_map_keys) #index GPI file  
    type="" 
    title="!gaf-version: %s\n"%(gaf_version)
    mult_parents=0
    log.write("==================\nAnnotated DB:Object_Form_ID  with multiple gene parents:\n")
    for key in gpad_map_keys:
        if len(gpad_map_keys[key])>1:
            parents="|".join(gpad_map_keys[key].keys())
            log.write(key+"\t%d [%s]"%(len(gpad_map_keys[key]),parents)+"\n")
            mult_parents+=1
    #initiate gaf object
    gaf._init(gaf_version)
    displayGFile_header(gafh,gaf_header,title,gaf.fields,GAF_FIELDS_LABEL,type,gaf_version)
    #process the gpad and corresponding gpi
    log.write("=================\nGPAD and GPI data log:\n")
    reader = csv.reader(open(gpad_file, 'rb'), dialect='excel-tab')
    row_count=0
    missFields=0
    object_missing=0
    fieldCountMis=0
    go_missing=0
    eco_w_mult_ev=0
    badDB=0

    for line in reader:
        row_count+=1
        if row_count%10000==0: print "%d lines processed"%(row_count)
        if not line[0].startswith("!"):
            gpad_row=[]
            for field in line:
                gpad_row.append(field)
            #
            #Data validation step
            #if number of fields does not match, store line in log
            if len(gpad_row) != len(GPAD_FIELDS):
                log.write("\t".join(gpad_row)+"\tFields count mismatch:%d - %d" %(len(gpad_row),len(GPAD_FIELDS))+"\n")
                fieldCountMis+=1
                continue
            #if missing required fields, store line in log
            field_missing_index=gaf.gpad_has_missing_fields(gpad_row)
            if field_missing_index :
                log.write("\t".join(gpad_row)+"\tThe first GPAD missing field is at index %d"%(field_missing_index)+"\n")
                missFields+=1
                continue 
            db=gpad_row[gpad_db_index]
            #if bad database name
            #if bad database abbreviation, store line in log
            if db not in goref.GO_DATABASES:
                log.write("\t".join(gpad_row)+"\tThe DB field has an invalid value\n")
                badDB+=1
                continue
            object_id= gpad_row[gpad_db_object_index]
            db_object_key= db+":="+object_id
            object_form_id=""
            if db_object_key not in gpad_map_keys:
                log.write(gpad_row[gpad_db_object_index]+" Not in GPI file\n")
                object_missing +=1
                continue
            #get associated gpi row(s) for this gpad line 
            for parent_id in gpad_map_keys[db_object_key]:
                gaf_row=[]
                for i in range(len(gaf.fields)):gaf_row.append("") #initiate gaf row with empty fields
                gpi_row=gpi_map[db_object_key+":="+parent_id]
                if parent_id:
                    fields=parent_id.split(":")
                    object_form_id=db+":"+object_id
                    if len(fields)>1:
                        db=fields[0]
                        object_id= "".join(fields[1:len(fields)])
                taxon=gpi_row[gpi_taxon_index]
                gaf_row[gaf.db_index]=db
                gaf_row[gaf.object_index]=object_id
                if gpad_row[gpad_taxon_index]:taxon+="|"+gpad_row[gpad_taxon_index]
                gaf_row[gaf.taxon_index]=taxon
                #Set the product_form_id only for gaf 2.0
                if gaf.product_form_id_index in range(len(gaf.fields)):
                    if object_id == object_form_id: object_form_id=""
                    gaf_row[gaf.product_form_id_index]=object_form_id
                if "|" in gpad_row[gpad_relationship_index]:
                    fields= gpad_row[gpad_relationship_index].split("|")
                    gaf_row[gaf.qual_index]="|".join(fields[0:len(fields)-1])
                #
                #Set gaf.evidence_code to gpad.annotation_properties or 
                #get it from eco.ECO_EVIDENCE_MAP - ECO_GOREF2EVIDENCE_MAP
                #
                gaf_row[gaf.annotation_extension_index]=gpad_row[gpad_annot_extension_index]
                eco_code= gpad_row[gpad_evidence_index]
                go_ref=gpad_row[gpad_goref_index]
                evidence_code= eco.getECO_code(go_ref,eco_code,goref.COLLECTION,eco.ECO_GOREF2EVIDENCE_MAP)
                
                gaf_row[gaf.evidence_index]=evidence_code
                if eco_code not in gaf_row[gaf.annotation_extension_index]:
                    gaf_row[gaf.annotation_extension_index]+="|"+eco_code
                #set GAF.Annotation_Extension = 
                #gpad_annot_extension_index=GPAD_FIELDS.index("Annotation_Extension")
                #gaf.annotation_extension_index
                if evidence_code not in GAF_EVIDENCE_CODES:
                   log.write("\t".join(gpad_row)+" -Ambiguous Evidence code:"+evidence_code+"\n") 
                   eco_w_mult_ev+=1

                if gpad_row[gpad_goid_index] not in go_id2aspect:
                    #This means we can't compute the GAF.Aspect field - since we are using the go_id2aspect map
                    log.write(gpad_row[gpad_goid_index]+" not found in MGI GO_terms report\n")
                    go_missing+=1
                    continue
                else:
                    gaf_row[gaf.aspect_index]=go_id2aspect[gpad_row[gpad_goid_index]]
                #now map the rest of empty gaf fields to corresponding gpad and gpi fields
                for i in range(len(gaf.fields)):
                    if i == gaf.db_index : continue
                    elif i == gaf.object_index : continue
                    elif i == gaf.taxon_index: continue
                    elif i == gaf.evidence_index: continue
                    elif i== gaf.aspect_index: continue
                    elif i== gaf.qual_index:continue
                    elif i== gaf.product_form_id_index: continue
                    elif i== gaf.annotation_extension_index: continue
                    if not gaf_row[i]:
                        if gaf.fields[i] in GPAD_FIELDS:
                            gpad_index=GPAD_FIELDS.index(gaf.fields[i])
                            gaf_row[i]=gpad_row[gpad_index]
                        elif gaf.fields[i] in GPI_FIELDS:
                            gpi_index=GPI_FIELDS.index(gaf.fields[i])
                            gaf_row[i]=gpi_row[gpi_index]
              
                # Display row 
                gafh.write("\t".join(gaf_row)+"\n")
    #              
    log.write("\nTotal rows in GPAD file: %d " %(row_count))
    log.write("\nTotal rows with Fields count mismatch : %d " %(fieldCountMis))
    log.write("\nTotal rows with missing required fields : %d " %(missFields))
    log.write("\nTotal Object_id in GPAD but not in GPI : %d " %(object_missing))
    log.write("\nTotal Annotated DB:Object_Form_ID with multiple gene parents : %d " %(mult_parents))
    log.write("\nTotal rows with ECO code mapping to multiple base Evidence code : %d " %(eco_w_mult_ev))
    log.write("\nTotal rows with invalid DB name : %d " %(badDB))
    log.write("\nTotal GPAD rows with GO_ID not in MGI GO_terms report : %d \n" %(go_missing))

    gafh.close()

#
# Create the expected directory structure and
# downloads dependencies (gaf-eco-map, mrk_list2.rpt) if needed
#
def setup_depends(log):
    try:
        right_now = datetime.now().strftime("%Y%m%d")
        if not os.path.isdir('%s' % (GAF_CONVERTER_LOG_BASE)):
            os.system('mkdir %s' % (GAF_CONVERTER_LOG_BASE))
        if not os.path.isdir('%s' % (GAF_CONVERTER_DATA_BASE)):
            os.system('mkdir %s' % (GAF_CONVERTER_DATA_BASE))
        if not os.path.isdir('%s' % (GAF_CONVERTER_DATA_TEMP)):
            os.system('mkdir %s' % (GAF_CONVERTER_DATA_TEMP))
        #
        # Download reports as needed
        #
        if not os.path.isfile('%s' % (LOCAL_MRK_REPORT_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_MRK_REPORT_FILE,REMOTE_MRK_REPORT_FILE))
        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_MRK_REPORT_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_MRK_REPORT_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_MRK_REPORT_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_MRK_REPORT_FILE,REMOTE_MRK_REPORT_FILE))

        if not os.path.isfile('%s' % (LOCAL_ECO_MAP_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_ECO_MAP_FILE,REMOTE_ECO_MAP_FILE))
        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_ECO_MAP_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_ECO_MAP_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_ECO_MAP_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_ECO_MAP_FILE,REMOTE_ECO_MAP_FILE))

        if not os.path.isfile('%s' % (LOCAL_GO_REF_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_GO_REF_FILE,REMOTE_GO_REF_FILE))
        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_GO_REF_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_GO_REF_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_GO_REF_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_GO_REF_FILE,REMOTE_GO_REF_FILE))

        if not os.path.isfile('%s' % (LOCAL_GO_XREF_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_GO_XREF_FILE,REMOTE_GO_XREF_FILE))
        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_GO_XREF_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_GO_XREF_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_GO_XREF_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_GO_XREF_FILE,REMOTE_GO_XREF_FILE))

        if not os.path.isfile('%s' % (LOCAL_GO_REPORT_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_GO_REPORT_FILE,REMOTE_GO_REPORT_FILE))
        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_GO_REPORT_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_GO_REPORT_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_GO_REPORT_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_GO_REPORT_FILE,REMOTE_GO_REPORT_FILE))
        if not os.path.isfile('%s' % (LOCAL_MGI_GO_ECO_REPORT_FILE)):
            os.system('wget -q -O %s "%s"' %(LOCAL_MGI_GO_ECO_REPORT_FILE,REMOTE_MGI_GO_REPORT_FILE))
            #Only get the three fields from mgi_eco (evidenceCode, mgiRef, ecoCode)
            temp_file=LOCAL_MGI_GO_ECO_REPORT_FILE+".temp"
            os.system("cut -f '6 8 12' %s > %s"%(LOCAL_MGI_GO_ECO_REPORT_FILE,temp_file))
            if os.path.isfile('%s' % (temp_file)):
                os.system("mv %s %s"%(temp_file,LOCAL_MGI_GO_ECO_REPORT_FILE))

        else:
            local_f_mtime = datetime.fromtimestamp(os.path.getmtime(LOCAL_MGI_GO_ECO_REPORT_FILE))
            local_f_mtime=local_f_mtime.strftime("%Y%m%d")
            diff_time= int(right_now) - int(local_f_mtime)
            if diff_time > update_index[LOCAL_MGI_GO_ECO_REPORT_FILE]:
                os.system('cp -Rp %s %s' %(LOCAL_MGI_GO_ECO_REPORT_FILE,GAF_CONVERTER_DATA_TEMP))
                os.system('wget -q -O %s "%s"' %(LOCAL_MGI_GO_ECO_REPORT_FILE,REMOTE_MGI_GO_REPORT_FILE))
                #Only get the three fields from mgi_eco (evidenceCode, mgiRef, ecoCode)
                temp_file=LOCAL_MGI_GO_ECO_REPORT_FILE+".temp"
                os.system("cut -f '6 8 12' %s > %s"%(LOCAL_MGI_GO_ECO_REPORT_FILE,temp_file))
                if os.path.isfile('%s' % (temp_file)):
                    os.system("mv %s %s"%(temp_file,LOCAL_MGI_GO_ECO_REPORT_FILE))
    except OSError, e:
        errStr= "Error %d: %s" % (e.args[0], e.args[1])
        log.write("Converter Initiation failed:" + errStr+"\n")
        #restore the files
        os.system("mv %s/* %s"%(GAF_CONVERTER_DATA_TEMP,GAF_CONVERTER_DATA_BASE))
        sys.exit(1)
    return 0
#
# Initiate converter
#  
def converter_init(log):
    setup_depends(log)
    #Load dictionaries 
    map_log_file=GAF_CONVERTER_LOG_BASE+"/maps.log"
    map_log=open(map_log_file,"w")
    goref._init(LOCAL_GO_REF_FILE,LOCAL_GO_REPORT_FILE,LOCAL_GO_XREF_FILE,map_log)
    eco._init(LOCAL_ECO_MAP_FILE,LOCAL_MGI_GO_ECO_REPORT_FILE,map_log)

#
#returns true if one of the key maps is empty
#
def empty_goref_containers():
    is_empty=0
    if len(goref.GO_ONTOLOGY_MAP)<=0:is_empty=1
    if len(goref.COLLECTION)<=0:is_empty=1
    return  is_empty
#
#returns true if one of the key maps is empty
#
def empty_eco_containers():
    is_empty=0
    if len(eco.GAF_ECO_MAP)<=0: is_empty=1
    if len(eco.ECO_EVIDENCE_MAP)<=0: is_empty=1
    if len(eco.MGI_ECO_MAP)<=0: is_empty=1
    return  is_empty
#
#Generates gpad and gpi files 
#from a gaf file using goa specification
#Both the gpi and gpad file will have the same
#version as the input gaf file
#
def generateGpiGpad(gaf_file,gpad_file,gpi_file,log,filtermgi):
    gafh=open(gaf_file)
    gpad=open(gpad_file,"w")
    gpi=open(gpi_file,"w")
    gaf_header=[]
    gaf_version=getGaf_version(gafh)
    #initiate gaf object
    gaf._init(gaf_version)
    getGFile_header(gafh,gaf_header)
    gafh.close()
    gpad_row_displayed={} #structure to filter gpad duplicate rows if any
    gpi_row_displayed={}  #structure to filter gpi duplicate rows if any
    feature_type_map={}
    protein_map={}        #stores protein-gene mapping
    totalEmptyProt=0      #Keeps track of entries where the protein field is empty
    fieldCountMis=0       #Keeps track of entrie where the field count !=17
    missFields=0          #Keeps track of entries where at least one required field is empty
    badEvCode=0           #Keeps track of entries where the provided Evidence code is invalid
    badDB=0               #Keeps track of entries where the provided DB is invalid
    badAspect=0           #Keeps track of entries where the GAF Aspect field is invalid
    badIsoform=0           #Keeps track of entries where the GAF Gene form id field is invalid
    row_count=0
    req_gaf_fields={}
    gpad_title="!gpad-version: %s\n"%(gaf_version)
    gpi_title="!gpi-version: %s\n"%(gaf_version)
    type="gpad"
    displayGFile_header(gpad,gaf_header,gpad_title,GPAD_FIELDS,GPAD_FIELDS_LABEL,type,gaf_version)
    type="gpi"
    displayGFile_header(gpi,gaf_header,gpi_title,GPI_FIELDS,GPI_FIELDS_LABEL,type,gaf_version)
    
    log.write("GAF file data log:\n")
    reader = csv.reader(open(gaf_file, 'rb'), dialect='excel-tab')
    for line in reader:
        row_count+=1
        if row_count%10000==0: print "%d lines processed"%(row_count)
        if not line[0].startswith("!"):
           gaf_row=[]
           for field in line:
               gaf_row.append(field)
           #if number of fields does not match, store line in log
           if len(gaf_row) != len(gaf.fields): 
              log.write("\t".join(gaf_row)+"\tFields count mismatch:%d - %d" %(len(gaf_row),len(gaf.fields))+"\n")
              fieldCountMis+=1
              continue       
           #if missing required fields, store line in log
           field_missing_index=gaf.has_missing_fields(gaf_row)
           if field_missing_index :
              log.write("\t".join(gaf_row)+"\tThe first GAF missing field is at index %d"%(field_missing_index)+"\n")
              missFields+=1
              continue
           #if bad database abbreviation, store line in log
           if gaf_row[0] not in goref.GO_DATABASES:
              log.write("\t".join(gaf_row)+"\tThe DB field has an invalid value\n")
              badDB+=1
              continue 
           #
           # Set GPAD ECO code  
           # eco_code= evidence2eco(gaf.evidence_code+gaf.goref)
           #
           go_ref=gaf_row[gaf.goref_index]
           evidence_code=gaf_row[gaf.evidence_index]
           # if invalid evidence code, skip and store line in log
           eco_code= eco.getECO_code(go_ref,evidence_code,goref.COLLECTION,eco.GAF_ECO_MAP)
           if not eco_code:
              log.write("\t".join(gaf_row)+"\tBad Evidence code\n")
              badEvCode +=1
              continue
           # Skip if bad Aspect field
           relationship=gaf.getGPAD_relationship(gaf_row)
           if not relationship:
              log.write("\t".join(gaf_row)+"\tBad Aspect field\n")
              badAspect+=1
              continue
           #Initiate GPI parent id in case this is a gene variant
           parent_gp_id=gaf.getParent_gp_id(gaf_row)
           #set if this is an annotation of specific variants of the gene or gene product
           is_g_variant=gaf.is_gene_variant(gaf_row)
           #
           #Check if protein object type and set BD, and DB_Object_ID as needed - only for gaf 2.*
           #Will overwrite the DB and DB_Object_ID :
           if  is_g_variant:
               protein=gaf_row[gaf.product_form_id_index] #GAF - field 17
               #if "|" in protein or len(protein.split(":"))>2:
               if "|" in protein:
                   log.write("\t".join(gaf_row)+"\t --- Bad Gene Product Form ID field\n")
                   badIsoform+=1
                   continue
               totalEmptyProt+=gaf.setProtein(protein_map,gaf_row)
           #overwrite the evidence code in the gaf line with the translated ECO code
           gaf_row[gaf.evidence_index]=eco_code
           #Store the ECO code in  gaf.annotation_extension field used in evidence property note
           if gaf_row[gaf.annotation_extension_index]:
               gaf_row[gaf.annotation_extension_index]+="|"+eco_code
           else:
               gaf_row[gaf.annotation_extension_index]=eco_code

           #Keep tally of feature types
           feature_type=gaf_row[gaf.feature_type_index]
           if feature_type in feature_type_map:
              feature_type_map[feature_type]+=1
           else:
              feature_type_map[feature_type]=1
           #
           # Set gpad relationship field
           gaf_row[gaf.qual_index]= gaf.getGPAD_relationship(gaf_row)
           #Display gpad  row - filter duplicates
           is_gpad=1
           writeGP_row(gaf_row,gpad_row_displayed,gpad,GPAD_FIELDS,is_gpad,is_g_variant,parent_gp_id,evidence_code)
           #
           #Display gpi  row - filter duplicates
           is_gpad=0
           if not filtermgi:
              writeGP_row(gaf_row,gpi_row_displayed,gpi,GPI_FIELDS,is_gpad,is_g_variant,parent_gp_id,evidence_code)
    mrkCountMis=0
    if filtermgi:
       mrkCountMis+=generateGPI_mgi(gpi,protein_map,log)
    log.write("***************\nGAF file features tally\n*************\n\n")
    log.write("feature_type\tRow Count\n")
    for feature_type in feature_type_map:
        log.write(feature_type+"\t%d"%(feature_type_map[feature_type])+"\n")
   
    log.write("\nTotal rows with empty Protein field: %d " %(totalEmptyProt))
    log.write("\nTotal rows with invalid Evidence Code : %d " %(badEvCode))
    log.write("\nTotal Protein ID count: %d " %(len(protein_map)))
    log.write("\nTotal GAF rows with field count !=%d : %d" %(len(gaf.fields),fieldCountMis))
    log.write("\nTotal GAF rows with invalid DB : %d" %(badDB))
    log.write("\nTotal GAF rows with invalid Aspect : %d" %(badAspect))
    log.write("\nTotal GAF rows with invalid Gene form ID : %d" %(badIsoform))
    if filtermgi: 
       log.write("\nTotal MGI Marker report: rows with field count !=%d : %d" %(len(MRK_FIELDS),mrkCountMis))
    log.write("\nProgram Complete")
    gpad.close()
    gpi.close()

#
#generates the GPI file using MGI MRK_List2.rpt and 
#gpi_type as defined by Mary Dolan
#
def generateGPI_mgi(gpi,protein_map,log):
     reader = csv.reader(open(LOCAL_MGI_REPORT_FILE, 'rb'), dialect='excel-tab')
     log.write("\n************\nData log for MGI marker report :"+LOCAL_MGI_REPORT_FILE+"\n************\n") 
     total=0
     mgi2symbol={} #maps mgi accession id to marker symbol
     mgi2name={}   #maps mgi accession id to marker name
     mgi2syn={}    #maps mgi accession id to marker synonym
     gpi_row_displayed={}
     taxon = "taxon:10090"
     db="MGI"
     for line in reader:
         if not line[3].startswith("genome"):
            mrk_row=[]
            for field in line:
                mrk_row.append(field)
            if len(mrk_row) != len(MRK_FIELDS):
               log.write(line+" -------- %d" %(len(mrk_row))+"\n")
               total+=1
               continue
            mgi_id=mrk_row[mrk_object_id_index]
            mgi2symbol[mgi_id]=mrk_row[mrk_object_index]
            mgi2name[mgi_id]=mrk_row[mrk_object_name_index]
            mgi2syn[mgi_id]=mrk_row[mrk_object_syn_index]
            ftype=mrk_row[mrk_ftype_index]
            if ftype in goa_specs.gpi_type:
               mrk_row[mrk_ftype_index]=goa_specs.gpi_type[ftype]
               writeGPI_MGI_row(mrk_row,gpi_row_displayed,gpi,taxon,db) 
     #Now display proteins from protein_map
     writeGPI_PROT_row(protein_map,gpi_row_displayed,gpi,mgi2symbol,mgi2name,mgi2syn,taxon,log,db) 
     return total
#
# Write MGI genes
#
def writeGPI_MGI_row(mrk_row,gpi_row_displayed,gpi,taxon,db):
    gpi_row=[]
    for i in range(len(GPI_FIELDS)):
        if GPI_FIELDS[i] in MRK_FIELDS:
           mrk_index=MRK_FIELDS.index(GPI_FIELDS[i])
           gpi_row.append(mrk_row[mrk_index])
        elif i == gpi_taxon_index:
           gpi_row.append(taxon)
        elif i ==gpi_db_index:
           gpi_row.append(db) 
        else:
           gpi_row.append("")
    line="\t".join(gpi_row)
    if not line in gpi_row_displayed:
       gpi_row_displayed[line]=1
       gpi.write(line+"\n")
#
# Write assocated protein 
#
def writeGPI_PROT_row(protein_map,gpi_row_displayed,gpi,mgi2symbol,mgi2name,mgi2syn,taxon,log,db):
    for protein in protein_map:
        mgi_id=protein_map[protein]
        fields=protein.split(":")
        prot_db=fields[0]
        protein_id=fields[1]
        if not mgi_id in mgi2symbol:
            log.write(protein_id+"\n"+mgi_id+" - not in MGI\n")
        else:
            gpi_row=[]
            for i in range(len(GPI_FIELDS)):
                if i == gpi_db_index:
                    gpi_row.append(prot_db)
                elif i==gpi_db_object_index:
                    gpi_row.append(protein_id)
                elif i== gpi_object_symbol_index:
                    gpi_row.append(mgi2symbol[mgi_id])
                elif i == gpi_object_name_index:
                    gpi_row.append(mgi2name[mgi_id])
                elif i == gpi_object_syn_index:
                    gpi_row.append(mgi2syn[mgi_id])
                elif i== gpi_object_type_index:
                    gpi_row.append("protein")
                elif i == gpi_taxon_index:
                    gpi_row.append(taxon)
                elif i== gpi_parent_object_id_index:
                    gpi_row.append(db+":"+mgi_id)
                else:
                    gpi_row.append("")
            line="\t".join(gpi_row)
            if not line in gpi_row_displayed:
                gpi_row_displayed[line]=1
                gpi.write(line+"\n")
               
#
