#!/usr/bin/env python
import getopt, sys
import os,csv
import re

#
# Create data containers for fields mapping in the GAF, GPAD,GPI, and MRK report
# files 
# Fields name are normalized (simmilar fields are assigned the same label) 
# and field labels created
#
#Assumption:
# 1. The order of the first 17/15 fields of a GAF 2.0/1.0 file does not change
# 2. The order of the first 11 fields of a GPAD file does not change
# 3. The order of the first 9 fields of  a GPI file does not change
# 4. The order of the first 12 fields of the MRK_List2.rpt  file does not change
#
# History:
#
# lnh	TR11630
#	- created - 01/2015
#
#
#GO ID GO id is a root term:  GO:0003674, GO:0008150, GO:0005575
#
GOID_IS_ROOT=['GO:0003674', 'GO:0008150', 'GO:0005575']

#GAF 2.0 fields setup - at least the following first 17 fields
#From http://geneontology.org/page/go-annotation-file-gaf-format-20
GAF_FIELDS=['DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'Qualifier',
            'GO_ID',
            'DB_Reference',
            'Evidence_Code',
            'With',
            'Aspect',
            'DB_Object_Name',
            'DB_Object_Synonym',
            'DB_Object_Type',
            'Taxon',
            'Date',
            'Assigned_By',
            'Annotation_Extension',
            'Gene_Product_Form_ID']

#GAF 1.0 fields setup - at least the following first 15 fields
#http://geneontology.org/page/go-annotation-file-gaf-format-10
#
GAF1_FIELDS=['DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'Qualifier',
            'GO_ID',
            'DB_Reference',
            'Evidence_Code',
            'With',
            'Aspect',
            'DB_Object_Name',
            'DB_Object_Synonym',
            'DB_Object_Type',
            'Taxon',
            'Date',
            'Assigned_By']


REQUIRED_GAF_FIELDS=[
            'DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'GO_ID',
            'DB_Reference',
            'Evidence_Code',
            'Aspect',
            'DB_Object_Type',
            'Taxon',
            'Date',
            'Assigned_By']

GAF_FIELDS_LABEL={'DB':'DB',
             'DB_Object_ID':'DB Object ID',
             'DB_Object_Symbol':'DB Object Symbol',
             'Qualifier':'Qualifier',
             'GO_ID':'GO ID',
             'DB_Reference':'DB:ReferenceDB(|DB:Reference)',
             'Evidence_Code':'Evidence code',
             'With':'With (or) From',
             'Aspect':'Aspect',
             'DB_Object_Name':'DB Object Name',
             'DB_Object_Synonym':'DB Object Synonym(|Synonym)',
             'DB_Object_Type':'DB Object Type',
             'Taxon':'Taxon(|taxon)',
             'Date':'Date',
             'Assigned_By':'Assigned By',
             'Annotation_Extension':'Annotation Extension',
             'Gene_Product_Form_ID':'Gene Product Form ID'}
#
#GAP fields setup - at least the following first 12 fields
#
GPAD_FIELDS=['DB',
             'DB_Object_ID',
             'Qualifier',
             'GO_ID',
             'DB_Reference',
             'Evidence_Code',
             'With',
             'Taxon',
             'Date',
             'Assigned_By',
             'Annotation_Extension',
             'Annotation_Properties']

GPAD_FIELDS_LABEL={'DB':'DB',
             'DB_Object_ID':'DB_Object_ID',
             'Qualifier':'Relationship',
             'GO_ID':'GO_ID',
             'DB_Reference':'DB:Reference(s)',
             'Evidence_Code':'Evidence Code',
             'With':'With (or) From',
             'Taxon':'Interacting taxon ID',
             'Date':'Date',
             'Assigned_By':'Assigned_by',
             'Annotation_Extension':'Annotation XP','Annotation_Properties':'Annotation Properties'}

REQUIRED_GPAD_FIELDS=[
            'DB',
            'DB_Object_ID',
            'Qualifier',
            'GO_ID',
            'DB_Reference',
            'Evidence_Code',
            'Date',
            'Assigned_By']
#
#GPI fields setup - at least the following first  9 fields
#I ran into a case where a gene product was assigned to
# more than one gene on different dates. 
# Without the date field in the GPI, it was hard to
# reconstruct the gaf row using the gpad and gpi files
# in those cases 
#
GPI_FIELDS=['DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'DB_Object_Name',
            'DB_Object_Synonym',
            'DB_Object_Type',
            'Taxon',
            'Parent_Object_ID',
            'DB_Xref'] #,'Date']

GPI_FIELDS_LABEL={'DB':'DB',
                  'DB_Object_ID':'DB_Object_ID',
                  'DB_Object_Symbol':'DB_Object_Symbol',
                  'DB_Object_Name':'DB_Object_Name',
                  'DB_Object_Synonym':'DB_Object_Synonym(s)',
                  'DB_Object_Type':'DB_Object_Type',
                  'Taxon':'Taxon',
                  'Parent_Object_ID':'Parent_Object_ID',
                  'DB_Xref':'DB_Xref(s)'} #,'Date':'Date'}

REQUIRED_GPI_FIELDS=[
            'DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'DB_Object_Type',
            'Taxon'
]
#
#
#pub/reports/MRK_List2.rpt report fields setup - at least the following first   fields
#
MRK_FIELDS=['DB_Object_ID',
            'chr',
            'cm',
            'start',
            'end',
            'strand',
            'DB_Object_Symbol',
            'status',
            'DB_Object_Name',
            'mrk_type',
            'DB_Object_Type',
            'DB_Object_Synonym']

MRK_FIELDS_MAP={'MGI Accession ID':'DB_Object_ID',
                'Chr':'chr',
                'cM Position':'cm',
                'genome coordinate start':'start',
                'genome coordinate end':'end',
                'strand':'strand',
                'Marker Symbol':'DB_Object_Symbol',
                'Status':'status',
                'Marker Name':'DB_Object_Name',
                'Marker Type':'mrk_type',
                'Feature Type':'DB_Object_Type',    
                'Marker Synonyms (pipe-separated)':'DB_Object_Synonym'}
#
#Predefined feature types in the GPI FILE - as defined by Mary Dolan
gpi_type = {
        "protein coding gene":"protein coding gene",
        "non-coding RNA gene":"ncRNA",
        "rRNA gene":"rRNA",
        "tRNA gene":"tRNA",
        "snRNA gene":"snRNA",
        "snoRNA gene":"snoRNA",
        "miRNA gene":"miRNA",
        "scRNA gene":"scRNA",
        "lincRNA gene":"lincRNA",
        "SRP RNA gene":"SRP_RNA",
        "RNase P RNA gene":"RNase_P_RNA",
        "RNase MRP RNA gene":"RNase_MRP_RNA",
        "telomerase RNA gene":"telomerase_RNA"
        }

#GO Reference Collection
#from: http://www.geneontology.org/cgi-bin/references.cgi
go_ref_collection={}
'''
       "MGI:2152098":"GO_REF:0000002",
       "J:72247":"GO_REF:0000002",
       "ZFIN:ZDB-PUB-020724-1":"GO_REF:0000002",
       "FB:FBrf0174215":"GO_REF:0000002",
       "dictyBase_REF:10157":"GO_REF:0000002",
       "SGD_REF:S000124036":"GO_REF:0000002",
       "MGI:2152096":"GO_REF:0000003",
       "J:72245":"GO_REF:0000003",
       "ZFIN:ZDB-PUB-031118-3":"GO_REF:0000003",
       "SGD_REF:S000124037":"GO_REF:0000003",
       "MGI:1354194":"GO_REF:0000004",
       "J:60000":"GO_REF:0000004",
       "ZFIN:ZDB-PUB-020723-1":"GO_REF:0000004",
       "SGD_REF:S000124038":"GO_REF:0000004",
       "MGI:2152097":"GO_REF:0000006",
       "J:72246":"GO_REF:0000006",
       "MGI:2154458":"GO_REF:0000008",
       "J:73065":"GO_REF:0000008",
       "MGI:1347124":"GO_REF:0000010",
       "J:56000":"GO_REF:0000010",
       "MGI:4459044":"GO_REF:0000033",
       "J:161428":"GO_REF:0000033",
       "SGD_REF:S000146947":"GO_REF:0000033",
       "J:104715":"GO_REF:0000024",
       "J:73065":"GO_REF:0000024",
       "MGI:2156816":"GO_REF:0000015",
       "RGD:1598407":"GO_REF:0000015",
   }
'''

qualifier_order = ['not',
	           'colocalizes_with',
	           'contributes_to']
#hash that maps relation to ontology and (where appropriate) the GAF 2.0 equivalent
#compute the GAF qualifier field using GPAD relationship field
rln2qualifier = {
	'contributes_to':'contributes_to',
	'colocalizes_with':'colocalizes_with',
	'not':'not'}

# hash containing the default relations for each ontology
# maps GAF Aspect field to
aspect2rln = {
	'P':'actively_participates_in',
	'F':'actively_participates_in',
	'C':'part_of',
	'default':'annotated_to'}

# hash that maps GPAD relation to GAF ontology and (where appropriate) the GAF 2.0 equivalent
relations = {
        # cellular_component
	'part_of-Aspect-C':'C',
	'part_of-Qualifier-C':'',
	'colocalizes_with-Aspect-C':'C',
	'colocalizes_with-Qualifier-C':'colocalizes_with',
	'active_in-Aspect-C':'C' ,
	'transported_by-Aspect-C':'C',
	'posttranslationally_modified_in-Aspect-C':'C',
	'located_in_other_organism-Aspect-C':'C',
	'located_in_host-Aspect-C':'C',
	'member_of-Aspect-C':'C',
	'intrinsic_to-Aspect-C':'C',
	'extrinsic_to-Aspect-C':'C',
	'spans-Aspect-C':'C',
	'partially_spans-Aspect-C':'C',
        # molecular_function
        'actively_participates_in-Aspect-F':'F',  
        'actively_participates_in-Qualifier-F':'',
	'contributes_to-Aspect-F':'F',
	'contributes_to-Qualifier-F':'contributes_to',
	'functions_in_other_organism-Aspect-F':'F',
	'functions_in_host-Aspect-F':'F',
	'substrate_of-Aspect-F':'F',
        # biological_process
        'actively_participates_in-Aspect-P':'P',   
}
aspect2ontology={
     'C':"Cellular Component",
     'F':"Molecular Funtion",
     'P':"Biological Process"
}

ontology2aspect={
    "Cellular Component":"C",
    "Molecular Function":"F",
    "Biological Process":"P"
}
#
#Default list of gaf evidence codes
# This is updated at runtime if new code is 
# found from goa eco_map file
#
evidence_codes=[
  #Experimental Evidence codes 
  "EXP",
  "IDA",
  "IPI",
  "IMP",
  "IGI",
  "IEP",
  #Computational Analysis evidence codes 
  "ISS",
  "ISO",
  "ISA",
  "ISM",
  "IGC",
  "IBA",
  "IBD",
  "IKR",
  "IRD",
  "RCA",
  #Author Statement evidence codes 
  "NAS",
  "TAS",
  #Curatorial Statement codes 
  "IC",
  "ND",
  #Automatically-Assigned evidence code
  "IEA"]

#global variables
gpad_db_index=GPAD_FIELDS.index("DB")
gpad_db_object_index=GPAD_FIELDS.index("DB_Object_ID")
gpad_evidence_index=GPAD_FIELDS.index("Evidence_Code")
gpad_goref_index=GPAD_FIELDS.index("DB_Reference")
gpad_relationship_index=GPAD_FIELDS.index("Qualifier")
gpad_taxon_index=GPAD_FIELDS.index("Taxon")
gpad_goid_index=GPAD_FIELDS.index("GO_ID")
gpad_annot_properties_index=GPAD_FIELDS.index("Annotation_Properties")
gpad_annot_extension_index=GPAD_FIELDS.index("Annotation_Extension")

gpi_db_index=GPI_FIELDS.index("DB")
gpi_db_object_index=GPI_FIELDS.index("DB_Object_ID")
gpi_parent_index=GPI_FIELDS.index("Parent_Object_ID")
gpi_taxon_index=GPI_FIELDS.index("Taxon")
gpi_object_symbol_index=GPI_FIELDS.index("DB_Object_Symbol")
gpi_object_name_index=GPI_FIELDS.index("DB_Object_Name")
gpi_object_syn_index=GPI_FIELDS.index("DB_Object_Synonym")
gpi_object_type_index=GPI_FIELDS.index("DB_Object_Type")

mrk_ftype_index=MRK_FIELDS.index("DB_Object_Type")
mrk_object_index=MRK_FIELDS.index("DB_Object_Symbol")
mrk_object_name_index=MRK_FIELDS.index("DB_Object_Name")
mrk_object_syn_index=MRK_FIELDS.index("DB_Object_Synonym")
mrk_object_id_index=MRK_FIELDS.index("DB_Object_ID")

#############################################################
#
# Classes Definition
#
#
##############################################################
# goref class maps: 
# A) GO Reference Collection - external ID mapping
#    1) External ids such as J#, MGI ids, pubmed ids to GO_Ref
#    2) Alternate GO_Ref ids to GO_Ref
# B) GO database cross-references: maps known database ids to database name
#
# Source: http://www.geneontology.org/cgi-bin/references.cgi
# http://www.geneontology.org/doc/GO.references
# http://geneontology.org/doc/GO.xrf_abbs
class Goref:
    #####################
    # Global variables ##
    #####################
    ALT_ID_MAP={}       #maps alternate GO_Ref ids of a given GO_Ref number
    COLLECTION={}       #maps external ids such as J#, MGI ids, pubmed ids... to goref number
    GO_ONTOLOGY_MAP={}  #maps GO_ID to ontology term (C=> Cellular Component, 
                        #F=> Molecular Function, P=> Biological Process) 
    GO_DATABASES={}     #maps of current GO databases 

    ##########################
    # class constructor    ###
    ##########################
    def _init(self,local_goref_file,local_go_file,local_goxref_file,log):
       #
       # Index GO_ID TO ontology
       if os.path.isfile('%s' % (local_go_file)):
          efile=open(local_go_file)
          lines=efile.readlines()
          for line in lines:
              line=line.rstrip()
              fields= line.split("\t")
              if len(fields) == 3 :
                 ontology=fields[0].strip()
                 go_id=fields[1].strip()
                 Goref.GO_ONTOLOGY_MAP[go_id]=ontology2aspect[ontology]

          log.write("\n====================\nGO_ID - GAF.Aspect mapping\n")
          for go_id in Goref.GO_ONTOLOGY_MAP:
              log.write(go_id+"\t"+Goref.GO_ONTOLOGY_MAP[go_id]+"\n")

       if os.path.isfile('%s' % (local_goref_file)):
          efile=open(local_goref_file)
          #remove header
          for line in efile:
              line=line.rstrip()
              if not line.startswith("!"): break
          goref_id=""
          for line in efile:
              if not line: continue
              line=line.strip()
              fields=line.split(":")
              if len(fields)>=2:
                 label=fields[0]
                 value=line.replace(label+":","").strip()
                 if "citation" in label: label="external_accession"
                 elif "abstract" in label: continue
                 elif "year" in label: continue
                 elif "authors" in label: continue
                 elif "title" in label: continue
                 elif "comment" in label: continue
                 if "go_ref_id" in label: goref_id=value
                 else:
                     if "external_accession" in label:
                         Goref.COLLECTION[value]=goref_id
                     elif "alt_id" in label:
                         Goref.ALT_ID_MAP[value]=goref_id

       if os.path.isfile('%s' % (local_goxref_file)):
           efile=open(local_goxref_file)
           for line in efile:
               line=line.rstrip()
               if not line.startswith("!"): break
           db_abbrev=""
           db_name=""
           for line in efile:
               if "abbreviation:" in line:
                   db_abbrev=line.replace("abbreviation:","").strip()   
               elif "database:" in line:
                   db_name=line.replace("database:","").strip()
                   Goref.GO_DATABASES[db_abbrev]=db_name
                   db_abbrev=""
                   db_name=""
               else: continue
          #remove header
       #Now display map content
       log.write("\n====================\nGO_REF External Accession IDs mapping\n")
       for extaccid in Goref.COLLECTION:
           log.write(extaccid+"\t"+Goref.COLLECTION[extaccid]+"\n")
       log.write("\n====================\nGO_REF Alternate IDs \n")
       for altid in Goref.ALT_ID_MAP:
           log.write(altid+"\t"+Goref.ALT_ID_MAP[altid]+"\n") 
       log.write("\n====================\nGO database cross-reference \n")
       for db_abbrev in Goref.GO_DATABASES:
           log.write(db_abbrev+"\t"+Goref.GO_DATABASES[db_abbrev]+"\n")

###################################################################
# eco class maps:
# 1) ECO to Evidence code
# 2) Evidence code to ECO
#
# Source: http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt
##################################################################
 
class Eco:
 
    #####################
    # Global variables ##
    #####################
    GAF_ECO_MAP={}       #maps evidenceCode-goref to ECO code
    ECO_EVIDENCE_MAP={}  #maps eco to evidence code (gpad ECO code to gaf evidence code) 
    MGI_ECO_MAP={}       #maps evidenceCode-MGIRef accession to ECO code  
    ECO_GOREF2EVIDENCE_MAP={} #maps ECO-GOREF to GAF Evidence code

    ##########################
    # class constructor    ### 
    def _init(self,local_eco_file,local_mgi_eco_file,log):
    
        #
        #Index gaf-eco mapping from genontology public file 
        #
        ecoref_ev_map={}
        if os.path.isfile('%s' % (local_eco_file)):
            efile=open(local_eco_file)
            lines=efile.readlines()
            for line in lines:
                line=line.rstrip()
                if not line.startswith("#"):
                    fields= line.split("\t")
                    if len(fields) == 3 :
                        evidencecode=fields[0].strip()
                        goref=fields[1].strip()
                        ecocode=fields[2].strip()
                        Eco.GAF_ECO_MAP[evidencecode+"-"+goref]=ecocode
                        Eco.ECO_EVIDENCE_MAP[ecocode]=evidencecode
                        eco_key=ecocode+"-"+goref
                        if eco_key not in ecoref_ev_map: ecoref_ev_map[eco_key]={}
                        ecoref_ev_map[eco_key][evidencecode]=1
                        #update evidence_codes list
                        if evidencecode not in evidence_codes: 
                            evidence_codes.append(evidencecode)
            
        #
        #Index ECO code using mgi-eco mapping file
        #
        
        if os.path.isfile('%s' % (local_mgi_eco_file)):
            reader=open(local_mgi_eco_file)
            mgi_eco_map={}
            eco_evidence_map={} #test the assumption that 1 eco code is assigned to only 1 evidence code
            eco_p= re.compile('ECO:\d+')
            for line in reader:
                fields=line.split("\t")
                if len(fields) == 3 :
                    evidencecode=fields[0].strip()
                    goref=fields[1].strip()
                    temp_line=line
                    #filter ecocode ()
                    ecolist=eco_p.findall(temp_line)
                    ecokey=evidencecode+"-"+goref
                    for ecocode in ecolist:
                        if not ecokey in mgi_eco_map:
                            mgi_eco_map[ecokey]={}
                            mgi_eco_map[ecokey][ecocode]=1
                        else: mgi_eco_map[ecokey][ecocode]=1
                        eco_key=ecocode+"-"+goref
                        if eco_key not in ecoref_ev_map: ecoref_ev_map[eco_key]={}
                        ecoref_ev_map[eco_key][evidencecode]=1
                        if ecocode not in eco_evidence_map:eco_evidence_map[ecocode]={}
                        eco_evidence_map[ecocode][evidencecode]=1
        # now generate Eco.MGI_ECO_MAP from mgi_eco_map
        for ecokey in mgi_eco_map:
            Eco.MGI_ECO_MAP[ecokey]="|".join(mgi_eco_map[ecokey].keys())
        for ecocode in eco_evidence_map:
            Eco.ECO_EVIDENCE_MAP[ecocode]="|".join(eco_evidence_map[ecocode].keys()) 
        for ecocode in ecoref_ev_map:
            Eco.ECO_GOREF2EVIDENCE_MAP[ecocode]="|".join(ecoref_ev_map[ecocode].keys())
        log.write("\n====================\nEco code - evidence code mapping\n")
        for ecocode in Eco.ECO_EVIDENCE_MAP:
            log.write(ecocode+"\t"+Eco.ECO_EVIDENCE_MAP[ecocode]+"\n")
        log.write("\n====================\nEco code - Goref to evidence code mapping\n")
        for ecocode_key in Eco.ECO_GOREF2EVIDENCE_MAP:
            log.write(ecocode_key+"\t"+Eco.ECO_GOREF2EVIDENCE_MAP[ecocode_key]+"\n")

        #print evidence code + GOF_REF => ECO code
        log.write("\n====================\n")
        log.write("Gaf Evidence code - GOF_REF => ECO code mapping\n")
        for key in sorted(Eco.GAF_ECO_MAP.iterkeys()):
            log.write(key+"\t"+Eco.GAF_ECO_MAP[key]+"\n")
        log.write("\n====================\n")
        log.write("Gaf Evidence code - MGI_REF => ECO code mapping\n")
        for key in sorted(Eco.MGI_ECO_MAP.iterkeys()):
            log.write(key+"\t"+Eco.MGI_ECO_MAP[key]+"\n") 
        log.write("\n====================\n")
        log.write("Current Gaf Evidence codes List\n")
        for evidencecode in evidence_codes: log.write(evidencecode+"\n")
          

    #
    # Returns the Evidence code associated with a given ECO code
    #
    # 
    def getEvidence_code(self,eco_code):
        eco_code=eco_code.strip()
        evidence_code=""
        if eco_code in Eco.ECO_EVIDENCE_MAP: 
            evidence_code=Eco.ECO_EVIDENCE_MAP[eco_code]
        return evidence_code

    #
    # Returns the Eco code in the gaf row given the evidence code and a GO_Ref
    # using the GAF_ECO_MAP: eco code=GAF_ECO_MAP[evidence_code+"-"+goref]
    # using the GAF_ECO_MAP: eco code=GAF_ECO_MAP[evidence_code+"-"+goref]
    #
    def getECO_code(self,goref,evidence_code,goref_collection={},evidence_map={}):
        evidence_code=evidence_code.strip()
        #check for invalid evidence codes
        ##if evidence_code not in evidence_codes: return ""
        goref_ids=goref.split("|")
        mgi_ref=""
        j_ref=""
        goref=""
        pubmed_ref=""
        eco_code=""
        for token in goref_ids:
            token=token.strip()
            if "GO_REF:" in token: goref=token
            elif "MGI:" in token:
                fields=token.split(":");
                if len(fields)<=2:mgi_ref=token
                else: mgi_ref=":".join(fields[1:2])
            elif "PMID:" in token:pubmed_ref=token
            elif "J:" in token:j_ref=token
        if not goref:
            if mgi_ref: goref=mgi_ref
            elif pubmed_ref: goref=pubmed_ref
            elif j_ref: goref=j_ref
            else: goref=goref_ids[0]
        #check if this is an external id
        if goref in goref_collection:
            goref= goref_collection[goref]

        goref=goref.strip()
        key=evidence_code+"-"+goref
        #if key in Eco.GAF_ECO_MAP: eco_code=Eco.GAF_ECO_MAP[key]
        if key in evidence_map: eco_code=evidence_map[key]
        else:
            key=evidence_code+"-"+"Default"
            eco_code=evidence_map[key]

        return eco_code

  
#class goref:

class gaf:

 #####################
 # Global variables ##
 #####################
 fields=[]
 fields_index={}
 taxon_index=-1
 aspect_index=-1
 goref_index=-1
 feature_type_index=-1
 evidence_index=-1
 db_index=-1
 object_index=-1
 qual_index=-1
 product_form_id_index=-1
 object_symbol_index=-1
 go_id_index=-1
 with_index=-1
 object_name_index=-1
 object_synonym_index=-1
 date_index=-1
 assigned_by_index=-1
 annotation_extension_index=-1

 ##########################
 # class constructor    ###
 # default version: 2.0 ###
 ##########################
 def _init(self,gaf_version="2.0"):
     
     if "1." in gaf_version:
        for i in range(len(GAF1_FIELDS)):
           self.fields_index[GAF1_FIELDS[i]]=i
           self.fields.append(GAF1_FIELDS[i])
     else:
        for i in range(len(GAF_FIELDS)):
           self.fields_index[GAF_FIELDS[i]]=i
           self.fields.append(GAF_FIELDS[i])

     self.taxon_index= gaf.fields.index("Taxon")
     self.aspect_index= gaf.fields.index("Aspect")
     self.goref_index=gaf.fields.index("DB_Reference")
     self.feature_type_index=gaf.fields.index("DB_Object_Type")
     self.evidence_index=gaf.fields.index("Evidence_Code")
     self.db_index=gaf.fields.index("DB")
     self.object_index=gaf.fields.index("DB_Object_ID")
     self.qual_index = gaf.fields.index("Qualifier")
     if "Gene_Product_Form_ID" in self.fields:
         self.product_form_id_index=self.fields.index("Gene_Product_Form_ID")
     self.object_symbol_index= self.fields.index('DB_Object_Symbol')
     self.go_id_index= self.fields.index('GO_ID')
     self.with_index=self.fields.index('With')
     self.object_name_index= self.fields.index('DB_Object_Name')
     self.object_synonym_index=self.fields.index('DB_Object_Synonym')
     self.date_index= self.fields.index('Date')
     self.assigned_by_index=self.fields.index('Assigned_By')
     if 'Annotation_Extension' in self.fields:
         self.annotation_extension_index=self.fields.index('Annotation_Extension')

 #############################################################
 ##
 # Is_geneVariant returns true if a gaf annotation is
 # a form of a gene (protein,...). Returns false otherwise
 #############################################################
 def is_gene_variant(self,gaf_row):
    is_g_variant=0
    if self.product_form_id_index in range(len(gaf_row)):
       if ":" in gaf_row[self.product_form_id_index]:
           is_g_variant=1
    return is_g_variant

 #
 #Return true/false whether or not a given gaf_row has 
 # missing mandatory fields (fields are empty)
 # False=0; True = index of the field + 1
 #
 def has_missing_fields(self,gaf_row):
    missing= 0
    for i in range(len(gaf_row)):
        field=self.fields[i]
        if field in REQUIRED_GAF_FIELDS:
           if not gaf_row[i]:
              missing=i+1
              break

    return missing
 #
 ##Return true/false whether or not a given gpi_row has 
 # missing mandatory fields (fields are empty)
 # False=0; True = index of the field + 1
 #
 def gpi_has_missing_fields(self,gpi_row):
    missing= 0
    for i in range(len(gpi_row)):
        field=GPI_FIELDS[i]
        if field in REQUIRED_GPI_FIELDS:
           if not gpi_row[i]:
              missing=i+1
              break

    return missing
 #
 ##Return true/false whether or not a given gpad_row has 
 # missing mandatory fields (fields are empty)
 # False=0; True = index of the field + 1
 #
 def gpad_has_missing_fields(self,gpad_row):
    missing= 0
    for i in range(len(gpad_row)):
        field=GPAD_FIELDS[i]
        if field in REQUIRED_GPAD_FIELDS:
           if not gpad_row[i]:
              missing=i+1
              break

    return missing
 #
 # Get GPI parent_id field from GAF.db and GAF.object
 # olny if the GAF.product_form_id field is not empty
 #
 def getParent_gp_id(self,gaf_row):
    return gaf_row[self.db_index]+":"+gaf_row[self.object_index]
 
 #
 #Compute the GPAD relationship field 
 # using gaf qualifier and the aspect fields 
 #
 def getGPAD_relationship(self,gaf_row):
    qualifiers= gaf_row[self.qual_index].lower().split("|")
    relations=[]
    for qualifier in qualifier_order:
        if qualifier in qualifiers:
           relations.append(qualifier)
    aspect="default"
    if gaf_row[self.aspect_index]:
       aspect=gaf_row[self.aspect_index]
    if aspect in aspect2rln:
        relations.append(aspect2rln[aspect])
        return "|".join(relations)
    else:
        return ""

 #
 #   If protein object: (only for gaf 2.0 version)
 #   set protein_map[protein]=object
 #   set the DB of this object to the protein db
 #   set the DB_Oject_ID of this object to the protein ID
 #   
 def setProtein(self,protein_map,gaf_row):
     total=0;
     if self.product_form_id_index not in range(len(gaf_row)):
        return total
     protein=gaf_row[self.product_form_id_index] #GAF - field 17
     if not protein:
        total=1
     else: #get the DB and object ID from field 17
        fields=protein.split(":")
        if len(fields)>=2:
           protein_map[protein]=gaf_row[self.object_index]      #DB_Object_ID - field #2
           gaf_row[self.db_index]=fields[0]                     #set the DB of this object to the protein db
           last=len(fields)
           pid=""
           if len(fields)==2:pid=fields[1]
           else:pid=":".join(fields[1:last])
           gaf_row[self.object_index]=pid                       #set the DB_Oject_ID of this object to the protein ID
     return total

 #
 #Returns columns mapping between GAF and GPAD/GPI
 #GAF_column#\tcolumn_label\tgpad_column#\tgpi_column#
 #
 def getgaf2gpadgpi(self):
    row=""
    for i in range(len(self.fields)):
        #get index of this field in gpad and gpi files
        gpi_index=-1
        gpad_index=-1
        if self.fields[i] in GPI_FIELDS: 
            gpi_index=GPI_FIELDS.index(self.fields[i])
        if self.fields[i] in GPAD_FIELDS:
            gpad_index=GPAD_FIELDS.index(self.fields[i])
        if i==self.product_form_id_index:
            gpi_index=gpi_parent_index
        gpad_label="N/A"
        gpi_label="N/A"
        if(gpad_index>=0):
            gpad_label="%d"%(gpad_index+1)
        if(gpi_index>=0):
            if (gpi_index==gpi_db_index) or (gpi_index==gpi_db_object_index):
                gpi_label="%d/%d"%(gpi_index+1,gpi_parent_index+1)
            else:
                gpi_label="%d"%(gpi_index+1)
        row+="\n!\t(%d) %s\t%s\t%s" %(i+1,GAF_FIELDS_LABEL[gaf.fields[i]],gpad_label,gpi_label)
    return row

 #
 #Returns columns mapping between GPI and GAF
 #
 def getgpi2gaf(self):
    rows=""
    for i in range(len(GPI_FIELDS)):
        if GPI_FIELDS[i] in self.fields:
           gaf_index=self.fields.index(GPI_FIELDS[i])
           if i == gpi_db_index :
              rows+="\n!\t(%d) %s\t%d/%d" %(i+1,GPI_FIELDS_LABEL[GPI_FIELDS[i]],gaf_index+1,len(self.fields))
           elif i == gpi_db_object_index:
              rows+="\n!\t(%d) %s\t%d/%d" %(i+1,GPI_FIELDS_LABEL[GPI_FIELDS[i]],gaf_index+1,len(self.fields))
           else: rows+="\n!\t(%d) %s\t%d" %(i+1,GPI_FIELDS_LABEL[GPI_FIELDS[i]],gaf_index+1)
        else:
           if i == gpi_parent_index:
              rows+="\n!\t(%d) %s\t%d + %d" %(i+1,GPI_FIELDS_LABEL[GPI_FIELDS[i]],self.db_index+1,self.object_index+1)
           else:
              rows+="\n!\t(%d)%s\tn/a" % (i+1,GPI_FIELDS_LABEL[GPI_FIELDS[i]])
    return rows

 #
 #Returns columns mapping between GPAD and GAF
 #
 def getgpad2gaf(self):
    rows=""
    for i in range(len(GPAD_FIELDS)):
           if GPAD_FIELDS[i] in self.fields:
              gaf_index=self.fields.index(GPAD_FIELDS[i])
              if i == gpad_evidence_index :
                 rows+="\n!\t(%d) %s\t%d + %d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],gaf_index+1,gaf_index)
              elif i == gpad_db_object_index :
                 rows+="\n!\t(%d) %s\t%d/%d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],gaf_index+1,len(self.fields))
              elif i == gpad_db_index:
                 rows+="\n!\t(%d) %s\t%d/%d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],gaf_index+1,len(self.fields))
              elif i==gpad_relationship_index:
                 rows+="\n!\t(%d) %s\t%d|%d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],self.qual_index+1,self.aspect_index+1)
              else: rows+="\n!\t(%d) %s\t%d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],gaf_index+1)
           else: 
              if i == gpad_annot_extension_index:
                  rows+="\n!\t(%d) %s\t%d" %(i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]],self.evidence_index+1)
              else:
                 rows+="\n!\t(%d) %s\tn/a" % (i+1,GPAD_FIELDS_LABEL[GPAD_FIELDS[i]])
    return rows

