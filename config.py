#!/usr/bin/env python
#
# The main config file for the file processor
#
import os
#
#This file contain global variables
GOA_DATA_BASE="ftp://ftp.ebi.ac.uk/pub/databases/GO/goa"
GOA_GPAD_SPEC_BASE="ftp://ftp.geneontology.org/go/www/GO.format.gpad.shtml"
GOA_GPI_SPEC_BASE="ftp://ftp.geneontology.org/go/www/GO.format.gpi.shtml"
#http://www.geneontology.org/software/reporting/GOBO/AnnotationFormats.pm
#
#
GAF_FIELDS_COUNT=17
GPAD_FIELDS_COUNT=12
GPI_FIELDS_COUNT=9
#
#local index file update setup - in days
#
DAILY_UPDATE=1
WEEKLY_UPDATE=7
BIWEEKLY_UPDATE=14
UPDATE_INDEX={}
#
# Converter base
#
GAF_CONVERTER_BASE=os.path.abspath(os.path.dirname(__file__))
GAF_CONVERTER_LOG_BASE=GAF_CONVERTER_BASE+"/log"
GAF_CONVERTER_DATA_BASE=GAF_CONVERTER_BASE+"/data"
GAF_CONVERTER_DATA_TEMP=GAF_CONVERTER_BASE+"/data/temp"
#
#For converter files header
#
GAF_VERSION="2.0"
GAF1_VERSION="1.0"
GPAD_VERSION="1.0"
GPI_VERSION="1.0"
CVS_VERSION="1.0"
DB_URL="http://www.informatics.jax.org/"
PROJECT_NAME="Mouse Genome Informatic (MGI)"
CONTACT_EMAIL="mgi-help@jax.org"

#
#Converter  GAF - ECO evidence map - dependency 
#
GAF_ECO_MAP="gaf-eco-mapping.txt"
GAF_ECO_MAP_URL="http://purl.obolibrary.org/obo/eco/"
REMOTE_MAP_FILE=GAF_ECO_MAP_URL+GAF_ECO_MAP
LOCAL_MAP_FILE=GAF_CONVERTER_DATA_BASE+"/"+GAF_ECO_MAP
UPDATE_INDEX[LOCAL_MAP_FILE]=BIWEEKLY_UPDATE
#
#
#
ECO_EVIDENCE="eco.obo"
ECO_EVIDENCE_URL="https://evidenceontology.googlecode.com/svn/trunk/"
REMOTE_ECO_EVIDENCE_FILE=ECO_EVIDENCE_URL+ECO_EVIDENCE
LOCAL_ECO_EVIDENCE_FILE=GAF_CONVERTER_DATA_BASE+"/"+ECO_EVIDENCE
UPDATE_INDEX[LOCAL_ECO_EVIDENCE_FILE]=BIWEEKLY_UPDATE
#
#GO Reference Collection
#
GO_REF_COLLECTION_URL="http://www.geneontology.org/cgi-bin/references.cgi"
GO_REF_COLLECTION="GO.references"
GO_REF_COLLECTION_BASE="http://www.geneontology.org/doc/"
REMOTE_GO_REF_FILE=GO_REF_COLLECTION_BASE+GO_REF_COLLECTION
LOCAL_GO_REF_FILE=GAF_CONVERTER_DATA_BASE+"/"+GO_REF_COLLECTION
UPDATE_INDEX[LOCAL_GO_REF_FILE]=BIWEEKLY_UPDATE
#
# GO database cross-reference
#
GO_XREF_DB="GO.xrf_abbs"
REMOTE_GO_XREF_FILE=GO_REF_COLLECTION_BASE+GO_XREF_DB
LOCAL_GO_XREF_FILE=GAF_CONVERTER_DATA_BASE+"/"+GO_XREF_DB
UPDATE_INDEX[LOCAL_GO_XREF_FILE]=DAILY_UPDATE

#Converter MGI report dependency
# MGI marker reports path
# This report is used to generate mouse gpi files
#
MGI_FTP_HOST="ftp://ftp.informatics.jax.org"
MGI_REPORT_PATH="/pub/reports/"
MGI_MRK_REPORT="MRK_List2.rpt"
REMOTE_MRK_REPORT_FILE=MGI_FTP_HOST+MGI_REPORT_PATH+MGI_MRK_REPORT
LOCAL_MRK_REPORT_FILE=GAF_CONVERTER_DATA_BASE+"/"+MGI_MRK_REPORT
UPDATE_INDEX[LOCAL_MRK_REPORT_FILE]=WEEKLY_UPDATE
#
# GO Ontology 
#
GO_TERM_REPORT="go_terms.mgi"
REMOTE_GO_REPORT_FILE=MGI_FTP_HOST+MGI_REPORT_PATH+GO_TERM_REPORT
LOCAL_GO_REPORT_FILE=GAF_CONVERTER_DATA_BASE+"/"+GO_TERM_REPORT
UPDATE_INDEX[LOCAL_GO_REPORT_FILE]=WEEKLY_UPDATE
#
#MGI GO_eco_association report
#
MGI_GO_ECO_REPORT="GO_eco_association.rpt"
REMOTE_MGI_GO_REPORT_FILE=MGI_FTP_HOST+MGI_REPORT_PATH+MGI_GO_ECO_REPORT
LOCAL_MGI_GO_ECO_REPORT_FILE=GAF_CONVERTER_DATA_BASE+"/"+MGI_GO_ECO_REPORT
UPDATE_INDEX[LOCAL_MGI_GO_ECO_REPORT_FILE]=DAILY_UPDATE

