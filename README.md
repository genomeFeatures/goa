Create standard python libraries for GO Annottation. 
These libraries can be used to parse GO files and convert one file type to another. 

Example: One could use the library- to convert GPAD/GPI files to GAF file and vise versa -

To generate the mapping between Evidence code <-> ECO code, the converter indexes
publicly available maps - for both the GO reference and the GAF-ECO mapping.
The files are fetched from the following sources and indexed locally for further use:
  1) GO Reference Collection - dependency :http://www.geneontology.org/doc/GO.references
  2) GAF - ECO evidence map - dependency :http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt
  3) GAF - ECO evidence map - dependency: ftp://ftp.informatics.jax.org/pub/reports/GO_eco_association.rpt
This allowed the converter to resolve MGI IDS, J#s, Pubmed Ids,... into GO_REF

To generate MGI specific GPI file,the converter indexes the following file:
  1) ftp://ftp.informatics.jax.org/pub/reports/MRK_List2.rpt

Each indexed file has an expiration time (update frequency) - for example 
MGI reports are generated daily so the associated indexed reports will be updated daily.
The update frequency of each file is defined in config.py configuration file.

Step 1: GAF2GPAD converter - gaf2gpad.py
  Usage: run gaf2gpad.py on the command line for usage

  The converter is written with the option --mgi that if present,
  the converter will generate the GPI file using additional
  annotions from MGI marker report - If this option is missing,
  the GPI file will be generated using only the information in the GAF

 Logs: This script generates three logs stored under the log/ directory relative
 to the converter base.
 a) gaf2gpad.py.log to keep track of the process logs, maps used
 b) gaf_file.log  to keep track of user input GAF file process QC and data tally
 c) maps.log to keep track of all the mapping generated and used

 data: The converter uses external files to create data maps
       these files are indexed locally under data/ directory relative 
   to the converter base.The following files are indexed locally 
 a) GO.references --> GO Reference Collection
 b) MRK_List2.rpt --> MGI marker reports
 c) gaf-eco-mapping.txt  --> GAF - ECO evidence map
 d) GO_eco_association.rpt  --> MGI GAF - ECO evidence map
 e) go_terms.mgi --> GO Ontology

Step 2: GPAD2GAF converter - gpad2gaf.py
  Usage: run gpad2gaf.py for usage

 

