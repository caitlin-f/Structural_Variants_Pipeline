# Structural_Variants_Pipeline
Updated 22/10/2017

Installation instructions.pdf
	- Installation instructions for MacOS for the tools used in largeSV_pipeline.sh 

LargeSV_pipeline.sh
	- Executable script for finding large structural variations and outputting a distance tree (computed in R)

collate_output.py
	- called within largeSV_pipeline.sh to collate output from the various tools into separate <samplename>.txt files as well as a combined all_data.txt file

dellytree.py
	- called within largeSV_pipeline.sh to output a matrix (distmat.csv) that is used by largeSVtree.R

largeSVtree.R
	- called within largeSV_pipeline.sh to produce a jpeg of the final tree calculated using distmat.csv
