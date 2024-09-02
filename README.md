####################
knotty.py 
####################

DESCRIPTION:
	Python tool for knot detection in protein structures using implementation of the Taylor Knot Detection Algorithm and protein structure visualization.
	
USAGE:
	knotty.py [FUNCTION] [OPTIONS]
	
FUNCTIONS: find, visualize

########################################
  
   FIND:
	Python Implementation of the Taylor Knot Detection Algorithm
 	
	OPTIONS:
	-i [ --input ] 		file of amino acid sequence coordinates (.crd file)

	-o [ --output ] 	Name of Output File

	-m [ --max-iterations ]         maximum number of iterations for smoothing algorithm (default: 250)
	
	-e [ --epsilon ]         threshold (in Angstroms) for approximating collinearity of alpha carbons (default: 0.25)


   VISUALIZE:
        3D interative visualization tool for amino acid chain structure. 

	OPTIONS:
        -i [ --input ]          file of amino acid sequence coordinates (.crd file)


########################################

Last Updated: March 13th, 2020
