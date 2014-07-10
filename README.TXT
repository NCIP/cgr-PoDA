Because there are many pathways to test and many permutations to run, 
the PoDA analysis is best run as a swarm of jobs on biowulf.  This
is designed to help you do that.  The scripts will break the pathways
up into 36 block of 40 pathways each, and will break the permutations
into 100 blocks of 10 runs each (1000 permutations, total).

To run the NASH analysis, you will need 4 files:

* nash-dataForPoDA.RData -- the NASH data
* PoDA.R -- the R PoDA script
* nashPoDA.skel.R -- a "skeleton" file for each chunk
* nashPoDA.swarmgen.sh -- a shell script to generate the swarm 
* getResults.R -- an R script to concatenate the swarm results

You will need to EDIT two of these files (nashPoDA.skel.R
and nashPoDA.swarmgen.sh) before you can submit the job.


STEP 1: look at the data!

In R,
	load("nash-dataForPoDA.RData")
	# what is here?
	ls()
	# nashguide will tell you what all the variables are;
	# notice we have the pathway/gene/snp mappings 
  # as well as the nash genotype & clinical data
	nashguide
	# here's a glimpse of the clinical data (1st 3 rows):
  nashclin[1:3,]

Ultimately, you will want to run PoDA using different clinical 
outcomes, so it is good to be familiar with what is in
nashclin.  For now, let us suppose we want to use nashdx
as the outcome, take "2" as a case and "0" as a control
(tossing the others):

	# nashdx is 0, 1, 1.5, or 2:
  table(nashclin$nashdx)
	

STEP 2: make a copy!

Because the output file names are all the same, you'll want to make a
(meaningfully-named) subdirectory for the particular outcome you're
running PoDA on.  For example, let's say we're looking at nashdx; at
the command line, in this directory, do:

	mkdir nashdx_poda # or some appropriate name
	cd nashdx_poda
	# link the ones we don't need to edit
	ln -s ../PoDA.R ../getResults.R ../nash-dataForPoDA.RData .
	# copy the ones that are run-specific
	cp ../nashPoDA.skel.R ../nashPoDA.swarmgen.sh .
	# make note of the directory!
	pwd

STEP 3: edit nashPoDA.skel.R !

nashPoDA.skel.R is an R script.  We need to tell it who is a case
and who is a control by editing lines 8 & 9.  If you look at it
with
  head nashPoDA.skel.R
you'll see that right now, we're using nashdx==2 as a case,
nashdx==0 as a control, and nashdx==1 or nashdx==1.5 are NA'd.
You can set it to whatever you need.

STEP 4: edit nashPoDA.swarmgen.sh !

You need to change where it says "CHANGE.THIS" to be
the output of the pwd command, ie, the full path of the
NASH PoDA subdir you just made.

STEP 5: start the run!

At a shell prompt in the NASH PoDA subdir, 

	# First, generate the swarm files:
  ./nashPoDA.swarmgen.sh

	# Then submit the swarm:
  swarm -f nashPoDA.swarm


STEP 6: Wait for it to finish and collect the results!

First, get an interactive node on biowulf with 

  qsub -I -l nodes=1

Then go to the appropriate directory and start R: 

  cd NAME-0F-THE-ANALYSIS-DIRECTORY
  R

In R, run

  source("getResults.R")

to create a table of the results as well as a data
file that contains the results as a data frame and
the PoDA "S" values for each sample/pathway as 
a matrix.  These will be saved in
* nash.PoDA_table.txt (tab-delimitted table of results)
* nash.PoDA_table_and_S_vals.Rdata (table and S-values as R data)





