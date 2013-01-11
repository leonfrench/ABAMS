Below are the steps to run the BAMS and Allen code created by Leon French. The performs various operations on data from the Allen Brain Atlas and Brain Architecture management system. Further information can be found at "http://www.chibi.ubc.ca/ABAMS":http://www.chibi.ubc.ca/ABAMS .

Tested on Ubuntu (desktop) and CentOS (server). Installation steps tested on pavnote-05 (ubuntu).

h3. Setup steps

* download ABAMSSupplements.zip file (151mb, 560mb uncompressed):  "https://docs.google.com/file/d/0B3w9lE7AjmJTb292NGpHYk12SGc/edit":https://docs.google.com/file/d/0B3w9lE7AjmJTb292NGpHYk12SGc/edit
* install maven
* install eclipse
** create M2_REPO variable to point to maven repository (usually ~/.m2)
** install egit (optional)
* download source (from github)
* run mvn compile in project home folder
** use Maven to download dependencies and eclipse for compilation
* edit ABAMS.properties (in source home folder)
** change/replace all occurances of the "/home/abams/ABAMSData/" folder to the location you unzipped the supplements
** set number of threads to use (whitetext.max_threads, used for mantel test oprimization)
* install ermineJ (optional)
** set ermineJ location in ABAMS.properties (abams.ermineJ.bin)



h3. Example run
* try to run ubic.BAMSandAllen.MatrixPairs/ConnectivityAndAllenExpressionMatrixPair.java in eclipse.
** it will recreate the global mantel test results from ""Relationships between gene expression and brain wiring in the adult rodent brain.":http://www.ncbi.nlm.nih.gov/pubmed/21253556?dopt=Abstract"
* also you can try ubic.BAMSandAllen.optimize.GreedyMultiThreaded.java to optimize by removing genes, don't forget to set the options like number of threads in the start of the main function.

h3. Contact

If you have any trouble running it, please contact me at leonfrench@gmail.com