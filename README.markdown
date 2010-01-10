
This is the source code implementing the Gemoda algorithm described in the following paper

[A generic motif discovery algorithm for sequential data.](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/bti745v1)  
[Jensen KL](http://www.epernicus.com/klj), Styczynski MP, Rigoutsos I, Stephanopoulos GN.  
[Bioinformatics](http://bioinformatics.oxfordjournals.org/). 2006 Jan 1;22(1):21-8. Epub 2005 Oct 27.

> *MOTIVATION:* Motif discovery in sequential data is a problem of great interest and with many applications. However, previous methods have been unable to combine exhaustive search with complex motif representations and are each typically only applicable to a certain class of problems.
>  
> *RESULTS:* Here we present a generic motif discovery algorithm (Gemoda) for sequential data. Gemoda can be applied to any dataset with a sequential character, including both categorical and real-valued data. As we show, Gemoda deterministically discovers motifs that are maximal in composition and length. As well, the algorithm allows any choice of similarity metric for finding motifs. Finally, Gemoda's output motifs are representation-agnostic: they can be represented using regular expressions, position weight matrices or any number of other models for any type of sequential data. We demonstrate a number of applications of the algorithm, including the discovery of motifs in amino acids sequences, a new solution to the (l,d)-motif problem in DNA sequences and the discovery of conserved protein substructures

Use the following sequence of commands to install the Gemoda implementation provided by this software:

./configure  
make   
make install

Note that the software requires the [GNU Scientific Library](http://www.gnu.org/software/gsl/).

To familiarize you with the software,  we included 3 examples in the examples/ subdirectory.  The README files in these folders details how to run the software.


Kyle  
<kljensen@gmail.com>