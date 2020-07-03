# 3StronglyRegularGraphs

This repository collects some research, code, and data I created while investigating a graph theoretic problem during the spring of 2020. Building off the concept of strongly regular graphs, which have been extensively researched, I extended the concept to a more general concept of three-strongly regular graphs, or 3SRGs. 

## Contents

The extent of my research is documented in ThreeStronglyRegularGraphs.pdf; this is an unpublished and unfinished piece of research, and as such has gaps in citations and detail. An unofficial list of useful sources, as a base of a bibliography, is included in biblio.txt. In graphTheoryFunctions.py, I've collected the python code I've used to find, categorize and examine 3SRGs, as well as useful helper functions for working with graph data in text or exporting to LaTeX. In importantGraphs.txt I've collected all examples of nontrivial 3SRGs that I've found so far. 

## Prerequisites

Prerequisites of the work are an understanding of basic graph theory, some algebraic graph theory, and some of the theory of combinatorial designs. I've mentioned the key textbooks I used to learn these topics in biblio.txt.

## Acknowledgements

Thanks to Professor Emily King of CSU for helping suggest this problem and teaching me how mathematical research works. I used data compiled by Brendan McKay of Australian National University to search for 3SRGs-his website is http://users.cecs.anu.edu.au/~bdm/. When dealing with representations of graphs in plaintext, I used the g6 format-full explanation of this format can be found at https://users.cecs.anu.edu.au/~bdm/data/formats.html. Helper functions for converting from g6 to numpy adjacency matrices and back are included in graphTheoryFunctions.py.


