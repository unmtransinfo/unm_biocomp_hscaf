# Welcome to UNM-Biocomp-HScaf #

## Introduction ##

Java package edu.unm.health.biocomp.hscaf, from the [University of New Mexico Translational Informatics Division](http://medicine.unm.edu/informatics/), is an open-source library built upon the JChem toolkit from [ChemAxon](http://www.chemaxon.com), for chemical scaffold analysis, for single molecules, datasets, or common-scaffold comparisons of pairs or groups of molecules.

## Features ##

  * The HierS hierarchical scaffolds algorithm published by Wilkens et al. (2005) is implemented.
  * This library is dependent on cheminformatics toolkit JChem from [ChemAxon Ltd.](http://www.chemaxon.com)
  * Individual molecules or large datasets can be analyzed.
  * Scaffold trees define a comprehensive, hierarchical set of scaffolds for a molecule, the largest of which corresponds to the Bemis-Murcko framework.
  * Common scaffolds between two molecules can be identified.
  * Analyze a molecule and identify all scaffolds, linkers and side-chains.
  * Scaffold attachment points can be identified.
  * Unique sets of scaffolds for a molecule or for a dataset of many molecules can be found.
  * The hierarchical scaffold tree for a molecule can be generated.
  * An N-attached variant of the HierS scaffold definition may be used, where atoms single-bonded to ring nitrogens are considered part of the scaffold.
  * Optionally scaffolds may be stereo or (by default) non-stereo.
  * Maximum common scaffold of two molecules, maximum common subtrees of two scaffold trees.

## Used By ##

  * The <a href='http://carlsbad.health.unm.edu'>CARLSBAD Project</a> is an NIH-funded bioactivity database project utilizing this library for scaffold analysis.
  * The public web app <a href='http://pasilla.health.unm.edu/tomcat/biocomp/badapple'>BADAPPLE</a> utilizes this library for scaffold analysis.  BADAPPLE = Bio Activity Data Associative Promiscuity Pattern Learning Engine.
  * The <a href='http://bard.nih.gov'>BARD Project</a> (BioAssay Research Database) includes a BADAPPLE plugin component which uses unm-biocomp-hscaf.

## Author ##

  * Jeremy Yang

## Prerequisites ##

  * Java 1.6 (SE6)
  * ChemAxon JChem 5.5.0+ [5.8.3 recommended - Oct. 2012]
  * BerkeleyDB (optional, faster for larger datasets)
  * PostgreSQL (optional, faster for largest datasets)

## Documentation ##

  * [API (Javadoc)](http://unm-biocomp-hscaf.googlecode.com/svn/trunk/doc/hscaf/index.html)
  * [Getting Started](GettingStarted.md)

## Support ##

This software is provided as-is, with no warranty and no support.

## Versions, Release Notes ##

| **date** | **description** |
|:---------|:----------------|
|Jun 30, 2011 | Initial release |
|Oct 27, 2011 | Several improvements including faster due to use of BerkeleyDB |
|Feb 26, 2012 | Several improvements.  Faster, partly by simply avoiding pathological mega-scaffolds. |
|Oct 2, 2012 | Much improved again.  Option for PostgreSQL back end has many advantages.  |
|May 31, 2013 | Mostly bugfixes.  Webapp included.    |

## License ##

This code released under terms of the [BSD 3-Clause License](http://www.opensource.org/licenses/BSD-3-Clause).

## Acknowledgements ##

  * This software project was partially funded via NIH award 1R21GM095952-01, "Chemical Pattern Detection and Visualization in Biological Networks", PI Tudor Oprea.
  * Includes source code for class [ShortestPath](http://unm-biocomp-hscaf.googlecode.com/svn/trunk/doc/edu/unm/health/biocomp/hscaf/ShortestPath.html), written by Miklos Vargyas, ChemAxon, and obtained from the ChemAxon forum.
  * Thanks to the fine work and 2005 publication of the HierS algorithm by Steven J. Wilkens, Jeff Janes and Andrew I. Su, of the Genomics Institute of the Novartis Research Foundation.
  * Thanks to colleagues at UNM Biocomputing for supporting this project, especially Cristian Bologa for conceiving the project.

## References ##

  1. HierS: Hierarchical Scaffold Clustering, Steven J. Wilkens, Jeff Janes, Andrew I. Su, J. Med. Chem. 2005, 48, 3182-3193.
  1. [ChemAxon Ltd.](http://www.chemaxon.com), Budapest, Hungary.