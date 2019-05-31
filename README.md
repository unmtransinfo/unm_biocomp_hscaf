# `UNM_BIOCOMP_HSCAF` <img align="right" src="/doc/images/hscaf_logo.png">

Implements the HierS hierarchical scaffold algorithm by Wilkens et al.
Designed for HierS scaffold analysis, for single molecules, datasets, or
common-scaffold comparisons between pairs or groups of molecules.

## Dependencies
* Java 8
* Maven 3.5+
* `unm_biocomp_depict`, `unm_biocomp_db`, `unm_biocomp_util`
* ChemAxon JChem (19.3.0 ok)
* Access to [ChemAxon Maven repository](https://hub.chemaxon.com) (see [documentation](https://docs.chemaxon.com/display/docs/Public+Repository)).
  * Requires ChemAxon-Hub Artifactory credentials.

## Compiling

```
mvn clean install
```

## Usage

```
mvn exec:java -Dexec.mainClass="edu.unm.health.biocomp.hscaf.hier_scaffolds"
```

Developed at the [UNM Translational Informatics Division](http://datascience.unm.edu).

* [Public web app](http://pasilla.health.unm.edu/tomcat/hscaf)

