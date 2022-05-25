# `UNM_BIOCOMP_HSCAF` <img align="right" src="/doc/images/hscaf_logo.png">

Implements the HierS hierarchical scaffold algorithm (Wilkens, 2005).
Designed for HierS scaffold analysis, for single molecules, datasets, or
common-scaffold comparisons between pairs or groups of molecules.

## Dependencies

* Java 8
* Maven 3.5+
* [`unm_biocomp_depict`](https://github.com/unmtransinfo/unm_biocomp_depict)
* [`unm_biocomp_util`](https://github.com/unmtransinfo/unm_biocomp_util)
* ChemAxon JChem (19.3.0 ok)
* Access to [ChemAxon Maven repository](https://hub.chemaxon.com) (see [documentation](https://docs.chemaxon.com/display/docs/Public+Repository)).
  * Requires ChemAxon-Hub Artifactory credentials.

## Compiling

```
mvn clean install
```

## Deploying `HSCAF`

```
mvn --projects hscaf_war tomcat7:deploy
```

or

```
mvn --projects hscaf_war tomcat7:redeploy
```

## Usage

```
mvn --projects unm_biocomp_hscaf exec:java -Dexec.mainClass="edu.unm.health.biocomp.hscaf.hier_scaffolds" -Dexec.args="-i hscaf_testset.smi -inc_mol -inc_scaf -o hscaf_testset_out.smi -out_scaf hscaf_testset_scaf.smi -vvv"
```

```
mvn --projects unm_biocomp_hscaf exec:java -Dexec.mainClass="edu.unm.health.biocomp.hscaf.bm_scaffolds" -Dexec.args="-i hscaf_testset.smi"
```

Developed at the [UNM Translational Informatics Division](https://datascience.unm.edu).

* [Public web app](https://datascience.unm.edu/tomcat/hscaf)
* HierS: hierarchical scaffold clustering using topological chemical graphs. Wilkens SJ(1), Janes J, Su AI. J Med Chem. 2005 May 5;48(9):3182-93.
