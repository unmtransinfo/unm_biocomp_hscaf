package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.SearchException;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents any single HierS sidechain
	<br />
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@author Jeremy J Yang
*/
public class Sidechain extends Molecule
{
  /** JChem format for canonical sidechain smiles */
  private static final String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static final String cansmifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w-a";
  private static final String cxsmifmt="cxsmiles:+L+l+e+d+D+p+R+f+w";
  /** used for equality comparison */
  private String cansmi;
  private String smi;
  /** unique ID for dataset scope */
  private Integer id;
  /////////////////////////////////////////////////////////////////////////////
  private Sidechain() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	When initialized from a molecule, sidechains are removed
	and junction pseudo atoms remain at each junction bond&#46;
  */
  public Sidechain(Molecule mol)
    throws MolFormatException,MolExportException,SearchException
  {
    mol.clonecopy(this);
    hscaf_utils.replaceJPseudoatomsWithJHydrogens(this);
    this.cansmi=this.exportToFormat(cansmifmt);
    this.smi=this.exportToFormat(smifmt);
    this.id=0;
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(Sidechain chain2)
  {
    return (this.cansmi.equals(chain2.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get unique ID for dataset scope.
  */
  public int getID() { return this.id; }
  /**	Set unique ID for dataset scope.
  */
  public void setID(int id2) { this.id=id2; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get canonical SMILES for this sidechain (as used for equality comparison).
  */
  public String getCansmi() { return this.cansmi; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this sidechain.
  */
  public String getSmi() { return this.smi; }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this sidechain with junction atoms as pseudo-atoms&#46;
	Note this is in ChemAxon cxsmiles format e&#46;g&#46;
	"*CCC |$p_J;;;|"&#46;
  */
  public String getJsmi()
    throws MolExportException
  {
    Molecule chainmol=this.cloneMolecule();
    hscaf_utils.replaceJHydrogensWithJPseudoatoms(chainmol);
    return chainmol.exportToFormat(cxsmifmt);
  }
  ///////////////////////////////////////////////////////////////////////////
}
