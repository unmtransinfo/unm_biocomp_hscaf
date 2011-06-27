package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.SearchException;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents any single HierS linker
	<br />
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@author Jeremy J Yang
*/
public class Linker extends Molecule
{
  /** JChem format for canonical linker smiles */
  private static final String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static final String cansmifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w-a";
  private static final String cxsmifmt="cxsmiles:+L+l+e+d+D+p+R+f+w";
  /** used for equality comparison */
  private String cansmi;
  private String smi;
  /** unique ID for dataset scope */
  private Integer id;
  /** number of junctions (2+) */
  private Integer degree;
  /////////////////////////////////////////////////////////////////////////////
  private Linker() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	When initialized from a molecule, sidechains are removed
	and junction pseudo atoms remain at each junction bond&#46;
  */
  public Linker(Molecule mol)
    throws MolFormatException,MolExportException,SearchException
  {
    mol.clonecopy(this);
    this.degree=pruneLinker(this);
    hscaf_utils.replaceJPseudoatomsWithJHydrogens(this);
    this.cansmi=this.exportToFormat(cansmifmt);
    this.smi=this.exportToFormat(smifmt);
    this.id=0;
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(Linker linker2)
  {
    return (this.cansmi.equals(linker2.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get unique ID for dataset scope.
  */
  public int getID() { return this.id; }
  /**	Set unique ID for dataset scope.
  */
  public void setID(int id2) { this.id=id2; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get degree, number of junctions.
  */
  public int getDegree() { return this.degree; }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get canonical SMILES for this linker (as used for equality comparison).
  */
  public String getCansmi() { return this.cansmi; }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this linker.
  */
  public String getSmi() { return this.smi; }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this linker with junction atoms as pseudo-atoms&#46;
	Note this is in ChemAxon cxsmiles format e&#46;g&#46;
	"*CCC* |$p_J;;;;p_J$|"&#46;
  */
  public String getJsmi()
    throws MolExportException
  {
    Molecule linkmol=this.cloneMolecule();
    hscaf_utils.replaceJHydrogensWithJPseudoatoms(linkmol);
    return linkmol.exportToFormat(cxsmifmt);
  }
  ///////////////////////////////////////////////////////////////////////////
  /**   Removes side-chains from linker, retains only shortest path and atoms multiply
        bonded to shortest path&#46; Must avoid matching pseudo-atoms: "[*] |$p_J$|"&#46;
        (But cxsmarts cannot represent the negation of a pseudo-atom!?)
        Kludge: First attach "protecting ring" say *1**1 at each pseudoatom
        to prevent deletion&#46;  Then run rmSideChains()&#46;  Then delete protecting ring&#46;
  */
  private static int pruneLinker(Molecule mol)
    throws SearchException,MolFormatException,MolExportException
  {
    int degree=0;
    for (MolAtom atom: mol.getAtomArray())
    {
      if (atom.getAtno()==MolAtom.PSEUDO)
      {
        MolAtom pa1 = new MolAtom(MolAtom.PSEUDO);
        MolAtom pa2 = new MolAtom(MolAtom.PSEUDO);
        mol.add(pa1);
        mol.add(pa2);
        mol.add(new MolBond(pa1,atom));
        mol.add(new MolBond(pa2,atom));
        mol.add(new MolBond(pa1,pa2));
      }
      ++degree;
    }
    int n_del=hscaf_utils.rmSideChains(mol,false);
    ArrayList<MolAtom> atoms = new ArrayList<MolAtom>();
    for (MolAtom atom: mol.getAtomArray())
    {
      if (atom.getAtno()==MolAtom.PSEUDO && atom.getBondCount()==2)
        atoms.add(atom);
    }
    for (MolAtom atom: atoms)
      mol.removeAtom(atom,MoleculeGraph.RMCLEANUP_EDGES);
    return degree;
  }
  ///////////////////////////////////////////////////////////////////////////
}
