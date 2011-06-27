package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents any single HierS scaf, which may be either
	<OL>
	<LI> a leaf scaffold in a HierS tree (no linkers) or
	<LI> a super-scaffold built from 2+ child scaffolds and 1+ linker[s]
	</OL>
	A super-scaffold may also contain one or more smaller super-scaffolds,
	hence each scaffold has zero or one parent scaffold and zero or more
	child scaffolds&#46;  The largest scaffold of a molecule is the root scaffold
	(a&#46;k&#46;a&#46; the Bemis-Murko framework)&#46;  For a single molecule
	a related set of Scaffold objects represents the HierS
	scaffold hierarchy and is implemented by class ScaffoldSet;
	similarly a scaffold hierarchy can be defined for a database&#46;
	Scaffolds are by default non-stereo but may be stereo&#46;
	<br />
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds_utils
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds
	@see edu.unm.health.biocomp.hscaf.hscaf_utils
	@author Jeremy J Yang
*/
public class Scaffold extends Molecule implements Comparable<Object>
{
  /** JChem format used for getCansmi() , for canonical scaffold smiles */
  public static String CANSMIFMT="cxsmiles:u-L-l-e-d-D-p-R-f-w-H+a_gen";
  /** JChem format used for getSmi() */
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w-a";
  /** JChem format used for getJsmi() */
  private static String cxsmifmt="cxsmiles:+L+l+e+d+D+p+R+f+w";
  /** parent scaffold in hierarchy for molecule[s] */
  private Scaffold parentscaf;
  /** child scaffold[s] in hierarchy for molecule[s] */
  private ArrayList<Scaffold> childscafs;
  /** canonical smiles used for equality comparison */
  private String cansmi;
  /** Kekule smiles (good for export) */
  private String smi;
  /** unique ID for dataset scope */
  private Integer id;
  /** flag applies to scaf and all children */
  private boolean keep_nitro_attachments;
  /** flag applies to scaf and all children */
  private boolean stereo;
  /////////////////////////////////////////////////////////////////////////////
  private Scaffold() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	When initialized from a molecule, sidechains are removed
	and explicit hydrogens are added at each junction bond&#46;
	If molecule has no rings, return empty scaffold&#46;
  */
  public Scaffold(Molecule mol,boolean keep_nitro_attachments,boolean stereo)
    throws MolFormatException,MolExportException,SearchException
  {
    if (!stereo)
    {
      this.CANSMIFMT+="+0"; //non-stereo
      this.smifmt+="+0"; //non-stereo
      this.cxsmifmt+="+0"; //non-stereo
    }
    this.parentscaf=null;
    this.id=0;
    this.childscafs = new ArrayList<Scaffold>();
    this.cansmi="";
    this.smi="";
    this.keep_nitro_attachments=keep_nitro_attachments;
    this.stereo=stereo;
    MolSearch patR = new MolSearch(); //ring matcher
    patR.setQuery(MolImporter.importMol("[R]","smarts:"));
    patR.setTarget(mol);
    boolean ok=patR.isMatching();
    if (ok)
    {
      mol.clonecopy(this);
      ok=this.dearomatize(); //avoid auto-H-adjustments by JChem
    }
    if (ok)
    {
      hscaf_utils.rmSideChains(this,keep_nitro_attachments);
      this.aromatize(MoleculeGraph.AROM_GENERAL);
      this.cansmi=this.exportToFormat(CANSMIFMT);
      this.smi=this.exportToFormat(smifmt);
    }
    else
      (new Molecule()).clonecopy(this); //empty molecule
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Copy constructor&#46;
  */
  public Scaffold(Scaffold scaf)
  {
    scaf.clonecopy(this);
    this.cansmi=scaf.cansmi;
    this.smi=scaf.smi;
    this.childscafs = new ArrayList<Scaffold>();
    this.parentscaf=null;
    this.id=scaf.id;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	For sorting by size&#46;
  */
  public int compareTo(Object o) throws ClassCastException
  {
    return ((Scaffold)o).getAtomCount()-this.getAtomCount();
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean isRoot() { return (this.parentscaf==null); }
  /////////////////////////////////////////////////////////////////////////////
  public boolean isLeaf() { return (this.childscafs.size()==0); }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(Scaffold scaf2)
  {
    return (this.cansmi.equals(scaf2.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Reject (1) benzene, and (2) empty scaffolds &#46;
	Note that benzene is identified by canonical smiles, unfortunately 
	hard-coded&#46;
  */
  public boolean isLegal()
  {
    if (this.getAtomCount()==0) return false;
    if (this.getCansmi().equals("c1ccccc1")) return false; //benzene disallowed
    if (this.getCansmi().equals("C1=CC=CC=C1")) return false; //benzene disallowed
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Adds child scaffold if new and legal&#46;
  */
  private boolean addChild(Scaffold scaf2)
  {
    if (!scaf2.isLegal()) return false;
    boolean is_present=false;
    for (Scaffold cscaf: this.getChildScaffolds())
    {
      if (cscaf.equals(scaf2)) { is_present=true; break; }
    }
    if (!is_present) {
      this.childscafs.add(scaf2);
      return true;
    }
    return false;
  }
  /////////////////////////////////////////////////////////////////////////////
  public void setParentScaffold(Scaffold scaf2) { this.parentscaf=scaf2; }
  public Scaffold getParentScaffold() { return this.parentscaf; }
  ///////////////////////////////////////////////////////////////////////////
  public Scaffold getRootScaffold()
  {
    if (this.isRoot()) return this;
    else return (this.getParentScaffold().getRootScaffold());
  }
  /////////////////////////////////////////////////////////////////////////////
  public int getChildCount() { return this.childscafs.size(); }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get unique ID for dataset scope.
  */
  public int getID() { return this.id; }
  /**	Set unique ID for dataset scope.
  */
  public void setID(int id2) { this.id=id2; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get canonical SMILES for this scaffold (as used for equality comparison)&#46;
	Canonical SMILES are generated using JChem format Scaffold.CANSMIFMT&#46;
  */
  public String getCansmi() { return this.cansmi; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get Kekule SMILES for this scaffold&#46;  Recommended for export
	as resulting SMILES will be more universally compatible with
	other software&#46;
  */
  public String getSmi() { return this.smi; }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this scaffold with junction atoms as pseudo-atoms&#46;
	Note this is in ChemAxon cxsmiles format e&#46;g&#46;
	"*C1CC1* |$p_J;;;;p_J$|"&#46;
  */
  public String getJsmi()
    throws MolExportException
  {
    Molecule scafmol=this.cloneMolecule();
    hscaf_utils.replaceJHydrogensWithJPseudoatoms(scafmol);
    return scafmol.exportToFormat(cxsmifmt);
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	This recursive algorithm finds all child scafs&#46;  For each junction bond,
	copy scaffold and cut bond, resulting in two parts which may be valid
	scaffolds (sidechains have been deleted, but one part may comprise
	invalid scaffold benzene)&#46; Special case is junction bonds to hydrogens,
	which denote attachment points&#46; Then recurse and find grandchild scaffolds
	for each child scafs&#46;  For simple 2-way linkers it is redundant to cut both
	junctions of a linker, but for (3+)-way linkers this allows enumeration of all
	combinations&#46;  For each scaffold found, if unique (via canonical smiles),
	add new child and specify its parentage&#46; If no junction bonds exist,
	findChildScaffolds() returns null, recursion terminates, and the current
	scaffold is a leaf&#46;
	<br />
	This method does the heavy-lifting&#46;
  */
  public int findChildScaffolds()
    throws MolExportException,SearchException,MolFormatException
  {
    int n_cscafs=0;
    for (MolBond bond: this.getBondArray())
    {
      if (bond.getSetSeq()>0)
      {
        int bidx=this.indexOf(bond); //Note: bond idxs must be same in mol clone.
        Molecule scafmol = this.cloneMolecule();
        MolBond jbond=scafmol.getBond(bidx);
        MolAtom a1=jbond.getAtom1();
        MolAtom a2=jbond.getAtom2();
        if (a1.getAtno()==1 || a2.getAtno()==1) continue; //junction-atoms
        MolAtom h_new1 = new MolAtom(1);
        MolAtom h_new2 = new MolAtom(1);
        scafmol.add(h_new1);
        scafmol.add(h_new2);
        MolBond b_new1 = new MolBond(h_new1,a1,jbond.getType());
        MolBond b_new2 = new MolBond(h_new2,a2,jbond.getType());
        scafmol.add(b_new1);
        scafmol.add(b_new2);
        b_new1.setSetSeq(jbond.getSetSeq());
        b_new2.setSetSeq(jbond.getSetSeq());
        scafmol.removeBond(jbond);
        for (Molecule partmol: scafmol.convertToFrags()) //should be 2
        {
          Scaffold cscaf = new Scaffold(partmol,this.keep_nitro_attachments,this.stereo);
          // Check: has scaf been seen already in this tree?  If yes save time.
          Scaffold foundscaf = this.getRootScaffold().findChildScaffold(cscaf);
          if (foundscaf!=null)
          {
            this.addChild(foundscaf);
            ++n_cscafs;
            n_cscafs+=foundscaf.getAllChildCount();
          }
          else if (this.addChild(cscaf))
          {
            ++n_cscafs;
            cscaf.setParentScaffold(this);
            n_cscafs+=cscaf.findChildScaffolds(); //size=0 if leaf
          }
        }
      }
    }
    return n_cscafs;
  }
  ///////////////////////////////////////////////////////////////////////////
  public ArrayList<Scaffold> getChildScaffolds()
  {
    return this.childscafs;
  }
  ///////////////////////////////////////////////////////////////////////////
  public int getAllChildCount()
  {
    if (this.childscafs==null) return 0;
    int n_cscaf=this.childscafs.size();
    for (Scaffold cscaf: this.childscafs)
      n_cscaf+=cscaf.getAllChildCount();
    return n_cscaf;
  }
  ///////////////////////////////////////////////////////////////////////////
  public ArrayList<Scaffold> getAllChildScaffolds()
  {
    ArrayList<Scaffold> scafs = new ArrayList<Scaffold>();
    if (this.childscafs==null) return scafs;
    scafs.addAll(this.childscafs);
    for (Scaffold cscaf: this.childscafs)
    {
      scafs.addAll(cscaf.getAllChildScaffolds());
    }
    return scafs;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Search for child scaffold equal to the specified scaffold; if not
	found return null&#46;
  */
  public Scaffold findChildScaffold(Scaffold scaf)
  {
    if (this.childscafs==null) return null;
    for (Scaffold cscaf: this.childscafs)
    {
      if (cscaf.equals(scaf)) return cscaf;
    }
    return null;
  }
  ///////////////////////////////////////////////////////////////////////////
}
