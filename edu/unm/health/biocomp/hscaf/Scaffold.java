package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.marvin.io.MolExportException;
 
/**	Represents any single HierS scaffold, which may contain zero,
	or two or more child scaffolds.  A scaffold with zero child
	scaffolds is called a leaf scaffold.  Each scaffold has a unique
	set of child scaffolds.  However, a scaffold may have multiple
	parent scaffolds, in a scaffold tree context or a scaffold
	database context, since the same scaffold may be present in
	multiple locations within one molecule or multiple molecules.
	The largest scaffold of a molecule is the root scaffold
	(a.k.a. the Bemis-Murko framework).  For a single molecule
	a related set of Scaffold objects represents the HierS
	scaffold hierarchy; similarly a scaffold hierarchy can be defined for a database.
	Scaffolds are by default non-stereo but may be stereo.
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
  /** JChem format used for getCansmi() , for canonical scaffold smiles.
      E.g. cansmi=MolExporter.exportToFormat(mol,Scaffold.CANSMIFMT)
  */
  public static final String CANSMIFMT="cxsmiles:u-L-l-e-d-D-p-R-f-w-H+a_gen+0";
  /** JChem format used for getCansmi() , for canonical scaffold smiles.
      E.g. cansmi=MolExporter.exportToFormat(mol,Scaffold.CANSMIFMT_STEREO)
  */
  public static final String CANSMIFMT_STEREO="cxsmiles:u-L-l-e-d-D-p-R-f-w-H+a_gen";
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
  private Long id;
  /** flag applies to scaf and all children */
  private boolean keep_nitro_attachments;
  /** flag applies to scaf and all children */
  private boolean stereo;
  /////////////////////////////////////////////////////////////////////////////
  private Scaffold() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	When initialized from a molecule, sidechains are removed
	and explicit hydrogens are added at each junction bond.
	If molecule has no rings, return empty scaffold.
  */
  public Scaffold(Molecule mol,boolean stereo,boolean keep_nitro_attachments)
    throws MolFormatException
  {
    if (!stereo)
    {
      this.smifmt+="+0"; //non-stereo
      this.cxsmifmt+="+0"; //non-stereo
    }
    this.parentscaf=null;
    this.id=0L;
    this.childscafs=null;
    this.cansmi=null; // Lazily evaluated
    this.smi=null; // Lazily evaluated
    this.keep_nitro_attachments=keep_nitro_attachments;
    this.stereo=stereo;
    int [][] sssr = mol.getSSSR();
    boolean ok=(sssr.length>0);
    if (ok)
    {
      mol.clonecopy(this);
      ok=this.dearomatize(); //avoid auto-H-adjustments by JChem
    }
    if (ok)
    {
      hscaf_utils.rmSideChains(this,keep_nitro_attachments);
      this.aromatize(MoleculeGraph.AROM_GENERAL);
    }
    else
      (new Molecule()).clonecopy(this); //empty molecule
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Initialized molecule from smiles and call usual constructor.
  */
  public Scaffold(String smiles,Boolean stereo,Boolean keep_nitro_attachments)
    throws MolFormatException
  {
    this(new Molecule(),stereo,keep_nitro_attachments);
    try { MolImporter.importMol(smiles.getBytes(),"smiles:",null,this); }
    catch (MolFormatException e) { } // should not happen!
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Copy constructor.
  */
  public Scaffold(Scaffold scaf)
  {
    scaf.clonecopy(this);
    this.cansmi=scaf.cansmi;
    this.smi=scaf.smi;
    this.childscafs=null;
    this.parentscaf=null;
    this.id=scaf.id;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	For sorting by size.
  */
  public int compareTo(Object o) throws ClassCastException
  {
    return ((Scaffold)o).getAtomCount()-this.getAtomCount();
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean isRoot() { return (this.parentscaf==null); }
  /////////////////////////////////////////////////////////////////////////////
  public boolean isLeaf()
  {
    return (this.childscafs==null || this.childscafs.size()==0);
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(Scaffold scaf2)
  {
    return (this.getCansmi().equals(scaf2.getCansmi()));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Reject (1) benzene, and (2) empty scaffolds .
	Note that benzene is identified by canonical smiles, unfortunately 
	hard-coded.
  */
  public boolean isLegal()
  {
    if (this.isEmpty()) this.decompress(); // isEmpty() from Molecule
    if (this.getAtomCount()==0) return false;
    if (this.getCansmi().equals("c1ccccc1")) return false; //benzene disallowed
    if (this.getCansmi().equals("C1=CC=CC=C1")) return false; //benzene disallowed
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Adds child scaffold if new and legal.
  */
  public boolean addChild(Scaffold scaf2)
  {
    if (this.childscafs==null) this.childscafs = new ArrayList<Scaffold>();
    if (!scaf2.isLegal()) return false;
    boolean is_already_child=false;
    for (Scaffold cscaf: this.getChildScaffolds())
    {
      if (cscaf.equals(scaf2)) { is_already_child=true; break; }
    }
    if (is_already_child) return false;
    this.childscafs.add(scaf2);
    return true;
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
  public int getChildCount()
  {
    if (this.childscafs==null) return 0;
    else return this.childscafs.size();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get unique ID for dataset scope.
  */
  public long getID() { return this.id; }
  /**	Set unique ID for dataset scope.
  */
  public void setID(long id2) { this.id=id2; }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get canonical SMILES for this scaffold (as used for equality comparison).
	Canonical SMILES are generated using JChem format Scaffold.CANSMIFMT.
	Lazy evaluation.
  */
  public String getCansmi()
  {
    if (this.cansmi==null)
    {
      try { this.cansmi=MolExporter.exportToFormat(this,(this.stereo?this.CANSMIFMT_STEREO:this.CANSMIFMT)); }
      catch (MolExportException e) { this.cansmi=""; }
      catch (IOException e) { this.cansmi=""; }
    }
    return this.cansmi;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Get Kekule SMILES for this scaffold.  Recommended for export
	as resulting SMILES will be more universally compatible with
	other software. Lazy evaluation.
  */
  public String getSmi()
  {
    if (this.smi==null)
    {
      this.decompress();
      try { this.smi=MolExporter.exportToFormat(this,this.smifmt); }
      catch (MolExportException e) { this.smi=""; }
      catch (IOException e) { this.smi=""; }
    }
    return this.smi;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Get SMILES for this scaffold with junction atoms as pseudo-atoms.
	Note this is in ChemAxon cxsmiles format e.g.
	"*C1CC1* |$p_J;;;;p_J$|".
  */
  public String getJsmi()
  {
    this.decompress();
    Molecule scafmol=this.cloneMolecule();
    hscaf_utils.replaceJHydrogensWithJPseudoatoms(scafmol);
    String jsmi=null;
    try { jsmi=MolExporter.exportToFormat(scafmol,cxsmifmt); }
    catch (MolExportException e) { jsmi=""; }
    catch (IOException e) { jsmi=""; }
    return jsmi;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	This recursive algorithm finds all child scafs.  For each junction bond,
	copy scaffold and cut bond, resulting in two parts which may be valid
	scaffolds (sidechains have been deleted, but one part may comprise
	invalid scaffold benzene). Special case is junction bonds to hydrogens,
	which denote attachment points. Then recurse and find grandchild scaffolds
	for each child scafs.  For simple 2-way linkers it is redundant to cut both
	junctions of a linker, but for (3+)-way linkers this allows enumeration of all
	combinations.  For each scaffold found, if unique (via canonical smiles),
	add new child and specify its parentage. If no junction bonds exist,
	findChildScaffolds() returns null, recursion terminates, and the current
	scaffold is a leaf.
	<br />
	This method does the heavy-lifting.
  */
  public int findChildScaffolds()
    throws SearchException,MolFormatException
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
    if (this.childscafs==null) this.childscafs = new ArrayList<Scaffold>();
    return this.childscafs;
  }
  ///////////////////////////////////////////////////////////////////////////
  public ArrayList<Long> getChildIDs()
  {
    ArrayList<Long> childids = new ArrayList<Long>();
    for (Scaffold cscaf: this.getChildScaffolds())
      childids.add(cscaf.getID());
    return childids;
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
	found return null.
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
  /**   Generates string representing the hierarchical scaffold sub tree
	rooted by this scaffold.
        e.g. "1:(2,3)" or "1:(2:(3,4,5),6:(4,7))"
  */
  public String subTreeAsString()
  {
    if (this==null) return "";
    String str=""+this.getID();
    if (this.getChildCount()>0)
    {
      str+=":(";
      int i=0;
      for (Scaffold cscaf: this.getChildScaffolds())
      {
        if (i>0) str+=",";
        str+=cscaf.subTreeAsString();
        ++i;
      }
      str+=")";
    }
    return str;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	For storage; delete molecule object and retain canonical SMILES.
  	Canonical SMILES used since it is an identifier.
  */
  public void compress()
  {
    //System.err.println("DEBUG: Scaffold.compress: [id="+this.id+"]");
    if (this.cansmi==null) this.cansmi=this.getCansmi();
    this.smi=null;
    if (!this.isEmpty()) this.clear();
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Reconstruct molecule object from stored canonical SMILES.
  */
  public void decompress()
  {
    if (!this.isEmpty()) return; // already decompressed
    if (this.cansmi==null)
    {
      System.err.println("ERROR (Scaffold.decompress): [id="+this.id+"] No cansmi.  Should not happen!");
    }
    try { MolImporter.importMol(this.cansmi.getBytes(),"smiles:",null,this); }
    catch (MolFormatException e) { } // should not happen!
  }
  ///////////////////////////////////////////////////////////////////////////
  public Boolean isStereo()
  {
    return this.stereo;
  }
  ///////////////////////////////////////////////////////////////////////////
  public Boolean isKeep_nitro_attachments()
  {
    return this.keep_nitro_attachments;
  }
  ///////////////////////////////////////////////////////////////////////////
}
