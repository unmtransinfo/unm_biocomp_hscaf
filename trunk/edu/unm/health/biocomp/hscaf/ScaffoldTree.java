package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import com.sleepycat.je.DatabaseException;
import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.marvin.io.MolExportException;
 
/////////////////////////////////////////////////////////////////////////////
/**	Comprised of a set of hierarchical scaffolds each a Scaffold object&#46;
	Each input molecule has an associated ScaffoldTree, a hierarchical
	tree of Scaffolds, the largest of which is the root scaffold
	(and the Bemis-Murko framework)&#46;  The ScaffoldTree object
	does not contain the Scaffold objects, which may be contained
	in a ScaffoldSet or re-constructed from ScaffoldStore
	ScaffoldEntitys&#46;
	The ScaffoldTree object also contains the Linkers and Sidechains
	for the input molecule, and a copy of the input molecule&#46;
	<br />
	[Implementation notes] Here are the steps in perceiving the scaffolds
	of an input molecule:
	<ol>
	<li> ScaffoldTree constructor called by application with input molecule
	as parameter, and either a ScaffoldSet or ScaffoldStore instance, which
	stores known scaffolds, for output to unique list, and to avoid
	re-analysis of previously encountered scaffolds, and thereby
	enhance performance&#46;
	<li> The application may elect to use ScaffoldStore, which uses
	BerkeleyDB, for large datasets&#46;  ScaffoldSet in-memory storage is
	faster but may exhaust available memory for large datasets&#46;
	<li> For each input molecule, a complete ScaffoldTree is built&#46;
	When ScaffoldSet is used, the ScaffoldTree references Scaffold instances 
	which are contained in the ScaffoldSet&#46;  When ScaffoldStore is used,
	all Scaffold instances must be re-constructed from ScaffoldEntity
	objects&#46;
	<li> The perception of scaffolds from a molecule is accomplished by:
	<ul>
	<li> hscaf_utils.tagJunctions() - tags junction bonds
	<li> hscaf_utils.findChildScaffolds() - breaks one junction bond, then
	analyzes each part as candidate scaffold, then recurse  for each part until
	no more junctions&#46; 
	</ul>
	</ol>
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hscaf_utils
	@see edu.unm.health.biocomp.hscaf.ScaffoldStore
	@see edu.unm.health.biocomp.hscaf.ScaffoldEntity
	@author Jeremy J Yang
*/
public class ScaffoldTree
{
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static String cxsmifmt="cxsmiles:+L+l+e+d+D+p+R+f+w";
  /** input molecule (not changed) */
  public Molecule inmol;
  private Molecule mol; // working molecule (tagged, etc.)
  private boolean keep_nitro_attachments;
  private Scaffold rootscaf;
  private ArrayList<Linker> linkers;
  private ArrayList<Sidechain> sidechains;
  /////////////////////////////////////////////////////////////////////////////
  private ScaffoldTree() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	Create ScaffoldTree by analyzing input molecule&#46;  All scaffolds will
	be found and associated with the root scaffold, a&#46;k&#46;a&#46; the 
	Bemis-Murko framework&#46;  If input molecule is multi-fragment,
	throw exception&#46;   If scafstore is provided, it will be used as a lookup
	to avoid re-analysis of previously analyzed scaffolds&#46;
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
	@param scafstore Berkeley DB for lookup table
  */
  public ScaffoldTree(Molecule mol,boolean stereo,boolean keep_nitro_attachments,
        ScaffoldSet scafset,ScaffoldStore scafstore)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException
  {
    //System.err.println("DEBUG: (ScaffoldTree) ...");
    if (mol.getFragCount()>1)
      throw new ScaffoldException("Cannot analyze multi-fragment molecule.");
    this.inmol=mol.cloneMolecule();
    this.mol=mol.cloneMolecule();
    this.keep_nitro_attachments=keep_nitro_attachments;
    if (!stereo)
    {
      this.smifmt+="+0";
      this.cxsmifmt+="+0";
    }
    hscaf_utils.tagJunctions(this.mol,this.keep_nitro_attachments);
    this.rootscaf = new Scaffold(this.mol,this.keep_nitro_attachments,stereo);
    this.linkers = new ArrayList<Linker>(); //default empty
    this.sidechains = new ArrayList<Sidechain>(); //default empty
    if (this.rootscaf==null) return;
    try
    {
      if (scafstore!=null && scafstore.scaffoldByCanSmi.contains(this.rootscaf.getCansmi()))
      {
	ScaffoldEntity scent = scafstore.scaffoldByCanSmi.get(this.rootscaf.getCansmi());
        long id = scent.getId();
        this.rootscaf.setID(id);
      }
      else if  (scafset!=null && scafset.containsScaffold(this.rootscaf)) // Scaf present in scafset?
      {
        // Replace root scaffold with stored instance.
        this.rootscaf=scafset.getScaffoldByID(scafset.getScaffoldID(this.rootscaf));
      }
      else
      {
        if (this.rootscaf.isLegal())
        {
          hscaf_utils.findChildScaffolds(this.rootscaf,scafset,scafstore);
          //System.err.println("DEBUG: (ScaffoldTree) back from findChildScaffolds  ...");
          rmJBonds(this.mol);
          this.linkers=findLinkers(this.mol);
          this.sidechains=findSidechains(this.mol);
        }
        else
          this.rootscaf=null;
      }
      if (this.rootscaf!=null)
      {
        if (scafset!=null) scafset.mergeScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
        if (scafstore!=null) scafstore.mergeScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
      }
    }
    catch (DatabaseException e)
    {
      System.err.println("Fatal database exception: "+e);
      System.exit(1);
    }
    //System.err.println("DEBUG: (ScaffoldTree) leaving ...");
  }
  /**   No ScaffoldStore&#46;
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldSet scafset)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException
  {
    this(mol,keep_nitro_attachments,stereo,scafset,null);
  }
  /**   No ScaffoldSet&#46;
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldStore scafstore)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException
  {
    this(mol,keep_nitro_attachments,stereo,null,scafstore);
  }
  /**   No ScaffoldStore, no ScaffoldSet&#46;
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException
  {
    this(mol,keep_nitro_attachments,stereo,null,null);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns list of all scaffolds, deduplicated&#46;  Note that although the child scaffolds
	for each Scaffold are unique, for a set there may be duplicates in the
	hierarchy, thus this method must deduplicate&#46;
  */
  public ArrayList<Scaffold> getScaffolds()
  {
    ArrayList<Scaffold> scafs = new ArrayList<Scaffold>();
    ArrayList<String> usmis = new ArrayList<String>(); //should be HashSet<String,Boolean>
    if (this.rootscaf==null) return scafs;
    scafs.add(this.rootscaf);
    usmis.add(this.rootscaf.getCansmi());
    for (Scaffold cscaf: this.rootscaf.getAllChildScaffolds())
    {
      if (usmis.contains(cscaf.getCansmi())) continue;
      usmis.add(cscaf.getCansmi());
      scafs.add(cscaf);
    }
    return scafs;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns list of all scaffolds, not deduplicated&#46;  Note that although the child scaffolds
	for each Scaffold are unique, for a set there may be duplicates in the
	hierarchy&#46;
  */
  public ArrayList<Scaffold> getAllScaffolds()
  {
    ArrayList<Scaffold> scafs = new ArrayList<Scaffold>();
    if (this.rootscaf==null) return scafs;
    scafs.add(this.rootscaf);
    scafs.addAll(this.rootscaf.getAllChildScaffolds());
    return scafs;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns count of all scaffolds, deduplicated.
  */
  public int getScaffoldCount()
  {
    if (this.rootscaf==null) return 0;
    return this.getScaffolds().size();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns count of all leaf scaffolds, deduplicated.
  */
  public int getLeafScaffoldCount()
  {
    if (this.rootscaf==null) return 0;
    int n_scaf=0;
    for (Scaffold scaf: this.getScaffolds())
    {
      if (scaf.getChildCount()==0) ++n_scaf;
    }
    return n_scaf;
  }
  /////////////////////////////////////////////////////////////////////////////
  public ArrayList<Scaffold> getLeafScaffolds()
  {
    ArrayList<Scaffold> scafs = new ArrayList<Scaffold>();
    if (this.rootscaf==null) return scafs;
    for (Scaffold scaf: this.getScaffolds())
    {
      if (scaf.getChildCount()==0) scafs.add(scaf);
    }
    return scafs;
  }
  /////////////////////////////////////////////////////////////////////////////
  public Scaffold getRootScaffold() { return this.rootscaf; }
  /////////////////////////////////////////////////////////////////////////////
  public ArrayList<Linker> getLinkers()
  {
    return this.linkers;
  }
  public ArrayList<Sidechain> getSidechains()
  {
    return this.sidechains;
  }
  public int getLinkerCount()
  {
    return this.linkers.size();
  }
  public int getSidechainCount()
  {
    return this.sidechains.size();
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	SMILES representing all scaffolds as disconnected fragments&#46;
  */
  public String getSmiForScaffoldsGroupmol(boolean show_js)
    throws MolExportException
  {
    Molecule groupmol = new Molecule();
    for (Scaffold cscaf: this.getScaffolds())
      groupmol.fuse(cscaf.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return groupmol.exportToFormat(cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	SMILES representing all linkers as disconnected fragments&#46;
  */
  public String getSmiForLinkersGroupmol(boolean show_js)
    throws MolExportException
  {
    Molecule groupmol = new Molecule();
    for (Linker link: this.getLinkers())
      groupmol.fuse(link.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return groupmol.exportToFormat(cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	SMILES representing all side-chains as disconnected fragments&#46;
  */
  public String getSmiForSidechainsGroupmol(boolean show_js)
    throws MolExportException
  {
    Molecule groupmol = new Molecule();
    for (Sidechain chain: this.getSidechains())
      groupmol.fuse(chain.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return groupmol.exportToFormat(cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Removes junction bonds in molecule previously tagged
	with hscaf_utils.tagJunctions()&#46;
  */
  private static int rmJBonds(Molecule mol)
  {
    int n_del=0;
    for (MolBond bond: mol.getBondArray())
    {
      if (bond.getSetSeq()>0)
      {
        MolAtom a1=bond.getAtom1();
        MolAtom a2=bond.getAtom2();
        //if (a1.getAtno()==1 || a2.getAtno()==1) continue; //junction-atoms
        MolAtom ja1 = new MolAtom(MolAtom.PSEUDO);
        MolAtom ja2 = new MolAtom(MolAtom.PSEUDO);
        ja1.setAliasstr("J");
        ja2.setAliasstr("J");
        mol.add(ja1);
        mol.add(ja2);
        MolBond b_new1 = new MolBond(ja1,a1,bond.getType());
        MolBond b_new2 = new MolBond(ja2,a2,bond.getType());
        mol.add(b_new1);
        mol.add(b_new2);
        b_new1.setSetSeq(bond.getSetSeq());
        b_new2.setSetSeq(bond.getSetSeq());
        mol.removeBond(bond);
        ++n_del;
      }
    }
    return n_del;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Finds linkers in fragmented molecule previously tagged
	with hscaf_utils.tagJunctions() and fragmented with rmJBonds()&#46;
	Scaffolds identified by ring atoms&#46;  2+ junction atoms
	and not scaf implies linker&#46;  By process of elimination
	the rest are side-chains&#46;  Returned sidechains 
	are deduplicated&#46;  Note that linkers must have
	side-chains removed&#46;
  */
  private static ArrayList<Linker> findLinkers(Molecule mol)
    throws SearchException,MolFormatException,MolExportException
  {
    ArrayList<Linker> links = new ArrayList<Linker>();
    MolSearch patScaf = new MolSearch(); //scaf pattern
    patScaf.setQuery(MolImporter.importMol("[R]","smarts:"));
    MolSearch patLink = new MolSearch(); //linker pattern (2+ Js)
    patLink.setQuery(MolImporter.importMol("[*].[*] |$J_p;J_p$|","cxsmarts:"));
    ArrayList<String> usmis = new ArrayList<String>();
    for (Molecule fragmol: mol.cloneMolecule().convertToFrags())
    {
      patScaf.setTarget(fragmol);
      patLink.setTarget(fragmol);
      if (patScaf.isMatching()) continue;
      else if (patLink.isMatching())
      {
        Linker newlink = new Linker(fragmol);
        if (usmis.contains(newlink.getCansmi())) continue;
        else
        {
          usmis.add(newlink.getCansmi());
          links.add(newlink);
        }
      }
      else continue;
    }
    return links;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Finds side-chains in fragmented molecule previously tagged
	with hscaf_utils.tagJunctions() and fragmented with rmJBonds()&#46;
	Scaffolds identified by ring atoms&#46;  2+ junction atoms
	and not scaf implies linker&#46;  By process of elimination
	the rest are side-chains&#46;  Returned sidechains 
	are deduplicated&#46;
  */
  private static ArrayList<Sidechain> findSidechains(Molecule mol)
    throws SearchException,MolFormatException,MolExportException
  {
    ArrayList<Sidechain> chains = new ArrayList<Sidechain>();
    MolSearch patScaf = new MolSearch(); //scaf pattern
    patScaf.setQuery(MolImporter.importMol("[R]","smarts:"));
    MolSearch patLink = new MolSearch(); //linker pattern (2+ Js)
    patLink.setQuery(MolImporter.importMol("[*].[*] |$J_p;J_p$|","cxsmarts:"));
    ArrayList<String> usmis = new ArrayList<String>();
    for (Molecule fragmol: mol.cloneMolecule().convertToFrags())
    {
      patScaf.setTarget(fragmol);
      patLink.setTarget(fragmol);
      if (patScaf.isMatching()) continue;
      else if (patLink.isMatching()) continue;
      else
      {
        Sidechain newchain = new Sidechain(fragmol);
        if (usmis.contains(newchain.getCansmi())) continue;
        else
        {
          usmis.add(newchain.getCansmi());
          chains.add(newchain);
        }
      }
    }
    return chains;
  }
  ///////////////////////////////////////////////////////////////////////////
  public String toString()
  {
    if (this.getRootScaffold()==null) return "";
    String str="";
    return this.getRootScaffold().subTreeAsString();
  }
  ///////////////////////////////////////////////////////////////////////////
}
