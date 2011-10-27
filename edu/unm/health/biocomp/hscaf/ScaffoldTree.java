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
/**	Comprised of a set of hierarchical scaffolds each a Scaffold
	object&#46;
	<br />
	To do:
	[ ] If scaffold already in tree, save CPU, do not re-analyze, copy sub-tree
	or maybe cross-link&#46;  Maybe eventually avoid re-analysis for entire
	dataset&#46;
	<br />
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hscaf_utils
	@see edu.unm.health.biocomp.hscaf.ScaffoldStore
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
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldStore scafstore)
    throws SearchException,MolFormatException,MolExportException,ScaffoldTreeException
  {
    if (mol.getFragCount()>1)
      throw new ScaffoldTreeException("Cannot analyze multi-fragment molecule.");
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
    try
    {
      if (scafstore!=null && scafstore.scaffoldByCanSmi.contains(this.rootscaf.getCansmi()))
      {
        this.rootscaf.setID(scafstore.scaffoldByCanSmi.get(this.rootscaf.getCansmi()).getId());
      }
      else
      {
        if (this.rootscaf!=null)
        {
          if (this.rootscaf.isLegal())
          {
            this.rootscaf.findChildScaffolds();
            rmJBonds(this.mol);
            this.linkers=findLinkers(this.mol);
            this.sidechains=findSidechains(this.mol);
          }
          else
            this.rootscaf=null;
        }
      }
      if (scafstore!=null) scafstore.mergeScaffoldTree(this); // Scaf IDs assigned here.
    }
    catch (DatabaseException e)
    {
      System.err.println("Fatal database exception: "+e);
      System.exit(1);
    }
  }
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo)
    throws SearchException,MolFormatException,MolExportException,ScaffoldTreeException
  {
    this(mol,keep_nitro_attachments,stereo,null);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns list of all scaffolds, deduplicated&#46;  Note that although the child scaffolds
	for each Scaffold are unique, for a set there may be duplicates in the
	hierarchy, thus this method must deduplicate&#46;
  */
  public ArrayList<Scaffold> getScaffolds()
  {
    ArrayList<Scaffold> scafs = new ArrayList<Scaffold>();
    ArrayList<String> usmis = new ArrayList<String>();
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
}
