package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;

import com.sleepycat.je.DatabaseException;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.marvin.io.MolExportException;
 

/**	Comprised of a set of hierarchical scaffolds each a {@link Scaffold} object.
	Each input molecule has an associated ScaffoldTree, a hierarchical
	tree of {@link Scaffold}s, the largest of which is the root scaffold
	(and the Bemis-Murko framework).  The ScaffoldTree object
	does not contain the {@link Scaffold} objects, which may be contained
	in a {@link ScaffoldSet} or re-constructed from
	{@link ScaffoldStore} / {@link ScaffoldStoreEntity}s or
	{@link ScaffoldDB} / {@link ScaffoldDBRecord}s .
	The ScaffoldTree object also contains the {@link Linker}s and {@link Sidechain}s
	for the input molecule, and a copy of the input molecule.
	<br>
	[Implementation notes] Here are the steps in perceiving the scaffolds
	of an input molecule:
	<ol>
	<li> ScaffoldTree constructor called by application with input molecule
	as parameter, and either a {@link ScaffoldSet}, {@link ScaffoldStore}
	or {@link ScaffoldDB} instance, which
	stores known scaffolds, for output to unique list, and to avoid
	re-analysis of previously encountered scaffolds, and thereby
	enhance performance.
	<li> The application may elect to use {@link ScaffoldStore}, which uses
	BerkeleyDB, for large-ish datasets.  {@link ScaffoldSet} in-memory storage is
	faster for small datasetsbut may exhaust available memory for large datasets.
	For large (&gt; 100k compounds) jobs, or for directly populating a
	database, {@link ScaffoldDB} can be used, but requires an available
	database server.
	<li> For each input molecule, a complete ScaffoldTree is built.
	When {@link ScaffoldSet} is used, the ScaffoldTree references {@link Scaffold} instances 
	which are contained in the {@link ScaffoldSet}.
	When {@link ScaffoldStore} is used, all {@link Scaffold} instances must be
	re-constructed from {@link ScaffoldStoreEntity} objects.
	When {@link ScaffoldDB} is used, all {@link Scaffold} instances must be
	re-constructed from {@link ScaffoldDBRecord} objects.
	<li> The perception of scaffolds from a molecule is accomplished by:
	<ul>
	<li> hscaf_utils.tagJunctions() - tags junction bonds
	<li> hscaf_utils.findChildScaffolds() - breaks one junction bond, then
	analyzes each part as candidate scaffold, then recurse  for each part until
	no more junctions. 
	</ul>
	</ol>
	<br/>
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hscaf_utils
	@see edu.unm.health.biocomp.hscaf.ScaffoldStore
	@see edu.unm.health.biocomp.hscaf.ScaffoldStoreEntity
	@see edu.unm.health.biocomp.hscaf.ScaffoldDB
	@see edu.unm.health.biocomp.hscaf.ScaffoldDBRecord
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
  /**	Create ScaffoldTree by analyzing input molecule.  All scaffolds will
	be found and associated with the root scaffold, a.k.a. the 
	Bemis-Murko framework.  If input molecule is multi-fragment,
	throw exception.   Either (1) {@link ScaffoldSet}, (2) {@link ScaffoldStore},
	or (3) {@link ScaffoldDB} used as a lookup
	to avoid re-analysis of previously analyzed scaffolds.
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
	@param scafset in-memory lookup table
	@param scafstore BerkeleyDB for lookup table
	@param scafdb RDB for lookup table
  */
  public ScaffoldTree(Molecule mol,boolean stereo,boolean keep_nitro_attachments,
        ScaffoldSet scafset,ScaffoldStore scafstore,ScaffoldDB scafdb)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException,IOException,SQLException,DatabaseException
  {
    if (mol.getFragCount(MoleculeGraph.FRAG_BASIC)>1)
      throw new ScaffoldException("Cannot analyze multi-fragment molecule.");
    this.inmol=mol.cloneMolecule();
    this.mol=mol.cloneMolecule();
    this.keep_nitro_attachments=keep_nitro_attachments;
    if (!stereo)
    {
      this.smifmt+="+0";
      this.cxsmifmt+="+0";
    }
    //System.err.println("DEBUG (ScaffoldTree): tagJunctions...");
    hscaf_utils.tagJunctions(this.mol,this.keep_nitro_attachments);
    this.rootscaf = new Scaffold(this.mol,this.keep_nitro_attachments,stereo);
    this.linkers = new ArrayList<Linker>(); //default empty
    this.sidechains = new ArrayList<Sidechain>(); //default empty
    if (this.rootscaf==null) return;

    if  (scafset!=null)
    {
      if  (scafset.containsScaffold(this.rootscaf)) // Scaf present in scafset?
      {
        // Replace rootscaf with stored instance.
        this.rootscaf=scafset.getScaffoldByID(scafset.getScaffoldID(this.rootscaf));
      }
      else // New scaf.
      {
        if (this.rootscaf.isLegal())
        {
          //System.err.println("DEBUG (ScaffoldTree): findChildScaffolds...");
          hscaf_utils.findChildScaffolds(this.rootscaf,scafset);
          //System.err.println("DEBUG (ScaffoldTree): rmJBonds...");
          rmJBonds(this.mol);

          //COMMENTED OUT DUE TO ERRORS ASSOCIATED WITH JChem 6+:
          //System.err.println("DEBUG (ScaffoldTree): findLinkers...");
          //this.linkers=findLinkers(this.mol);

          //System.err.println("DEBUG (ScaffoldTree): findSidechains...");
          this.sidechains=findSidechains(this.mol);
          //System.err.println("DEBUG (ScaffoldTree): mergeScaffoldTree...");
          scafset.mergeScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
        }
        else
        {
          this.rootscaf=null;
          return;
        }
      }
    }
    else if (scafstore!=null)
    {
      if (scafstore.scaffoldByCanSmi.contains(this.rootscaf.getCansmi())) // Scaf present in scafstore?
      {
        ScaffoldStoreEntity scent=scafstore.scaffoldByCanSmi.get(this.rootscaf.getCansmi());
        long id=scent.getId();
        this.rootscaf.setID(id);
        scafstore.populateScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
      }
      else // New scaf.
      {
        if (this.rootscaf.isLegal())
        {
          hscaf_utils.findChildScaffolds(this.rootscaf,scafstore);
          rmJBonds(this.mol);
          this.linkers=findLinkers(this.mol);
          this.sidechains=findSidechains(this.mol);
          scafstore.mergeScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
        }
        else
        {
          this.rootscaf=null;
          return;
        }
      }
    }
    else if (scafdb!=null)
    {
      if (scafdb.containsScaffoldByCansmi(this.rootscaf.getCansmi())) // Scaf present in scafdb?
      {
        ScaffoldDBRecord scafrec=scafdb.getScaffoldByCansmi(this.rootscaf.getCansmi());
        long id=scafrec.getID();
        this.rootscaf.setID(id);
        scafdb.populateScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
      }
      else // New scaf.
      {
        if (this.rootscaf.isLegal())
        {
          hscaf_utils.findChildScaffolds(this.rootscaf,scafdb);
          rmJBonds(this.mol);
          this.linkers=findLinkers(this.mol);
          this.sidechains=findSidechains(this.mol);
          scafdb.mergeScaffoldTree(this.rootscaf); // Scaf IDs assigned here.
        }
        else
        {
          this.rootscaf=null;
          return;
        }
      }
    }
    else
    {
      System.err.println("ERROR (ScaffoldTree): Should not happen; no set, store or db.");
    }
  }
  /**   {@link ScaffoldSet} used for storage.
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
	@param scafset in-memory lookup table
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldSet scafset)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException,IOException,SQLException,DatabaseException
  {
    this(mol,keep_nitro_attachments,stereo,scafset,null,null);
  }
  /**   {@link ScaffoldStore} used for storage.
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
	@param scafstore Berkeley DB for lookup table
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldStore scafstore)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException,IOException,SQLException,DatabaseException
  {
    this(mol,keep_nitro_attachments,stereo,null,scafstore,null);
  }
  /**   {@link ScaffoldDB} used for storage.
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
	@param scafdb RDB for lookup table
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo,ScaffoldDB scafdb)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException,IOException,SQLException,DatabaseException
  {
    this(mol,keep_nitro_attachments,stereo,null,null,scafdb);
  }
  /**   No ScaffoldStore, no ScaffoldSet, no ScaffoldDB.
	@param mol input molecule
	@param keep_nitro_attachments true for N-attachments scaf definition
	@param stereo stereo scaffolds (default non-stereo)
  */
  public ScaffoldTree(Molecule mol,boolean keep_nitro_attachments,boolean stereo)
    throws SearchException,MolFormatException,MolExportException,ScaffoldException,IOException,SQLException,DatabaseException
  {
    this(mol,keep_nitro_attachments,stereo,null,null,null);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Returns list of all scaffolds, deduplicated.  Note that although the child scaffolds
	for each Scaffold are unique, for a set there may be duplicates in the
	hierarchy, thus this method must deduplicate.
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
  /**	Returns list of all scaffolds, not deduplicated.  Note that although the child scaffolds
	for each Scaffold are unique, for a set there may be duplicates in the
	hierarchy.
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
  /**	SMILES representing all scaffolds as disconnected fragments.
  */
  public String getSmiForScaffoldsGroupmol(boolean show_js)
    throws MolExportException,IOException
  {
    Molecule groupmol = new Molecule();
    for (Scaffold cscaf: this.getScaffolds())
      groupmol.fuse(cscaf.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return MolExporter.exportToFormat(groupmol,cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	SMILES representing all linkers as disconnected fragments.
  */
  public String getSmiForLinkersGroupmol(boolean show_js)
    throws MolExportException,IOException
  {
    Molecule groupmol = new Molecule();
    for (Linker link: this.getLinkers())
      groupmol.fuse(link.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return MolExporter.exportToFormat(groupmol,cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	SMILES representing all side-chains as disconnected fragments.
  */
  public String getSmiForSidechainsGroupmol(boolean show_js)
    throws MolExportException,IOException
  {
    Molecule groupmol = new Molecule();
    for (Sidechain chain: this.getSidechains())
      groupmol.fuse(chain.cloneMolecule());
    if (show_js)
      hscaf_utils.replaceJHydrogensWithJPseudoatoms(groupmol);
    else
      groupmol.implicitizeHydrogens(MolAtom.ALL_H);
    return MolExporter.exportToFormat(groupmol,cxsmifmt);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Removes junction bonds in molecule previously tagged
	with hscaf_utils.tagJunctions().
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
	with hscaf_utils.tagJunctions() and fragmented with rmJBonds().
	Scaffolds identified by ring atoms.  2+ junction atoms
	and not scaf implies linker.  By process of elimination
	the rest are side-chains.  Returned sidechains 
	are deduplicated.  Note that linkers must have
	side-chains removed.

	THIS CODE FAILS WITH JChem versions 6+ ?
	E.g.

	Query: *-*
	Target: O(*1**1)*1**1 |$;J_p;;;J_p;;$|
	Caused by:
	java.lang.NullPointerException
	Mar 18, 2016 11:28:26 AM chemaxon.sss.search.SearchUtil logException
	INFO: An error occurred during search: null
  */
  private static ArrayList<Linker> findLinkers(Molecule mol)
    throws SearchException,MolFormatException,MolExportException,IOException
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
	with hscaf_utils.tagJunctions() and fragmented with rmJBonds().
	Scaffolds identified by ring atoms.  2+ junction atoms
	and not scaf implies linker.  By process of elimination
	the rest are side-chains.  Returned sidechains 
	are deduplicated.
  */
  private static ArrayList<Sidechain> findSidechains(Molecule mol)
    throws SearchException,MolFormatException,MolExportException,IOException
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
