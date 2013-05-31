package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.marvin.io.MolExportException;
import chemaxon.calculations.TopologyAnalyser;

import com.sleepycat.je.DatabaseException;
 
/**	Contains internal static functions used in hierarchical scaffold analysis. 
	Not normally recommended for use by applications.
	<br />
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.ScaffoldSet
	@see edu.unm.health.biocomp.hscaf.ScaffoldStore
	@see edu.unm.health.biocomp.hscaf.ScaffoldDB
	@author Jeremy J Yang
*/
public class hscaf_utils
{
  private static final String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private hscaf_utils() {} //disallow default constructor
  private static final int HSCAF_MAX_COMPUTE_TIME=600;	//Max compute time in sec.
  ///////////////////////////////////////////////////////////////////////////
  public static int findChildScaffolds(Scaffold scaf,ScaffoldSet scafset)
    throws SearchException,MolFormatException,DatabaseException,SQLException,IOException
  {
    return findChildScaffolds(scaf,scafset,null,null);
  }
  ///////////////////////////////////////////////////////////////////////////
  public static int findChildScaffolds(Scaffold scaf,ScaffoldStore scafstore)
    throws SearchException,MolFormatException,DatabaseException,SQLException,IOException
  {
    return findChildScaffolds(scaf,null,scafstore,null);
  }
  ///////////////////////////////////////////////////////////////////////////
  public static int findChildScaffolds(Scaffold scaf,ScaffoldDB scafdb)
    throws SearchException,MolFormatException,DatabaseException,SQLException,IOException
  {
    return findChildScaffolds(scaf,null,null,scafdb);
  }
  ///////////////////////////////////////////////////////////////////////////
  /**   This recursive algorithm finds all child scaffolds.  For each junction bond,
        copy scaffold and cut bond, resulting in two parts which may be valid
        scaffolds (sidechains have been deleted, but one part may comprise
        invalid scaffold benzene). Special case is junction bonds to hydrogens,
        which denote attachment points. Then recurse and find grandchild scaffolds
        for each child scaffold.  For simple 2-way linkers it is redundant to cut both
        junctions of a linker, but for (3+)-way linkers this allows enumeration of all
        combinations.  For each scaffold found, if unique (via canonical smiles),
        add new child and specify its parentage. If no junction bonds exist,
        findChildScaffolds() returns null, recursion terminates, and the current
        scaffold is a leaf. The heavy-lifting goes through this method.
	ScaffoldSet and ScaffoldStore store known scaffolds to avoid costly
	re-analysis.  With ScaffoldSet, if scaffold already known, use
	existing scaffold object, not new child candidate object.  With
	ScaffoldStore, new object must be used.
  */
  public static int findChildScaffolds(Scaffold scaf,ScaffoldSet scafset,ScaffoldStore scafstore,ScaffoldDB scafdb)
    throws SearchException,MolFormatException,DatabaseException,SQLException,IOException
  {
    int n_cscafs=0;
    //System.err.println("DEBUG: (findChildScaffolds) ...");
    for (MolBond bond: scaf.getBondArray())
    {
      if (bond.getSetSeq()>0)
      {
        int bidx=scaf.indexOf(bond); //Note: bond idxs must be same in mol clone.
        Molecule scafmol = scaf.cloneMolecule();
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
          Scaffold cscaf = new Scaffold(partmol,scaf.isStereo(),scaf.isKeep_nitro_attachments());
          // Check: has scaf been seen already in this tree?  If yes save time, use foundscaf.
          Scaffold foundscaf=scaf.getRootScaffold().findChildScaffold(cscaf);
          if (foundscaf!=null)
          {
            boolean ok=scaf.addChild(foundscaf);
            ++n_cscafs;
            n_cscafs+=foundscaf.getAllChildCount();
          }
          // Check: is scaf already in scafset?  If yes save time, use foundscaf.
          else if (scafset!=null && scafset.containsScaffold(cscaf))
          {
            foundscaf=scafset.findScaffold(cscaf);
            //System.err.println("DEBUG: (findChildScaffolds) found in scafset: "+cscaf.getCansmi());
            boolean ok=scaf.addChild(foundscaf);
            //System.err.println("DEBUG: (findChildScaffolds) found in scafset: ID="+foundscaf.getID()+" "+foundscaf.getCansmi()+" addChild="+ok);
            ++n_cscafs;
            n_cscafs+=foundscaf.getAllChildCount();
          }
          // Check: is scaf already in scafstore?  If yes save time.
          // Surprisingly this does not save time!?
          else if (scafstore!=null && scafstore.scaffoldByCanSmi.contains(cscaf.getCansmi()))
          {
            ScaffoldEntity scent = scafstore.scaffoldByCanSmi.get(cscaf.getCansmi());
            cscaf.setID(scent.getId());
            boolean ok=scaf.addChild(cscaf);
            ++n_cscafs;
            cscaf.setParentScaffold(scaf);
            scafstore.populateScaffoldTree(cscaf);
            n_cscafs+=cscaf.getAllChildCount();
          }
          else if (scafdb!=null && scafdb.containsScaffoldByCansmi(cscaf.getCansmi()))
          {
            ScaffoldRecord scafrec=scafdb.getScaffoldByCansmi(cscaf.getCansmi());
            cscaf.setID(scafrec.getID());
            boolean ok=scaf.addChild(cscaf);
            ++n_cscafs;
            cscaf.setParentScaffold(scaf);
            scafdb.populateScaffoldTree(cscaf);
            n_cscafs+=cscaf.getAllChildCount();
          }
          else
          {
            boolean ok=scaf.addChild(cscaf);
            ++n_cscafs;
            cscaf.setParentScaffold(scaf);
            n_cscafs+=findChildScaffolds(cscaf,scafset,scafstore,scafdb); //recurse; size=0 if leaf
          }
        }
      }
    }
    //System.err.println("DEBUG: (findChildScaffolds) leaving; n_cscafs: "+n_cscafs);
    return n_cscafs;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Re-create Scaffold and all offspring as stored for use as sub-tree,
        etc.
  */
  public static Scaffold entity2Scaffold(ScaffoldStore scafstore,ScaffoldEntity scent)
        throws MolFormatException,DatabaseException
  {
    Scaffold scaf = new Scaffold(scent.getCanSmi(),scafstore.isStereo(),scafstore.isKeep_nitro_attachments());
    scaf.setID(scent.getId());
    for (long id: scent.getChildIds())
    {
      ScaffoldEntity cscent = scafstore.scaffoldById.get(id);
      Scaffold cscaf = entity2Scaffold(scafstore,cscent);       //recurse
      boolean ok=scaf.addChild(cscaf);
      cscaf.setParentScaffold(scaf);
    }
    return scaf;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Tag junction bonds joining scafs, linkers and sidechains, using 
	MolBond.setSetseq().  Set to unique connection ID &gt;1.  
	Not tagged are junctions connecting linkers and sidechains.
	@param mol molecule to be tagged
	@param keep_nitro_attachments true for N-attachments scaf definition
	@return all junction bond indexes
  */
  public static ArrayList<Integer> tagJunctions(Molecule mol,
	boolean keep_nitro_attachments)
    throws MolExportException,SearchException,MolFormatException
  {
    for (MolBond bond: mol.getBondArray()) bond.setSetSeq(0); //clear all

    // (first two atoms in smarts are junction)
    ArrayList<String>  junc_smarts = new ArrayList<String>();
    ArrayList<MolSearch>  junc_pats = new ArrayList<MolSearch>();
    junc_smarts.add("[!R](-!@[!R]=,#[R])=,#[R]");       //scaf-scaf
    junc_smarts.add("[!R](-!@[*])=,#[R]");      //scaf-nonscaf
    if (keep_nitro_attachments)
    {
      junc_smarts.add("[R&!#7]~!@[R&!#7]");     //scaf-scaf
      junc_smarts.add("[R&!#7]~!@[!R]=,#[R]");  //scaf-scaf
      junc_smarts.add("[R&!#7]-!@[!R;!$([R&#7])]");     //scaf-nonscaf
      junc_smarts.add("[!R]~!@[!R]-!@[R&#7]");  //scaf-nonscaf
    }
    else
    {
      junc_smarts.add("[R]-!@[R]");     //scaf-scaf
      junc_smarts.add("[R]-!@[!R]=,#[R]");      //scaf-scaf
      junc_smarts.add("[!R]-!@[R]");    //scaf-nonscaf
    }
    for (String smarts: junc_smarts)
    {
      MolSearch junc_pat = new MolSearch();
      junc_pat.setQuery(MolImporter.importMol(smarts,"smarts:"));
      junc_pat.setTarget(mol);
      junc_pats.add(junc_pat);
    }

    int con_id=0;       // junction idx
    ArrayList<Integer> jbonds = new ArrayList<Integer>();
    for (MolSearch junc_pat: junc_pats)
    {
      int[][] matches=junc_pat.findAll();
      if (matches==null) continue;
      for (int i=0;i<matches.length;++i)
      {
        if (matches[i].length<2) { continue; } //ERROR!
        MolAtom a1=mol.getAtom(matches[i][0]);
        MolAtom a2=mol.getAtom(matches[i][1]);
        MolBond b=a1.getBondTo(a2);
        if (jbonds.contains(mol.indexOf(b))) continue; //already-tagged
        jbonds.add(mol.indexOf(b));
        ++con_id;
        b.setSetSeq(con_id); //tag bond
      }
    }
    return jbonds;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Remove side chains leaving only scaffold.
	Removes, iteratively, terminal atom by terminal atom.
	Scaffolds and linkers should include atoms multiply-bonded to the scaffold
	or linker.  
	Arguably the algorithm used here involving smarts is a kludge, and
	should simply rely on the junction bonds and atoms to figure out
	which atoms to remove.  However, I currently do not see how
	to do that, since, for example, a given atom might need removal
	if part of a "dangling" linker, but not if the linker is involved
	in the scaffold.  Maybe the tagJunctions() method should also
	annotate all atoms as scaffold, linker or sidechain.
	@param mol input molecule or scaffold
	@param keep_nitro_attachments true for N-attachments scaf definition
	@return number of atoms removed
  */
  public static int rmSideChains(Molecule mol,boolean keep_nitro_attachments)
    throws MolFormatException
  {
    int n_del=0;
    ArrayList<MolSearch> pats = new ArrayList<MolSearch>();
    if (keep_nitro_attachments)
    {
      MolSearch pat;
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]-[!#1;!$([#7&R])]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D2;!R;$(*-[!$([#7&R])])]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D2;!R]~[D2;!R;$(*-[!$([#7&R])])]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D3;!R](~[D1])-[!$([R&#7])]","smarts:"));
      pats.add(pat);
    }
    else
    {
      MolSearch pat;
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]-[!#1]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D2;!R;$(*-*)]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D2;!R]~[D2;!R;$(*-*)]","smarts:"));
      pats.add(pat);
      pat = new MolSearch();
      pat.setQuery(MolImporter.importMol("[D1&!#1]~[D3;!R](~[D1])-[*]","smarts:"));
      pats.add(pat);
    }
    while (true)
    {
      int n_del_this=0;
      for (MolSearch pat: pats)
      {
        pat.setTarget(mol); //must re-initialize after mol change
        int[] match;
        try { match=pat.findFirst(); }
        catch (SearchException e) {
          System.err.println(e.getMessage());
          continue;
        }
        catch (Exception e) {
          System.err.println(e.getMessage()); //ArrayIndexOutOfBoundsException: -1
          break;
        }
        if (match==null) continue;
        MolAtom atom=mol.getAtom(match[0]);
        MolAtom nbr=mol.getAtom(match[1]);
        MolBond bond=atom.getBondTo(nbr);
        if (bond.getSetSeq()>0) //may signify junction
        {
          MolAtom h_new = new MolAtom(1);
          mol.add(h_new);
          MolBond b_new = new MolBond(h_new,nbr);
          b_new.setSetSeq(bond.getSetSeq());
          mol.add(b_new);
        }
        else
        {
          nbr.setImplicitHcount(1+nbr.getImplicitHcount());
        }
        mol.removeAtom(atom,MoleculeGraph.RMCLEANUP_EDGES);
        ++n_del_this;
      }
      if (n_del_this==0) break;
      n_del+=n_del_this;
    }
    rmPart(mol,"[H+]"); //clean up explicit-Hs
    rmPart(mol,"[H]"); //clean up explicit-Hs
    return n_del;
  }
  ///////////////////////////////////////////////////////////////////////////
  /** removes all fragments which match the specified smiles
	@param mol input molecule consisting of fragments
	@param usmi canonical smiles of frags to be removed
	@return number of fragments removed
  */
  public static int rmPart(Molecule mol,String usmi)
  {
    int n_del=0;
    Molecule[] partmols=mol.convertToFrags();
    for (Molecule partmol: partmols)
    {
      String usmi_this=null;
      try { usmi_this=MolExporter.exportToFormat(partmol,smifmt); }
      catch (MolExportException e) { usmi_this=""; }
      catch (IOException e) { usmi_this=""; }
      if (usmi_this.equals(usmi)) { ++n_del; }
      else { mol.fuse(partmol,true); }
    }
    return n_del;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Replaces junction hydrogens with junction pseudo atoms (labeled "J") .
	@param mol input molecule
	@return number of junction hydrogens replaced
  */
  public static int replaceJHydrogensWithJPseudoatoms(Molecule mol)
  {
    int n_j=0;
    for (MolBond bond:mol.getBondArray())
    {
      if (bond.getSetSeq()>0)
      {
        for (MolAtom atom: Arrays.asList(bond.getAtom1(),bond.getAtom2()))
        {
          if (atom.getAtno()==1)
          {
            atom.setAtno(MolAtom.PSEUDO);
            atom.setAliasstr("J");
            ++n_j;
          }
        }
      }
    }
    return n_j;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Replaces junction pseudo atoms (labeled "J") with junction hydrogens .
	@param mol input molecule
	@return number of junction pseudo atoms replaced
  */
  public static int replaceJPseudoatomsWithJHydrogens(Molecule mol)
  {
    int n_j=0;
    for (MolBond bond:mol.getBondArray())
    {
      if (bond.getSetSeq()>0)
      {
        for (MolAtom atom: Arrays.asList(bond.getAtom1(),bond.getAtom2()))
        {
          if (atom.getAtno()==MolAtom.PSEUDO)
          {
            atom.setAtno(1);
            ++n_j;
          }
        }
      }
    }
    return n_j;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	*****************EXPERIMENTAL*****************
	Called after tagJunctions, tags every atom identifying it as
	scaffold, linker, or sidechain.  Note that sidechains may be
	joined to scaffolds or linkers.
	<br />
	Contains experimental code in development not
	yet utilized, eventually to improve efficiency and reduce 
	reliance on smarts patterns.
	@param mol molecule to be tagged
	@param jbonds junction bond list from tagJunctions()
  */
  public static void tagAllRoles(Molecule mol,ArrayList<Integer> jbonds)
    throws MolExportException,SearchException,MolFormatException
  {
    // 1. Tag all ring atoms "s" for scaffold.
    TopologyAnalyser topan = new TopologyAnalyser();
    topan.setMolecule(mol);
//    for (MolAtom atom: mol.getAtomArray())
//    {
//      if (topan.isRingAtom(mol.indexOf(atom))) //OK in JChem 5.4, deprecated in 5.5
//      //if (Ring.isRingAtom(mol.indexOf(atom)))  //In JChem 5.5 should work, but where is Ring class?
//      {
//        atom.setAliasstr("s"); //s = scaffold
//      }
//    }
    // 2. For each atom scaffold tagged, depth first graph search terminating
    //    at junction bonds and tag all atoms as scaffold.  Also  terminate at
    //    atoms already tagged.
    for (MolAtom atom: mol.getAtomArray())
    {
      if (atom.getAliasstr().equals("s"))
        tagRole(atom,"s");
    }
    // 3. For all pairwise combinations of junction bonds, find the shortest
    //    path.  If that path does not include any scaffold atoms,
    //    along that path tag atoms as linker ("l");
    for (int bidxA: jbonds)
    {
      MolBond bondA = mol.getBond(bidxA);
      MolAtom atomA=null;
      for (MolAtom atom: Arrays.asList(bondA.getAtom1(),bondA.getAtom2()))
      {
        if (!atom.getAliasstr().equals("s")) { atomA=atom; }
      }
      if (atomA==null) continue;
      for (int bidxB: jbonds)
      {
        MolBond bondB = mol.getBond(bidxB);
        MolAtom atomB=null;
        for (MolAtom atom: Arrays.asList(bondB.getAtom1(),bondB.getAtom2()))
        {
          if (!atom.getAliasstr().equals("s")) { atomB=atom; }
        }
        if (atomB==null) continue;
        ShortestPath sp = new ShortestPath();
        sp.calculate(mol);

        //jj: what are the path atoms?
        int [] aidxs = new int [mol.getAtomCount()];
        int len = sp.getPath(mol.indexOf(atomA),mol.indexOf(atomB),aidxs);
      }
    }
  }
  /**	From specified atom, depth first graph search terminating
	at junction bonds and tag all atoms.  Also  terminate at
	atoms already tagged.
	<br />
	Contains experimental code in development not
	yet utilized.
  */
  private static int tagRole(MolAtom atom,String tag)
  {
    int n_tag=0;
    for (int i=0;i<atom.getBondCount();++i)
    {
      MolBond bond = atom.getBond(i);
      MolAtom a2 = bond.getOtherAtom(atom);
      if (bond.getSetSeq()==0 && !a2.getAliasstr().equals(tag))
      {
        a2.setAliasstr(tag);
        ++n_tag;
        n_tag+=tagRole(atom,tag);
      }
    }
    return n_tag;
  }

  public static String cansmi(String smi,boolean stereo)
    throws IOException,MolFormatException
  {
    Molecule mol=MolImporter.importMol(smi,"smiles:");
    String cansmi=MolExporter.exportToFormat(mol,(stereo?Scaffold.CANSMIFMT_STEREO:Scaffold.CANSMIFMT));
    return cansmi;
  }
}
