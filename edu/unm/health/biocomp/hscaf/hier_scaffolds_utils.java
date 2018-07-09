package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.*;
import chemaxon.sss.search.*;
import chemaxon.util.MolHandler;
import chemaxon.marvin.io.MolExportException;

import com.sleepycat.je.DatabaseException;
 
import edu.unm.health.biocomp.hscaf.*;
import edu.unm.health.biocomp.db.*;

/////////////////////////////////////////////////////////////////////////////
/**	Contains static functions used in hierarchical scaffold analysis.
	These public API functions for use in applications.
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@author Jeremy J Yang
*/
public class hier_scaffolds_utils
{
  private hier_scaffolds_utils() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	Finds maximum common scaffold in two ScaffoldTrees.
	If none returns null.
  */
  public static Scaffold MaxCommonScaffold(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB)
  {
    ArrayList<Scaffold> scafsA = scaftreeA.getScaffolds();
    ArrayList<Scaffold> scafsB = scaftreeB.getScaffolds();
    Collections.sort(scafsA,Collections.reverseOrder()); //sort by decreasing size
    Collections.sort(scafsB,Collections.reverseOrder()); //sort by decreasing size
    int iA=0;
    int iB=0;
    while (iA<scafsA.size() && iB<scafsB.size())
    {
      Scaffold scafA = scafsA.get(iA);
      Scaffold scafB = scafsB.get(iB);
      if (scafA.equals(scafB))
        return new Scaffold(scafA);
      else if (scafA.compareTo(scafB)<0)
      {
        ++iA;
      }
      else if (scafA.compareTo(scafB)>0)
      {
        ++iB;
      }
      else
      {
        ++iA;
        ++iB;
      }
    }
    return null;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Finds maximum common scaffold in two molecules.
	If none returns null.
  */
  public static Scaffold MaxCommonScaffold(Molecule molA,Molecule molB,
	boolean stereo,boolean keep_nitro_attachments,ScaffoldSet scafset)
	throws IOException,SearchException,MolFormatException,MolExportException,ScaffoldException,SQLException,DatabaseException
  {
    ScaffoldTree scaftreeA = new ScaffoldTree(molA,stereo,keep_nitro_attachments,scafset);
    ScaffoldTree scaftreeB = new ScaffoldTree(molB,stereo,keep_nitro_attachments,scafset);
    return MaxCommonScaffold(scaftreeA,scaftreeB);
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	@return	largest fragment of fragmented molecule
  */
  public static Molecule LargestPart(Molecule mol)
  {
    Molecule[] partmols=mol.convertToFrags();
    int i_largest=0;
    for (int i=0;i<partmols.length;++i)
    {
      if (i!=i_largest && partmols[i].getAtomCount()>partmols[i_largest].getAtomCount())
        i_largest=i;
    }
    partmols[i_largest].setName(mol.getName());
    return partmols[i_largest];
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Similarity metric designed to reflect the significance of
	the maximum common scaffold between two molecules.  Given by the
	formula:
	<br><tt>
	&nbsp; <b>sim = nc / (na + nb - nc) </b><br>
	where <br>
	&nbsp; nc = # atoms in common scaffold <br>
	&nbsp; na = # atoms in max scaffold, mol A <br>
	&nbsp; nb = # atoms in max scaffold, mol B <br>
	</tt>
	This is an experimental functionality.  For use in combination
	with a standard whole-molecule similarity metric.
  	@return	similarity metric
  */
  public static float CommonScaffoldTanimoto(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB)
  {
    Scaffold mcScaf=MaxCommonScaffold(scaftreeA,scaftreeB);
    if (mcScaf==null) return 0.0f;
    int nc=mcScaf.getAtomCount();
    int nA=scaftreeA.getRootScaffold().getAtomCount();
    int nB=scaftreeB.getRootScaffold().getAtomCount();
    float sim = (float)nc / (nA + nB - nc);
    return sim;
  }
  /**	
  */
  public static float CommonScaffoldTanimoto(Molecule molA,Molecule molB,
	boolean stereo,boolean keep_nitro_attachments,ScaffoldSet scafset)
	throws IOException,SearchException,MolFormatException,MolExportException,ScaffoldException,SQLException,DatabaseException
  {
    ScaffoldTree scaftreeA = new ScaffoldTree(molA,stereo,keep_nitro_attachments,scafset);
    ScaffoldTree scaftreeB = new ScaffoldTree(molB,stereo,keep_nitro_attachments,scafset);
    return CommonScaffoldTanimoto(scaftreeA,scaftreeB);
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Similarity metric designed to reflect the significance of
	the maximum common scaffold between two molecules.  Given by the
	formula:
	<br><tt>
	&nbsp; <b>sim = nc / (na + nb - nc) </b><br>
	where <br>
	&nbsp; nc = # atoms in common scaffold <br>
	&nbsp; na = # atoms in max scaffold, mol A <br>
	&nbsp; nb = # atoms in max scaffold, mol B <br>
	</tt>
	This is an experimental functionality.  For use in combination
	with a standard whole-molecule similarity metric.
	Use this faster method if mcScaf is already calculated (with MaxCommonScaffold()).
  	@return	similarity metric
  */
  public static float CommonScaffoldTanimoto(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB,
    Scaffold mcScaf)
  {
    if (mcScaf==null) return 0.0f;
    int nc=mcScaf.getAtomCount();
    int nA=scaftreeA.getRootScaffold().getAtomCount();
    int nB=scaftreeB.getRootScaffold().getAtomCount();
    float sim = (float)nc / (nA + nB - nc);
    return sim;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Returns simple count of all ringsystems.  Equivalent to count of
	all disconnected fragments after removing linkers.  Used to 
	estimate computational demands and identify high-cost
	"pathological" molecules for special handling.
  	@return simple count of all ringsystems
  */
  public static int RawRingsystemCount(Molecule mol)
  {
    Molecule xmol = mol.cloneMolecule();
    int [][] sssratoms=xmol.getSSSR();
    HashMap<Integer,Boolean> ringatoms = new HashMap<Integer,Boolean>();
    for (int [] sssratom: sssratoms) { for (int ratomidx: sssratom) { ringatoms.put(ratomidx,true); } }
    ArrayList<MolAtom> nonringatoms = new ArrayList<MolAtom>();
    for (MolAtom atom: xmol.getAtomArray())
    {
      if (!(ringatoms.containsKey(xmol.indexOf(atom)))) nonringatoms.add(atom);
    }
    for (MolAtom atom: nonringatoms) { xmol.removeAtom(atom); }
    //Still need to remove non-ring single bonds linking rings.
    int [][] sssrbonds=xmol.getSSSRBonds();
    HashMap<Integer,Boolean> ringbonds = new HashMap<Integer,Boolean>();
    for (int [] sssrbond: sssrbonds) { for (int rbondidx: sssrbond) { ringbonds.put(rbondidx,true); } }
    ArrayList<MolBond> nonringbonds = new ArrayList<MolBond>();
    for (MolBond bond: xmol.getBondArray())
    {
      if (!(ringbonds.containsKey(xmol.indexOf(bond)))) nonringbonds.add(bond);
    }
    for (MolBond bond: nonringbonds) { xmol.removeBond(bond); }
    int count=xmol.getFragCount(MoleculeGraph.FRAG_BASIC);
    return count;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	For convenience depicting molecules with highlighted atoms
	corresponding to a scaffold.
	@param	mol	molecule containing scaffold
	@param	qmol	molecule representing scaffold
  */
  public static boolean SubMatchMap(Molecule mol,Molecule qmol)
	throws MolExportException,SearchException,MolFormatException
  {
    if (mol==null || mol.getAtomCount()==0) return false;
    if (qmol==null || qmol.getAtomCount()==0) return false;

    boolean is_match=false;
    MolSearch search=new MolSearch();
    search.setTarget(mol);
    // Remove exp-H's from query molecule, ignore imp-H's for permissive match.
    qmol.implicitizeHydrogens(MolAtom.ALL_H);
    MolSearchOptions searchOpts = new MolSearchOptions(SearchConstants.SUBSTRUCTURE); 
    searchOpts.setImplicitHMatching(SearchConstants.IMPLICIT_H_MATCHING_DISABLED);
    search.setSearchOptions(searchOpts);
    search.setQuery(qmol);
    boolean ok=false;
    int[][] matchs=null;
    ok=search.isMatching();
    matchs=search.findAll(); //may throw SearchException

    if (!ok || matchs==null) return false;

    //Each match is an array of atom indices.
    for (int[] match: matchs)
    {
      for (int ia=0; ia<match.length; ++ia)
      {
        if (match[ia]>=0)
        {
          mol.getAtom(match[ia]).setAtomMap(1);
          is_match=true;
        }
        else
          mol.getAtom(match[ia]).setAtomMap(0);
      }
    }
    return is_match;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	For convenience depicting molecules with highlighted atoms
	corresponding to a scaffold.
	@param	mol	molecule containing scaffold
	@param	qmol	molecule representing scaffold
  */
  public static String SubMatchMapSmiles(Molecule mol,Molecule qmol)
	throws MolExportException,SearchException,MolFormatException,IOException
  {
    if (mol==null || qmol==null) return "";
    if (mol.getAtomCount()==0) return "";

    boolean is_match=SubMatchMap(mol,qmol);
    String smifmt="smiles:u-H-a";
    String smi=MolExporter.exportToFormat(mol,smifmt);
    return smi;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	For convenience depicting molecules with highlighted atoms
	corresponding to a scaffold.
	@param	smi	smiles containing scaffold
	@param	qsmi	smiles representing scaffold
  */
  public static String SubMatchMapSmiles(String smi,String qsmi)
	throws MolExportException,SearchException,MolFormatException,IOException
  {
    Molecule mol=MolImporter.importMol(smi,"smiles:");
    Molecule qmol=MolImporter.importMol(qsmi,"smiles:");
    return SubMatchMapSmiles(mol,qmol);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Minimal sanity check whether DB exists ("scaffold" and "scaf2scaf"
	tables) to support {@link ScaffoldDB} object.
  */
  public static boolean CheckExistsDB(
	String dbhost,
	int dbport,
	String dbname,
	String dbschema,
	String dbusr,
	String dbpw,
	String dbtableprefix)
        throws SQLException
  {
    Connection dbcon=pg_utils.DBConnect(dbhost,dbport,dbname,dbusr,dbpw);
    boolean ok=true;
    if (dbcon==null) return false;
    String sql="SELECT DISTINCT table_name FROM information_schema.tables WHERE table_schema='"+dbschema+"' AND table_name IN ('"+dbtableprefix+"scaffold','"+dbtableprefix+"scaf2scaf')";
    ResultSet rset=pg_utils.ExecuteSql(dbcon,sql);
    long rowcount=0L;
    while (rset.next()) ++rowcount;
    rset.getStatement().close();
    dbcon.close();
    if (rowcount==2) return true;
    else return false;
  }
}
