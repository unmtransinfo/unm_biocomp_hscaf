package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;
 
import edu.unm.health.biocomp.hscaf.*;

/////////////////////////////////////////////////////////////////////////////
/**	Contains static functions used in hierarchical scaffold analysis&#46;
	These public API functions for use in applications&#46;
	@see edu.unm.health.biocomp.hscaf.ScaffoldTree
	@see edu.unm.health.biocomp.hscaf.Scaffold
	@see edu.unm.health.biocomp.hscaf.Linker
	@see edu.unm.health.biocomp.hscaf.Sidechain
	@see edu.unm.health.biocomp.hscaf.hier_scaffolds
	@author Jeremy J Yang
*/
public class hier_scaffolds_utils
{
  private hier_scaffolds_utils() {} //disallow default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	Finds maximum common scaffold in two ScaffoldTrees&#46;
	If none returns null&#46;
  */
  public static Scaffold maxCommonScaffold(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB)
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
  ///////////////////////////////////////////////////////////////////////////
  /**	Returns largest fragment of fragmented molecule&#46;
  */
  public static Molecule largestPart(Molecule mol)
  {
    int acount_largest=0;
    int i_largest=0;
    Molecule[] partmols=mol.convertToFrags();
    for (int i=0;i<partmols.length;++i)
    {
      if (partmols[i].getAtomCount()>acount_largest) i_largest=i;
    }
    return partmols[i_largest];
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Returns similarity metric designed to reflect the significance of
	the maximum common scaffold between two molecules&#46;  Given by the
	formula:
	<br><tt>
	&nbsp; <b>sim = nc / (na + nb - nc) </b><br>
	where <br>
	&nbsp; nc = # atoms in common scaffold <br>
	&nbsp; na = # atoms in max scaffold, mol A <br>
	&nbsp; nb = # atoms in max scaffold, mol B <br>
	</tt>
	This is an experimental functionality&#46;  Intended for use in combination
	with a standard whole-molecule similarity metric&#46;
  */
  public static float commonScaffoldTanimoto(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB)
  {
    Scaffold mcScaf=maxCommonScaffold(scaftreeA,scaftreeB);
    if (mcScaf==null) return 0.0f;
    int nc=mcScaf.getAtomCount();
    int nA=scaftreeA.getRootScaffold().getAtomCount();
    int nB=scaftreeB.getRootScaffold().getAtomCount();
    float sim = (float)nc / (nA + nB - nc);
    return sim;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Returns similarity metric designed to reflect the significance of
	the maximum common scaffold between two molecules&#46;  Given by the
	formula:
	<br><tt>
	&nbsp; <b>sim = nc / (na + nb - nc) </b><br>
	where <br>
	&nbsp; nc = # atoms in common scaffold <br>
	&nbsp; na = # atoms in max scaffold, mol A <br>
	&nbsp; nb = # atoms in max scaffold, mol B <br>
	</tt>
	This is an experimental functionality&#46;  Intended for use in combination
	with a standard whole-molecule similarity metric&#46;
	Use this faster method if mcScaf is already calculated (with maxCommonScaffold())&#46;
  */
  public static float commonScaffoldTanimoto(ScaffoldTree scaftreeA,ScaffoldTree scaftreeB,
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
  /**	Returns simple count of all ringsystems&#46;  Equivalent to count of
	all disconnected fragments after removing linkers&#46;  Used to 
	estimate computational demands and identify high-cost
	"pathological" molecules for special handling&#46;
  */
  public static int rawRingsystemCount(Molecule mol)
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
    int count=xmol.getFragCount();
    //int count=xmol.getFragCount(MoleculeGraph.FRAG_BASIC);
    return count;
  }
  ///////////////////////////////////////////////////////////////////////////
}
