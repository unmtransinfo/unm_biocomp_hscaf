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
  /**	Assigns scaffold IDs to all scaffolds in tree&#46;  If scaffold is new,
	also new ID and updates scaf_usmis_global
  */
  public static void assignScaffoldIDs(HashMap<String,Integer> scaf_usmis_global,
	ScaffoldTree scaftree)
  {
    for (Scaffold scaf: scaftree.getAllScaffolds())
    if (scaf_usmis_global.containsKey(scaf.getCansmi()))
    {
      scaf.setID(scaf_usmis_global.get(scaf.getCansmi())); //existing scafID
    }
    else
    {
      scaf.setID(scaf_usmis_global.size()+1); //new scafID
      scaf_usmis_global.put(scaf.getCansmi(),scaf.getID());
    }
    return;
  }
  ///////////////////////////////////////////////////////////////////////////
  /**	Generates string representing the hierarchical scaffold tree&#46;
	e&#46;g&#46; "1:(2,3)" or "1:(2:(3,4,5),6:(4,7))"
  */
  public static String scafTreeAsString(Scaffold scaf)
  {
    String str="";
    if (scaf==null) return str;
    str+=scaf.getID();
    //if (scaf.getID()==0) System.err.println("DEBUG(scafTreeAsString): scafid==0 for scaf: "+scaf.getCansmi());
    if (scaf.getChildCount()>0)
    {
      str+=":(";
      int i=0;
      for (Scaffold cscaf: scaf.getChildScaffolds())
      {
        if (i>0) str+=",";
        str+=scafTreeAsString(cscaf);
        ++i;
      }
      str+=")";
    }
    return str;
  }
  ///////////////////////////////////////////////////////////////////////////
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
}
