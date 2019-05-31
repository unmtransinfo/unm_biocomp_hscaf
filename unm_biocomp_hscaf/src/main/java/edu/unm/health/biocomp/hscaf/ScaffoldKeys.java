package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*; //MolSearch
import chemaxon.marvin.io.MolExportException;
import chemaxon.marvin.calculations.*; //TopologyAnalyserPlugin
import chemaxon.marvin.plugin.*; //PluginException
 
/**	Represents Scaffold Keys descriptors by Peter Ertl.
	"Intuitive Ordering of Scaffolds and Scaffold Similarity Searching Using Scaffold Keys",
	Peter Ertl, J. Chem Inf Mod, 54, 2014, pp.1617-1622.
	<br>
	@author Jeremy J Yang
*/
public class ScaffoldKeys implements Comparable<Object>
{
  public final int N_KEYS = 32;
  private Molecule mol;
  private Molecule[] qmols = new Molecule[N_KEYS]; //parsed smarts
  private Integer[] kvals = new Integer[N_KEYS];
  private MolSearch msrch = null;
  private TopologyAnalyserPlugin topoPlugin  = null;
  
  private final String[] tags = new String[]{
	"RING_AND_LINK_ATOMS",		//1. ring and linker atoms
	"LINK_ATOMS",			//2. linker atoms
	"LINK_BONDS",			//3. linker bonds
	"RINGS",			//4. rings
	"SPIRO_ATOMS",			//5. spiro atoms
	"LARGEST_RING_SIZE",		//6. size of largest ring
	"CONJ_RING_BONDS",		//7. bonds in fully conjugated rings
	"UNCONJ_RING_MBONDS",		//8. multiple bonds in not fully conjugated rings
	"RING_HETEROATOMS",		//9. heteroatoms in rings
	"RING_HETEROATOMS_NOT_NSO",	//10. heteroatoms in rings, not N, S, O
	"S_RING_ATOMS",			//11. S ring atoms
	"O_RING_ATOMS",			//12. O ring atoms
	"N_RING_ATOMS",			//13. N ring atoms
	"HETEROATOMS",			//14. heteroatoms
	"HETEROATOMS_NOT_NSO",		//15. heteroatoms, not N, S, O
	"S_ATOMS",			//16. S atoms
	"O_ATOMS",			//17. O atoms
	"N_ATOMS",			//18. N atoms
	"LINK_MBONDS",			//19. multiple linker bonds
	"ADJACENT_HETEROATOMS2",	//20. 2 adjacent heteroatoms
	"ADJACENT_HETEROATOMS3",	//21. 3 adjacent heteroatoms
	"CSEPARATED_HETEROATOMS",	//22. 2 heteroatoms separated by a single C
	"C2SEPARATED_HETEROATOMS",	//23. 2 heteroatoms separated by a 2 Cs
	"HETEROATOM_DBONDS",		//24. double bonds with at least one heteroatom
	"ADJACENT_DBOND_HETEROATOMS",	//25. heteroatoms adjacent to a double (non-aromatic) bond
	"CONJ_DBOND_PAIRS",		//26. pairs of conjugated double (non-aromatic) bonds
	"ADJACENT_BRANCHEDATOM_PAIRS",	//27. pairs of adjacent branched atoms
	"SEPARATED_BRANCHEDATOMS",	//28. branched atoms separated by a single nonbranched atom
	"ADJACENT_BRANCHEDATOMS3",	//29. 3 adjacent branched atoms
	"SEPARATED_BRANCHEDATOMS2",	//30. branched atoms separated by any 2 atoms
	"EXO_ATOMS",			//31. exocyclic and exolinker atoms
	"GT2BOND_HETEROATOMS"};		//32. heteroatoms with more than 2 bonds

  private final double[] avgs = new double[]{
	20.029, 2.518, 3.993, 3.348, 0.031, 6.241, 13.824, 0.112,
	2.177, 0.003, 0.143, 0.310, 1.721, 4.248, 0.009, 0.289,
	1.603, 2.347, 0.109, 0.575, 0.350, 1.804, 1.505, 1.235,
	1.439, 0.094, 2.860, 1.504, 1.734, 4.294, 1.170, 0.673
	};
  
  private final double[] stds = new double[]{
	7.556, 3.481, 3.897, 3.156, 0.193, 1.905, 6.201, 0.383,
	1.640, 0.061, 0.388, 0.703, 1.507, 2.921, 0.131, 0.540,
	1.695, 1.789, 0.351, 1.162, 1.169, 1.953, 2.564, 1.433,
	1.861, 0.380, 2.320, 1.467, 2.591, 4.409, 1.425, 0.840
  };
 
  /////////////////////////////////////////////////////////////////////////////
  private ScaffoldKeys()
  {
    this.mol=null;
    this.msrch = new MolSearch();
    this.topoPlugin  = new TopologyAnalyserPlugin();
    try { this.initializeQmols(); }
    catch (Exception e) { System.err.println(e.getMessage()); }
  }
  /////////////////////////////////////////////////////////////////////////////
  public void setMolecule(Molecule _mol)
  {
    this.mol=_mol;
    this.msrch.setTarget(this.mol);
    try {
      this.topoPlugin.setMolecule(mol);
      this.topoPlugin.run();
    } catch (Exception e) { System.err.println(e.getMessage()); }
  }
  public Molecule getMolecule() { return this.mol; }
  /////////////////////////////////////////////////////////////////////////////
  public String[] getTags() { return this.tags; }
  /////////////////////////////////////////////////////////////////////////////
  /**	idx start with 0 (as published minus 1).
  */
  public String getTag(int idx)
  {
    return (idx>=this.tags.length)?null:this.tags[idx];
  }
  /////////////////////////////////////////////////////////////////////////////
  public double getVal(int idx)
  {
    double val = (this.kvals[idx] - this.avgs[idx]) / this.stds[idx];
    return val;
  }
  /////////////////////////////////////////////////////////////////////////////
  public int getRawVal(int idx)
  {
    return this.kvals[idx];
  }
  /////////////////////////////////////////////////////////////////////////////
  private void initializeQmols()
	throws MolFormatException
  {
    this.qmols[0]=MolImporter.importMol("[!#1]","smarts:d"); //1."RING_AND_LINK_ATOMS"
    this.qmols[1]=MolImporter.importMol("[!R&!X1]","smarts:d"); //2. "LINK_ATOMS"
    this.qmols[2]=MolImporter.importMol("*~!@*","smarts:d"); //3. "LINK_BONDS"
    this.qmols[3]=null; //4. "RINGS"
    this.qmols[4]=MolImporter.importMol("[R2&X4]","smarts:d"); //5. "SPIRO_ATOMS"
    this.qmols[5]=null; //6. "LARGEST_RING_SIZE"
    this.qmols[6]=MolImporter.importMol("a:a","smarts:d"); //7. "CONJ_RING_BONDS"
    this.qmols[7]=MolImporter.importMol("*=,#*","smarts:d"); //8. "UNCONJ_RING_MBONDS"
    this.qmols[8]=MolImporter.importMol("[!#6&!#1&R]","smarts:d"); //9. "RING_HETEROATOMS"
    this.qmols[9]=MolImporter.importMol("[R&!#6&!#1&!#7&!#16&!#8]","smarts:d"); //10. "RING_HETEROATOMS_NOT_NSO"
    this.qmols[10]=MolImporter.importMol("[#16;R]","smarts:d"); //11. "S_RING_ATOMS"
    this.qmols[11]=MolImporter.importMol("[#8;R]","smarts:d"); //12. "O_RING_ATOMS"
    this.qmols[12]=MolImporter.importMol("[#7;R]","smarts:d"); //13. "N_RING_ATOMS"
    this.qmols[13]=MolImporter.importMol("[!#6&!#1]","smarts:d"); //14. "HETEROATOMS"
    this.qmols[14]=MolImporter.importMol("[!#6&!#1&!#7&!#16&!#8]","smarts:d"); //15. heteroatoms "HETEROATOMS_NOT_NSO"
    this.qmols[15]=MolImporter.importMol("[#16]","smarts:d"); //16. "S_ATOMS"
    this.qmols[16]=MolImporter.importMol("[#8]","smarts:d"); //17. "O_ATOMS"
    this.qmols[17]=MolImporter.importMol("[#7]","smarts:d"); //18. "N_ATOMS"
    this.qmols[18]=MolImporter.importMol("*~!@&!-*","smarts:d"); //19. "LINK_MBONDS"
    this.qmols[19]=MolImporter.importMol("[!#6&!#1]~[!#6&!#1]","smarts:d"); //20. "ADJACENT_HETEROATOMS2"
    this.qmols[20]=MolImporter.importMol("[!#6&!#1]~[!#6&!#1]~[!#6&!#1]","smarts:d"); //21. "ADJACENT_HETEROATOMS3"
    this.qmols[21]=MolImporter.importMol("[!#6&!#1]~[#6]~[!#6&!#1]","smarts:d"); //22. "CSEPARATED_HETEROATOMS"
    this.qmols[22]=MolImporter.importMol("[!#6&!#1]~[#6]~[#6]~[!#6&!#1]","smarts:d"); //23. "C2SEPARATED_HETEROATOMS"
    this.qmols[23]=MolImporter.importMol("[!#6&!#1]=[!#6&!#1]","smarts:d"); //24. "HETEROATOM_DBONDS"
    this.qmols[24]=MolImporter.importMol("[!#6&!#1]=*","smarts:d"); //25. "ADJACENT_DBOND_HETEROATOMS"
    this.qmols[25]=MolImporter.importMol("*=*-*=*","smarts:d"); //26. "CONJ_DBOND_PAIRS"
    this.qmols[26]=MolImporter.importMol("[X3,X4]~[X3,X4]","smarts:d"); //27. "ADJACENT_BRANCHEDATOM_PAIRS"
    this.qmols[27]=MolImporter.importMol("[X3,X4]~[X2]~[X3,X4]","smarts:d"); //28. "SEPARATED_BRANCHEDATOMS"
    this.qmols[28]=MolImporter.importMol("[X3,X4]~[X3,X4]~[X3,X4]","smarts:d"); //29. "ADJACENT_BRANCHEDATOMS3"
    this.qmols[29]=MolImporter.importMol("[X3,X4]~*~*~[X3,X4]","smarts:d"); //30. "SEPARATED_BRANCHEDATOMS2"
    this.qmols[30]=MolImporter.importMol("[X1]","smarts:d"); //31. "EXO_ATOMS"
    this.qmols[31]=MolImporter.importMol("[!#6;!#1;!X0;!X1;!X2]","smarts:d"); //32. "GT2BOND_HETEROATOMS"
  }
  /////////////////////////////////////////////////////////////////////////////
  // Must deal with Unique_Set-of_Atom (USA) counts.
  // Current kludge: divide by 2 for symmetric patterns
  // "conjugated" = "aromatic" ?
  /////////////////////////////////////////////////////////////////////////////
  public boolean computeKeys()
	throws SearchException,PluginException
  {
    this.msrch.setQuery(qmols[0]); this.kvals[0]=this.msrch.getMatchCount(); //1. "RING_AND_LINK_ATOMS"
    this.msrch.setQuery(qmols[1]); this.kvals[1]=this.msrch.getMatchCount(); //2. "LINK_ATOMS"
    this.msrch.setQuery(qmols[2]); this.kvals[2]=this.msrch.getMatchCount()/2; //3. "LINK_BONDS"
    this.kvals[3]=this.topoPlugin.getRingCount(); //4. "RINGS"
    this.msrch.setQuery(qmols[4]); this.kvals[4]=this.msrch.getMatchCount(); //5. "SPIRO_ATOMS"
    this.kvals[5]=this.topoPlugin.getLargestRingSize(); //6. "LARGEST_RING_SIZE"
    this.msrch.setQuery(qmols[6]); this.kvals[6]=this.msrch.getMatchCount()/2; //7. "CONJ_RING_BONDS"
    this.msrch.setQuery(qmols[7]); this.kvals[7]=this.msrch.getMatchCount()/2; //8. "UNCONJ_RING_MBONDS"
    this.msrch.setQuery(qmols[8]); this.kvals[8]=this.msrch.getMatchCount(); //9. "RING_HETEROATOMS"
    this.msrch.setQuery(qmols[9]); this.kvals[9]=this.msrch.getMatchCount(); //10. "RING_HETEROATOMS_NOT_NSO"
    this.msrch.setQuery(qmols[10]); this.kvals[10]=this.msrch.getMatchCount(); //11. "S_RING_ATOMS"
    this.msrch.setQuery(qmols[11]); this.kvals[11]=this.msrch.getMatchCount(); //12. "O_RING_ATOMS"
    this.msrch.setQuery(qmols[12]); this.kvals[12]=this.msrch.getMatchCount(); //13. "N_RING_ATOMS"
    this.msrch.setQuery(qmols[13]); this.kvals[13]=this.msrch.getMatchCount(); //14. "HETEROATOMS"
    this.msrch.setQuery(qmols[14]); this.kvals[14]=this.msrch.getMatchCount(); //15. "HETEROATOMS_NOT_NSO"
    this.msrch.setQuery(qmols[15]); this.kvals[15]=this.msrch.getMatchCount(); //16. "S_ATOMS"
    this.msrch.setQuery(qmols[16]); this.kvals[16]=this.msrch.getMatchCount(); //17. "O_ATOMS"
    this.msrch.setQuery(qmols[17]); this.kvals[17]=this.msrch.getMatchCount(); //18. "N_ATOMS"
    this.msrch.setQuery(qmols[18]); this.kvals[18]=this.msrch.getMatchCount()/2; //19. "LINK_MBONDS"
    this.msrch.setQuery(qmols[19]); this.kvals[19]=this.msrch.getMatchCount()/2; //20. "ADJACENT_HETEROATOMS2"
    this.msrch.setQuery(qmols[20]); this.kvals[20]=this.msrch.getMatchCount()/2; //21. "ADJACENT_HETEROATOMS3"
    this.msrch.setQuery(qmols[21]); this.kvals[21]=this.msrch.getMatchCount()/2; //22. "CSEPARATED_HETEROATOMS"
    this.msrch.setQuery(qmols[22]); this.kvals[22]=this.msrch.getMatchCount()/2; //23. "C2SEPARATED_HETEROATOMS"
    this.msrch.setQuery(qmols[23]); this.kvals[23]=this.msrch.getMatchCount()/2; //24. "HETEROATOM_DBONDS"
    this.msrch.setQuery(qmols[24]); this.kvals[24]=this.msrch.getMatchCount(); //25. "ADJACENT_DBOND_HETEROATOMS"
    this.msrch.setQuery(qmols[25]); this.kvals[25]=this.msrch.getMatchCount()/2; //26. "CONJ_DBOND_PAIRS"
    this.msrch.setQuery(qmols[26]); this.kvals[26]=this.msrch.getMatchCount()/2; //27. "ADJACENT_BRANCHEDATOM_PAIRS"
    this.msrch.setQuery(qmols[27]); this.kvals[27]=this.msrch.getMatchCount()/2; //28. "SEPARATED_BRANCHEDATOMS"
    this.msrch.setQuery(qmols[28]); this.kvals[28]=this.msrch.getMatchCount()/2; //29. "ADJACENT_BRANCHEDATOMS3"
    this.msrch.setQuery(qmols[29]); this.kvals[29]=this.msrch.getMatchCount()/2; //30. "SEPARATED_BRANCHEDATOMS2"
    this.msrch.setQuery(qmols[30]); this.kvals[30]=this.msrch.getMatchCount()/2; //31. "EXO_ATOMS"
    this.msrch.setQuery(qmols[31]); this.kvals[31]=this.msrch.getMatchCount()/2; //32. "GT2BOND_HETEROATOMS"
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Similarity/distance.
  */
  public double distanceTo(ScaffoldKeys sk2)
  {
    double d=0.0;
    for (int i=0;i<this.N_KEYS;++i)
      d+=(Math.pow(Math.abs(this.getVal(i) - sk2.getVal(i)),1.5) / (i+1));
    return d;
  }

  /////////////////////////////////////////////////////////////////////////////
  /**   For sorting.
  */
  public int compareTo(Object o) throws ClassCastException
  {
    int d=0;
    for (int i=0;i<N_KEYS;++i)
    {
      d = ((ScaffoldKeys)this).getRawVal(i)-((ScaffoldKeys)o).getRawVal(i);
      if (d!=0) break;
    }
    return d;
  }
  /////////////////////////////////////////////////////////////////////////////
  public boolean equals(ScaffoldKeys sks2)
  {
    return (this.compareTo(sks2)==0);
  }

  /////////////////////////////////////////////////////////////////////////////
  /**   For testing only.
  */
  public static void main(String[] args)
	throws Exception
  {
    if (args.length<2)
    {
      System.err.println("ERROR: syntax: ScaffoldKeys SCAFSMI_A SCAFSMI_B");
      System.exit(1);
    }

    Scaffold scafA = new Scaffold(args[0],false,false);
    System.err.println("A: "+scafA.getCansmi());
    Scaffold scafB = new Scaffold(args[1],false,false);
    System.err.println("B: "+scafB.getCansmi());

    ScaffoldKeys sksA = new ScaffoldKeys();
    ScaffoldKeys sksB = new ScaffoldKeys();
    sksA.setMolecule(scafA);
    sksB.setMolecule(scafB);
    try {
      sksA.computeKeys();
      sksB.computeKeys();
    } catch (Exception e) {
      System.err.println(e.getMessage());
      System.exit(1);
    }
    for (int i=0;i<sksA.N_KEYS;++i)
    {
      System.err.println(String.format("\t%2d. %28s: A: %8.3f ; B: %8.3f",i+1,sksA.getTag(i),sksA.getVal(i),sksB.getVal(i)));
    }
    System.err.println(String.format("\tA-B distance = %8.3f",sksA.distanceTo(sksB)));
    System.err.println(String.format("\tRanking: A.compareTo(B) = %d",sksA.compareTo(sksB)));
  }
}
