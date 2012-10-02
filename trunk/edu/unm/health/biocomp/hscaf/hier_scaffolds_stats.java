package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.text.DateFormat;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.*;

//// to do: [ ] -in_link, -in_chain, -out_link, -out_chain (but first must fix hier_scaffolds.java).

/**	HierS hierarchical scaffolds postprocess application&#46;  Run program with
	no args for command-line help&#46;  Key task is to generate scaffolds file
	annotated with occurance info including compound IDs&#46;
	<br />
	Molecule files normally should be SMILES format&#46;  The SMILES "name" is a space delimited
	set of fields: (1) compound ID (integer), (2) set of scaffold IDs&#46;  For example:
	FC1=CC=C(NC2=NC(=NC3=CC=CC=C23)C2=CN=CC=C2)C=C1 1877191 S:4,5,6,17,18
	<br />
	Input scaffold file also normally SMILES, with "name" delimited set of fields:
	(1) scaffold ID, and (2) scaffold tree string&#46;  For example:
	N(C1=CC=NC=C1)C1=C2C=CC=CC2=NC(=N1)C1=COC=C1    35 35:(38:(5,6),6,36:(5,37),37)
	<br />
	@author Jeremy J Yang
*/
public class hier_scaffolds_stats
{
  private static void help(String msg)
  {
    System.err.println(msg+"\n"
      +"hier_scaffolds_stats - postprocess HierS scaffolds output (from hier_scaffolds)\n"
      +"\n"
      +"From the University of New Mexico Division of Biocomputing\n"
      +"(http://biocomp.health.unm.edu/).\n"
      +"\n"
      +"usage: hier_scaffolds_stats [options]\n"
      +"  required:\n"
      +"    -i IFILE ................... input file of molecules annotated with scaf info\n"
      +"    -in_scaf INSCAF ............ input file of scafs, numbered sequentially\n"
      +"    -out_scaf OUTSCAF .......... same scafs; stats appended to titles\n"
      +"  options:\n"
      +"    -scaflist_title  ........... scaf/link/chain list is title\n"
      +"    -scaflist_append2title ..... scaf/link/chain is 2nd title field\n"
      +"    -scaflist_sdtag SDTAG ...... scaf list input SD dataitem\n"
      +"    -report_mol_lists_ids ...... after stats, append list of mol ids (maybe fat)\n"
      +"    -sort_by_frequency ......... output scaf/link/chain lists sorted\n"
      +"    -v ......................... verbose\n"
      +"    -vv ........................ very verbose\n"
      +"    -h ......................... this help\n");
    System.exit(1);
  }
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static int verbose=0;
  private static String ifile_mol=null;
  private static String ifile_scaf=null;
  private static String ofile_scaf=null;
  //private static String ifile_link=null;
  //private static String ofile_link=null;
  //private static String ifile_chain=null;
  //private static String ofile_chain=null;
  private static Boolean scaflist_title=false;
  private static Boolean scaflist_append2title=false;
  private static String scaflist_sdtag="HSCAF_SCAFLIST";
  //private static String linklist_sdtag="HSCAF_LINKLIST";
  //private static String chainlist_sdtag="HSCAF_CHAINLIST";
  private static Boolean report_mol_lists_ids=false;
  private static Boolean sort_by_frequency=false;

  /////////////////////////////////////////////////////////////////////////////
  private static void parseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) ifile_mol=args[++i];
      else if (args[i].equals("-in_scaf")) ifile_scaf=args[++i];
      //else if (args[i].equals("-in_link")) ifile_link=args[++i];
      //else if (args[i].equals("-in_chain")) ifile_chain=args[++i];
      else if (args[i].equals("-out_scaf")) ofile_scaf=args[++i];
      //else if (args[i].equals("-out_link")) ofile_link=args[++i];
      //else if (args[i].equals("-out_chain")) ofile_chain=args[++i];
      else if (args[i].equals("-scaflist_title")) scaflist_title=true;
      else if (args[i].equals("-scaflist_append2title")) scaflist_append2title=true;
      else if (args[i].equals("-scaflist_sdtag")) scaflist_sdtag=args[++i];
      //else if (args[i].equals("-linklist_sdtag")) linklist_sdtag=args[++i];
      //else if (args[i].equals("-chainlist_sdtag")) chainlist_sdtag=args[++i];
      else if (args[i].equals("-report_mol_lists_ids")) report_mol_lists_ids=true;
      else if (args[i].equals("-sort_by_frequency")) sort_by_frequency=true;
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-vvv") || args[i].equals("-debug")) verbose=3;
      else if (args[i].equals("-h")) help("");
      else help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException
  {
    parseCommand(args);
    if (ifile_mol==null || ifile_scaf==null) help("-i and -in_scaf required.");
    if (!(new File(ifile_mol).exists())) help("Non-existent input file: "+ifile_mol);
    MolImporter molReader_mol = new MolImporter(ifile_mol);
    if (!(new File(ifile_scaf).exists())) help("Non-existent input file: "+ifile_scaf);
    MolImporter molReader_scaf = new MolImporter(ifile_scaf);

    if (ofile_scaf==null) help("-out_scaf required.");
    String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
    if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
    MolExporter molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf),ofmt);

    //MolImporter molReader_link=null;
    //if (ifile_link!=null)
    //{
    //  if (!(new File(ifile_link).exists())) help("Non-existent input file: "+ifile_link);
    //  molReader_link = new MolImporter(ifile_link);
    //}
    //MolExporter molWriter_link=null;
    //if (ofile_link!=null)
    //{
    //  ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_link);
    //  if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
    //  molWriter_link=new MolExporter(new FileOutputStream(ofile_link),ofmt);
    //}
    //
    //MolImporter molReader_chain=null;
    //if (ifile_chain!=null)
    //{
    //  if (!(new File(ifile_chain).exists())) help("Non-existent input file: "+ifile_chain);
    //  molReader_chain = new MolImporter(ifile_chain);
    //}
    //MolExporter molWriter_chain=null;
    //if (ofile_chain!=null)
    //{
    //  ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_chain);
    //  if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
    //  molWriter_chain=new MolExporter(new FileOutputStream(ofile_chain),ofmt);
    //}

    if (verbose>1)
      System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);

    // Read scaf file (required).
    LinkedHashMap<Long,ScaffoldInfo> scafs = new LinkedHashMap<Long,ScaffoldInfo>();
    ArrayList<ScaffoldInfo> scafs_list = new ArrayList<ScaffoldInfo>(); //for sorting only
    Molecule scafmol;
    long n_scaf=0;
    int n_err=0;
    java.util.Date t_0 = new java.util.Date();
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(t_0));
    java.util.Date t_i = t_0;
    for (n_scaf=0;true;)
    {
      try { scafmol=molReader_scaf.read(); }
      catch (MolFormatException e)
      {
        System.err.println(e.getMessage());
        ++n_err;
        continue;
      }
      if (scafmol==null) break; //EOF
      ++n_scaf;

      String molname=scafmol.getName();
      String scafsmi=MolExporter.exportToFormat(scafmol,smifmt);
      if (verbose>1)
        System.err.println(""+n_scaf+". "+molname+"\n\t"+scafsmi);

      String[] fields=Pattern.compile("\\s+").split(molname);
      if (fields.length!=2)
      {
        System.err.println("ERROR: Bad line ["+n_scaf+"]; nfields!=2 ; "+molname);
        ++n_err;
        continue;
      }

      long scafid;
      try {
        scafid=Long.parseLong(fields[0]);
      }
      catch (Exception e) {
        ++n_err;
        System.err.println("ERROR: Bad line ["+n_scaf+"]; cannot parse ScafID; "+molname);
        continue;
      }

      ScaffoldInfo scafinfo = new ScaffoldInfo(n_scaf);
      scafinfo.setSmiles(scafsmi);
      scafinfo.setTreeString(fields[1]);
      scafs.put(n_scaf,scafinfo);
      scafs_list.add(scafinfo);
    }
    System.err.println("Total scafs: "+n_scaf);

    //// Read link file (optional).
    //LinkedHashMap<Long,LinkerInfo> links = null;
    //ArrayList<LinkerInfo> links_list = null;
    //long n_link=0;
    //Molecule linkmol=null;
    //if (ifile_link!=null)
    //{
    //  links = new LinkedHashMap<Long,LinkerInfo>();
    //  links_list = new ArrayList<LinkerInfo>(); //for sorting only
    //  if (verbose>0)
    //    System.err.println(DateFormat.getDateTimeInstance().format(t_0));
    //  for (n_link=0;true;)
    //  {
    //    try { linkmol=molReader_link.read(); }
    //    catch (MolFormatException e)
    //    {
    //      System.err.println(e.getMessage());
    //      ++n_err;
    //      continue;
    //    }
    //    if (linkmol==null) break; //EOF
    //    ++n_link;
    //
    //    String molname=linkmol.getName();
    //    String linksmi=MolExporter.exportToFormat(linkmol,smifmt);
    //    if (verbose>1)
    //      System.err.println(""+n_link+". "+molname+"\n\t"+linksmi);
    //
    //    long linkid;
    //    try {
    //      linkid=Long.parseLong(molname);
    //    }
    //    catch (Exception e) {
    //      ++n_err;
    //      System.err.println("ERROR: Bad line ["+n_link+"]; cannot parse LinkID; "+molname);
    //      continue;
    //    }
    //    LinkerInfo linkinfo = new LinkerInfo(n_link);
    //    linkinfo.setSmiles(linksmi);
    //    links.put(n_link,linkinfo);
    //    links_list.add(linkinfo);
    //  }
    //  System.err.println("Total linkers: "+n_link);
    //}
    //
    //// Read chain file (optional).
    //LinkedHashMap<Long,ChainInfo> chains = null;
    //ArrayList<ChainInfo> chains_list = null;
    //long n_chain=0;
    //Molecule chainmol=null;
    //if (ifile_chain!=null)
    //{
    //  chains = new LinkedHashMap<Long,ChainInfo>();
    //  chains_list = new ArrayList<ChainInfo>(); //for sorting only
    //  if (verbose>0)
    //    System.err.println(DateFormat.getDateTimeInstance().format(t_0));
    //  for (n_chain=0;true;)
    //  {
    //    try { chainmol=molReader_chain.read(); }
    //    catch (MolFormatException e)
    //    {
    //      System.err.println(e.getMessage());
    //      ++n_err;
    //      continue;
    //    }
    //    if (chainmol==null) break; //EOF
    //    ++n_chain;
    //
    //    String molname=chainmol.getName();
    //    String chainsmi=MolExporter.exportToFormat(chainmol,smifmt);
    //    if (verbose>1)
    //      System.err.println(""+n_chain+". "+molname+"\n\t"+chainsmi);
    //
    //    long chainid;
    //    try {
    //      chainid=Long.parseLong(molname);
    //    }
    //    catch (Exception e) {
    //      ++n_err;
    //      System.err.println("ERROR: Bad line ["+n_chain+"]; cannot parse ChainID; "+molname);
    //      continue;
    //    }
    //    ChainInfo chaininfo = new ChainInfo(n_chain);
    //    chaininfo.setSmiles(chainsmi);
    //    chains.put(n_chain,chaininfo);
    //    chains_list.add(chaininfo);
    //  }
    //  System.err.println("Total chains: "+n_chain);
    //}

    // Read mol file.  Parse scaflist; annotate scaffoldInfo list.

    Molecule mol;
    int n_mol=0;
    int n_noscaf=0;
    //int n_nolink=0;
    //int n_nochain=0;
    for (n_mol=0;true;)
    {
      try { mol=molReader_mol.read(); }
      catch (MolFormatException e)
      {
        System.err.println(e.getMessage());
        ++n_err;
        continue;
      }
      if (mol==null) break; //EOF
      ++n_mol;

      String molname=mol.getName();
      String smi=MolExporter.exportToFormat(mol,smifmt);
      if (verbose>1)
        System.err.println(""+n_mol+". "+molname+"\n\t"+smi);

      // molname format:
      // <CID><space>S:[ID1,ID2...]<space>L:[ID1,ID2...]<space>C:[ID1,ID2...]

      if (!molname.matches("[0-9]+\\s+S:.*$"))
      {
        if (verbose>1)
          System.err.println("ERROR: Bad line (no hier_scaffolds output) ["+n_mol+"] "+molname);
        ++n_noscaf;
        continue;
      }
      String[] fields=Pattern.compile("\\s+").split(molname);
      if (fields.length<2)
      {
        System.err.println("ERROR: Bad line (no ScafIDs field) ["+n_mol+"]; nfields<2; "+molname);
        ++n_noscaf;
        continue;
      }

      Long cid=null;
      try {
        cid=Long.parseLong(fields[0]);
      }
      catch (Exception e) {
        ++n_err;
        System.err.println("ERROR: Bad line ["+n_mol+"]; cannot parse CID; "+molname);
        continue;
      }

      //scaffolds:
      String[] scafids_str=Pattern.compile(",").split(fields[1].replaceFirst("^S:",""));
      for (String scafid_str: scafids_str)
      {
        if (scafid_str.isEmpty()) continue;
        Long scafid=null;
        try {
          scafid=Long.parseLong(scafid_str);
          if (!scafs.containsKey(scafid))
          {
            System.err.println("ERROR: ScafID "+scafid+" not known.");
            ++n_err;
            continue;
          }
          scafs.get(scafid).addCID(cid);
        }
        catch (Exception e)
        { System.err.println("ERROR: ["+n_mol+"]; cannot parse ScafID (\""+scafid_str+"\"); "+molname); }
      }

    //  //linkers:
    //  if (ifile_link!=null)
    //  {
    //    if (fields.length<3)
    //    {
    //      System.err.println("ERROR: Bad line (no LinkIDs field) ["+n_mol+"]; nfields<3; "+molname);
    //      ++n_nolink;
    //      continue;
    //    }
    //    String[] linkids_str=Pattern.compile(",").split(fields[2].replaceFirst("^L:",""));
    //    for (String linkid_str: linkids_str)
    //    {
    //      if (linkid_str.isEmpty()) continue;
    //      Long linkid=null;
    //      try {
    //        linkid=Long.parseLong(linkid_str);
    //        if (!links.containsKey(linkid))
    //        {
    //          System.err.println("ERROR: LinkID "+linkid+" not known.");
    //          ++n_err;
    //          continue;
    //        }
    //        links.get(linkid).addCID(cid);
    //      }
    //      catch (Exception e)
    //      { System.err.println("ERROR: ["+n_mol+"]; cannot parse LinkID (\""+linkid_str+"\"); "+molname); }
    //    }
    //  }

    //  //chains:
    //  if (ifile_chain!=null)
    //  {
    //    if (fields.length<4)
    //    {
    //      System.err.println("ERROR: Bad line (no ChainIDs field) ["+n_mol+"]; nfields<4; "+molname);
    //      ++n_nochain;
    //      continue;
    //    }
    //    String[] chainids_str=Pattern.compile(",").split(fields[3].replaceFirst("^C:",""));
    //    for (String chainid_str: chainids_str)
    //    {
    //      if (chainid_str.isEmpty()) continue;
    //      Long chainid=null;
    //      try {
    //        chainid=Long.parseLong(chainid_str);
    //        if (!chains.containsKey(chainid))
    //        {
    //          System.err.println("ERROR: ChainID "+chainid+" not known.");
    //          ++n_err;
    //          continue;
    //        }
    //        chains.get(chainid).addCID(cid);
    //      }
    //      catch (Exception e)
    //      { System.err.println("ERROR: ["+n_mol+"]; cannot parse ChainID (\""+chainid_str+"\"); "+molname); }
    //    }
    //  }

    }
    System.err.println("Total mols: "+n_mol);

    if (sort_by_frequency)
    {
      Collections.sort(scafs_list);
    //  if (ifile_link!=null) Collections.sort(links_list);
    //  if (ifile_chain!=null) Collections.sort(chains_list);
    }

    // scaffolds:
    long scaf_freq_total=0;
    for (ScaffoldInfo scaf: scafs_list)
    {
      long scafid=scaf.getID();
      ArrayList<Long> cids = scaf.getCIDs();
      scaf_freq_total+=cids.size();

      try { scafmol=MolImporter.importMol(scaf.getSmiles(),"smiles:"); }
      catch (MolFormatException e) { }

      // format: <SCAFID> <TREESTR> <NCPD> [<CIDLIST>]
      String newname=(""+scafid+" "+scaf.getTreeString()+" "+cids.size());
      if (report_mol_lists_ids)
      {
        String cidlist_str="";
        if (cids.size()>0)
        {
          Collections.sort(cids);
          cidlist_str=(""+cids.get(0));
          for (int j=1;j<cids.size();++j) cidlist_str+=(","+cids.get(j));
        }
        newname+=(" "+cidlist_str);
      }
      scafmol.setName(newname);
      writeMol(molWriter_scaf,scafmol);
    }
    molWriter_scaf.close();

    //// linkers:
    //long link_freq_total=0;
    //if (ifile_link!=null) 
    //{
    //  for (LinkerInfo link: links_list)
    //  {
    //    long linkid=link.getID();
    //    ArrayList<Long> cids = link.getCIDs();
    //    link_freq_total+=cids.size();
    //
    //    try { linkmol=MolImporter.importMol(link.getSmiles(),"smiles:"); }
    //    catch (MolFormatException e) { }
    //
    //    // format: <LINKID> <NCPD> [<CIDLIST>]
    //    String newname=(""+linkid+" "+cids.size());
    //    if (report_mol_lists_ids)
    //    {
    //      String cidlist_str="";
    //      if (cids.size()>0)
    //      {
    //        Collections.sort(cids);
    //        cidlist_str=(""+cids.get(0));
    //        for (int j=1;j<cids.size();++j) cidlist_str+=(","+cids.get(j));
    //      }
    //      newname+=(" "+cidlist_str);
    //    }
    //    linkmol.setName(newname);
    //    if (ofile_link!=null) writeMol(molWriter_link,linkmol);
    //  }
    //  if (ofile_link!=null) molWriter_link.close();
    //}
    //
    //// chains:
    //long chain_freq_total=0;
    //if (ifile_chain!=null) 
    //{
    //  for (ChainInfo chain: chains_list)
    //  {
    //    long chainid=chain.getID();
    //    ArrayList<Long> cids = chain.getCIDs();
    //    chain_freq_total+=cids.size();
    //
    //    try { chainmol=MolImporter.importMol(chain.getSmiles(),"smiles:"); }
    //    catch (MolFormatException e) { }
    //
    //    // format: <LINKID> <NCPD> [<CIDLIST>]
    //    String newname=(""+chainid+" "+cids.size());
    //    if (report_mol_lists_ids)
    //    {
    //      String cidlist_str="";
    //      if (cids.size()>0)
    //      {
    //        Collections.sort(cids);
    //        cidlist_str=(""+cids.get(0));
    //        for (int j=1;j<cids.size();++j) cidlist_str+=(","+cids.get(j));
    //      }
    //      newname+=(" "+cidlist_str);
    //    }
    //    chainmol.setName(newname);
    //    if (ofile_chain!=null) writeMol(molWriter_chain,chainmol);
    //  }
    //  if (ofile_chain!=null) molWriter_chain.close();
    //}

    System.err.println("Total elapsed time: "+time_utils.timeDeltaStr(t_0,new java.util.Date()));
    if (verbose>0)
      System.err.println(DateFormat.getDateTimeInstance().format(new java.util.Date()));
    System.err.println("Total scaffold occurances: "+scaf_freq_total);
    //if (ifile_link!=null)
    //  System.err.println("Total linker occurances: "+link_freq_total);
    //if (ifile_chain!=null)
    //  System.err.println("Total chain occurances: "+chain_freq_total);
    System.err.println("Errors: "+n_err);
    System.err.println("Output scaffolds, annotated: "+ofile_scaf);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static boolean writeMol(MolExporter molWriter,Molecule mol)
  {
    try {
      molWriter.write(mol);
    }
    catch (Exception e) {
      System.err.println(e.getMessage());
      return false;
    }
    return true;
  }
}
/////////////////////////////////////////////////////////////////////////////
class ScaffoldInfo implements Comparable<Object>
{
  Long id;
  String smiles;
  String treestr;
  ArrayList<Long> cids;
  ScaffoldInfo(Long _id)
  {
    this.id=_id;
  }
  void setSmiles(String _smiles) { this.smiles=_smiles; }
  String getSmiles() { return this.smiles; }
  void setTreeString(String _treestr) { this.treestr=_treestr; }
  String getTreeString() { return this.treestr; }
  Long getID() { return this.id; }
  Integer addCID(Long _cid)
  {
    if (this.cids==null) this.cids = new ArrayList<Long>();
    this.cids.add(_cid);
    return this.cids.size();
  }
  ArrayList<Long> getCIDs()
  {
    return this.cids;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   For sorting by frequency&#46;
  */
  public int compareTo(Object o) throws ClassCastException
  {
    return ((ScaffoldInfo)o).cids.size()-this.cids.size();
  }
}
///////////////////////////////////////////////////////////////////////////////
//class LinkerInfo implements Comparable<Object>
//{
//  Long id;
//  String smiles;
//  ArrayList<Long> cids;
//  LinkerInfo(Long _id)
//  {
//    this.id=_id;
//  }
//  void setSmiles(String _smiles) { this.smiles=_smiles; }
//  String getSmiles() { return this.smiles; }
//  Long getID() { return this.id; }
//  Integer addCID(Long _cid)
//  {
//    if (this.cids==null) this.cids = new ArrayList<Long>();
//    this.cids.add(_cid);
//    return this.cids.size();
//  }
//  ArrayList<Long> getCIDs()
//  {
//    return this.cids;
//  }
//  /////////////////////////////////////////////////////////////////////////////
//  /**   For sorting by frequency&#46;
//  */
//  public int compareTo(Object o) throws ClassCastException
//  {
//    return ((LinkerInfo)o).cids.size()-this.cids.size();
//  }
//}
///////////////////////////////////////////////////////////////////////////////
//class ChainInfo implements Comparable<Object>
//{
//  Long id;
//  String smiles;
//  ArrayList<Long> cids;
//  ChainInfo(Long _id)
//  {
//    this.id=_id;
//  }
//  void setSmiles(String _smiles) { this.smiles=_smiles; }
//  String getSmiles() { return this.smiles; }
//  Long getID() { return this.id; }
//  Integer addCID(Long _cid)
//  {
//    if (this.cids==null) this.cids = new ArrayList<Long>();
//    this.cids.add(_cid);
//    return this.cids.size();
//  }
//  ArrayList<Long> getCIDs()
//  {
//    return this.cids;
//  }
//  /////////////////////////////////////////////////////////////////////////////
//  /**   For sorting by frequency&#46;
//  */
//  public int compareTo(Object o) throws ClassCastException
//  {
//    return ((ChainInfo)o).cids.size()-this.cids.size();
//  }
//}
