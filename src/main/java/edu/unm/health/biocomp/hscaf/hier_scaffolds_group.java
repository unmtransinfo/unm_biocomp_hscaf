package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.text.DateFormat;

import chemaxon.formats.*;
import chemaxon.sss.search.*;
import chemaxon.calculations.*; //Ring
import chemaxon.struc.*;
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.*;

/**	HierS hierarchical scaffolds postprocess application.  Run program with
	no args for command-line help.  
	Group molecules by scaffold membership (normally using largest scaffold only).
	<br>
	Molecule files normally should be SMILES format.  The SMILES "name" is a space delimited
	set of fields: (1) compound ID (integer), (2) set of scaffold IDs.  For example:
	FC1=CC=C(NC2=NC(=NC3=CC=CC=C23)C2=CN=CC=C2)C=C1 1877191 S:4,5,6,17,18
	<br>
	Input scaffold file also normally SMILES, with "name" delimited set of fields:
	(1) scaffold ID, and (2) scaffold tree string.  For example:
	N(C1=CC=NC=C1)C1=C2C=CC=CC2=NC(=N1)C1=COC=C1    35 35:(38:(5,6),6,36:(5,37),37)
	<br>
	Note: for NCI/CBC used hier_scaffolds_groupmolsbyscaf.py.
	<br>
	@author Jeremy J Yang
*/
public class hier_scaffolds_group
{
  private static int MIN_RINGS=1;

  private static void Help(String msg)
  {
    System.err.println(msg+"\n"
      +"hier_scaffolds_group - group mols by scaffold (for each mol largest scaf only)\n"
      +"\n"
      +"From the University of New Mexico Division Translational Informatics Division\n"
      +"(http://medicine.unm.edu/informatics/).\n"
      +"\n"
      +"usage: hier_scaffolds_group [options]\n"
      +"  required:\n"
      +"    -i IFILE ................... input file of molecules annotated with scaf info\n"
      +"    -in_scaf INSCAF ............ input file of scafs, numbered sequentially\n"
      +"    -out_scaf OUTSCAF .......... same scafs; mol name lists appended\n"
      +"    -out_scaf_reject OUTSCAF ... \n"
      +"  options:\n"
      +"    -scaflist_sdtag SDTAG ...... scaf list input SD dataitem\n"
      +"    -scaflist_title  ........... scaf/link/chain list is title\n"
      +"    -scaflist_append2title ..... scaf/link/chain is 2nd title field\n"
      +"    -min_rings N ............... min rings fused into scaffold ["+MIN_RINGS+"]\n"
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
  private static String ofile_scaf_reject=null;
  private static Boolean scaflist_title=false;
  private static Boolean scaflist_append2title=false;
  private static String scaflist_sdtag="HSCAF_SCAFLIST";

  /////////////////////////////////////////////////////////////////////////////
  private static void ParseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) ifile_mol=args[++i];
      else if (args[i].equals("-in_scaf")) ifile_scaf=args[++i];
      else if (args[i].equals("-out_scaf")) ofile_scaf=args[++i];
      else if (args[i].equals("-out_scaf_reject")) ofile_scaf_reject=args[++i];
      else if (args[i].equals("-scaflist_title")) scaflist_title=true;
      else if (args[i].equals("-scaflist_append2title")) scaflist_append2title=true;
      else if (args[i].equals("-scaflist_sdtag")) scaflist_sdtag=args[++i];
      else if (args[i].equals("-min_rings")) MIN_RINGS=Integer.parseInt(args[++i]);
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-vvv") || args[i].equals("-debug")) verbose=3;
      else if (args[i].equals("-h")) Help("");
      else Help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  private static Long LargestScaffold(ArrayList<Long> scafids,LinkedHashMap<Long,ScaffoldInfo2> scafs)
  {
    Long scafid_best=null;
    int smi_len_max=0;
    for (int i=0;i<scafids.size();++i)
    {
      ScaffoldInfo2 scaf=scafs.get(scafids.get(i));
      if (scaf!=null)
      {
        String smi=scaf.getSmiles();
        if (smi.length()>smi_len_max)
        {
          smi_len_max=smi.length();
          scafid_best=scafids.get(i);
        }
      }
    }
    return scafid_best;
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException
  {
    ParseCommand(args);
    if (ifile_mol==null || ifile_scaf==null) Help("-i and -in_scaf required.");
    if (!(new File(ifile_mol).exists())) Help("Non-existent input file: "+ifile_mol);
    MolImporter molReader_mol = new MolImporter(ifile_mol);
    if (!(new File(ifile_scaf).exists())) Help("Non-existent input file: "+ifile_scaf);
    MolImporter molReader_scaf = new MolImporter(ifile_scaf);

    if (ofile_scaf==null) Help("-out_scaf required.");
    String ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile_scaf);
    if (ofmt.equals("smiles")) ofmt="smiles:+n-a"; //Kekule for compatibility
    MolExporter molWriter_scaf=new MolExporter(new FileOutputStream(ofile_scaf),ofmt);

    //if (verbose>1)
      //System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.JCHEM_VERSION);
      //System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion()); //v6.3+

    // Read scaf file (required).
    // Format SCAFSMI SCAFID [SCAFTREE]
    // Maybe: Filter unqualified scaffolds (too small, too few rings?)
    LinkedHashMap<Long,ScaffoldInfo2> scafs = new LinkedHashMap<Long,ScaffoldInfo2>();
    Molecule scafmol;
    Ring ring = new Ring();	//ring analyzer
    long n_scaf=0;
    long n_scaf_ok=0;
    long n_scaf_minrings_notok=0;
    int n_err=0;
    long scafid_max=0;
    java.util.Date t_0 = new java.util.Date();
    if (verbose>1)
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

      ring.setMolecule(scafmol);
      if (ring.ringCount()<MIN_RINGS)
      {
        //System.err.println("DEBUG: min_rings violation; ["+ring.ringCount()+"<"+MIN_RINGS+"]");
        ++n_scaf_minrings_notok;
        continue;
      }
      ++n_scaf_ok;

      String scafsmi=MolExporter.exportToFormat(scafmol,smifmt);
      if (verbose>1)
        System.err.println(""+n_scaf+". "+molname+"\n\t"+scafsmi);

      String[] fields=Pattern.compile("\\s+").split(molname);
      if (fields.length<1)
      {
        System.err.println("ERROR: Bad line ["+n_scaf+"]; nfields<1 ; "+molname);
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

      if (scafid>scafid_max) scafid_max=scafid;

      ScaffoldInfo2 scafinfo = new ScaffoldInfo2(n_scaf);
      scafinfo.setSmiles(scafsmi);
      if (fields.length>1)
        scafinfo.setTreeString(fields[1]);
      scafs.put(scafid,scafinfo);
    }

    // Read mol file.  Parse scaflist; annotate scaffoldInfo list.
    Molecule mol;
    int n_mol=0;
    int n_noscaf=0;
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
      // 1st field must be a unique ID integer (CID).
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
        System.err.println("ERROR: Bad line ["+n_mol+"]; nfields<2; "+molname);
        ++n_noscaf;
        continue;
      }

      Long cid=null;
      try { cid=Long.parseLong(fields[0]); }
      catch (Exception e) {
        ++n_err;
        System.err.println("ERROR: Bad line ["+n_mol+"]; cannot parse CID; "+molname);
        continue;
      }

      // scaffolds: For each mol, we ignore all but the largest scaffold.
      // Note list in descending size (smiles length) order, so take 1st.
      String[] scafids_str=Pattern.compile(",").split(fields[1].replaceFirst("^S:",""));
      ArrayList<Long> scafids = new ArrayList<Long>();
      for (String scafid_str: scafids_str)
      {
        if (scafid_str.isEmpty()) continue;
        Long scafid=null;
        try { scafid=Long.parseLong(scafid_str); }
        catch (Exception e)
        { System.err.println("ERROR: ["+n_mol+"]; cannot parse ScafID (\""+scafid_str+"\"); "+molname); }
        scafids.add(scafid);
      }
      if (scafids.size()==0)
      {
        ++n_noscaf;
        continue;
      }
      Long scafid=LargestScaffold(scafids,scafs);
      if (scafid==null)
      {
        continue; //rejected scaf (min_rings violation)
      }
      else if (!scafs.containsKey(scafid))
      {
        System.err.println("ERROR: ScafID "+scafid+" not known.");
        ++n_err;
        continue;
      }
      scafs.get(scafid).addCID(cid);
    }

    //JJ: Write annotated scaf_groups file...
    // scaffolds:
    for (long scafid=1;scafid<=scafid_max;++scafid)
    {
      if (!scafs.containsKey(scafid)) continue;
      ScaffoldInfo2 scaf=scafs.get(scafid);
      ArrayList<Long> cids = scaf.getCIDs();

      try { scafmol=MolImporter.importMol(scaf.getSmiles(),"smiles:"); }
      catch (MolFormatException e) { }

      // format: <SCAFID> <TREESTR> <NCPD> [<CIDLIST>]
      String newname=(""+scafid+" "+scaf.getTreeString()+" "+cids.size());
      
      String cidlist_str="";
      if (cids.size()>0)
      {
        Collections.sort(cids);
        cidlist_str=(""+cids.get(0));
        for (int j=1;j<cids.size();++j) cidlist_str+=(","+cids.get(j));
      }
      newname+=(" "+cidlist_str);

      scafmol.setName(newname);
      WriteMol(molWriter_scaf,scafmol);
    }
    molWriter_scaf.close();

    System.err.println("total elapsed time: "+time_utils.TimeDeltaStr(t_0,new java.util.Date()));
    if (verbose>1)
      System.err.println(DateFormat.getDateTimeInstance().format(new java.util.Date()));
    System.err.println("total scafs: "+n_scaf);
    System.err.println("total ok scafs: "+n_scaf_ok);
    if (scafid_max!=n_scaf)
      System.err.println("Warning: n_scaf!=scafid_max ["+n_scaf+"!="+scafid_max+"]; possible numbering error.");
    System.err.println("total mols: "+n_mol);
    System.err.println("mols with no scaffolds: "+n_noscaf);
    System.err.println("output scaffolds, annotated: "+ofile_scaf);
    System.err.println("scaffolds rejected by min_rings [<"+MIN_RINGS+"] violation: "+n_scaf_minrings_notok);
    System.err.println("Errors: "+n_err);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static boolean WriteMol(MolExporter molWriter,Molecule mol)
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
class ScaffoldInfo2
{
  Long id;
  String smiles;
  String treestr;
  ArrayList<Long> cids;
  ScaffoldInfo2(Long _id)
  {
    this.cids = new ArrayList<Long>();
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
}
