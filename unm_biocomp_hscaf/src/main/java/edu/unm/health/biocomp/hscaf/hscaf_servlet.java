package edu.unm.health.biocomp.hscaf;

import java.io.*;
import java.lang.Math;
import java.net.*; //URLEncoder,InetAddress
import java.text.*;
import java.util.*;
import java.util.regex.*;
import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;

import com.oreilly.servlet.*; //MultipartRequest,Base64Encoder,Base64Decoder
import com.oreilly.servlet.multipart.DefaultFileRenamePolicy;

import chemaxon.formats.*;
import chemaxon.util.MolHandler;
import chemaxon.struc.*;
import chemaxon.sss.search.SearchException;
import chemaxon.license.*; //LicenseException,LicenseManager
import chemaxon.marvin.io.MolExportException;

import edu.unm.health.biocomp.util.*;
import edu.unm.health.biocomp.util.http.*;

/**	HScaf web app.  Hierarchical scaffold analysis (Wilkens, 2005).
	<br>
	@author Jeremy J Yang
*/
@WebServlet(urlPatterns = "/hscaf", name = "hscaf",
  initParams = {
    @WebInitParam(name = "PROXY_PREFIX", value = "/tomcat"), //Could be in web.xml?
    @WebInitParam(name = "APPNAME", value = "HScaf"),
    @WebInitParam(name = "UPLOADDIR", value = "/tmp"),
    @WebInitParam(name = "SCRATCHDIR", value = "/scratch"),
    @WebInitParam(name = "N_MAX", value = "200"),
    @WebInitParam(name = "MAX_POST_SIZE", value = "10485760"),
    @WebInitParam(name = "DEMOSMIFILE", value = "/data/hscaf_testset.smi")
	}
)
public class hscaf_servlet extends HttpServlet
{
  private static String SERVLETNAME=null;
  private static String CONTEXTPATH=null;
  private static String APPNAME=null;	// configured in web.xml
  private static String UPLOADDIR=null;	// configured in web.xml
  private static String SCRATCHDIR=null;      // configured in web.xml
  private static String SCRATCHDIR_PATH=null;
  private static String DEMOSMIFILE=null;      // configured in web.xml
  private static String PREFIX=null;
  private static int scratch_retire_sec=3600;
  private static int N_MAX=100; // configured in web.xml
  private static ServletContext CONTEXT=null;
  //private static ServletConfig CONFIG=null;
  private static ResourceBundle rb=null;
  private static PrintWriter out=null;
  private static ArrayList<String> outputs=null;
  private static ArrayList<String> errors=null;
  private static HttpParams params=null;
  private static ArrayList<Molecule> mols=null;
  private static LinkedHashMap<String,Integer> sizes_h=null;
  private static LinkedHashMap<String,Integer> sizes_w=null;
  private static int SERVERPORT=0;
  private static String SERVERNAME=null;
  private static String REMOTEHOST=null;
  private static String DATESTR=null;
  private static String color1="#EEEEEE";
  private static String smifmt="cxsmiles:u-L-l-e-d-D-p-R-f-w";
  private static String smifmt_dep="cxsmiles:u-L-e-d-D-p-R-f-w"; // "l" for labels/aliases ("J")
  private static File fout=null;
  private static File fout_scaf=null;
  private static String PROXY_PREFIX=null; // configured in web.xml
  private static String MOL2IMG_SERVLETURL=null;

  /////////////////////////////////////////////////////////////////////////////
  public void doPost(HttpServletRequest request,HttpServletResponse response)
      throws IOException,ServletException
  {
    SERVERPORT=request.getServerPort();
    SERVERNAME=request.getServerName();
    if (SERVERNAME.equals("localhost")) SERVERNAME=InetAddress.getLocalHost().getHostAddress();
    REMOTEHOST=request.getHeader("X-Forwarded-For"); // client (original)
    if (REMOTEHOST!=null)
    {
      String[] addrs=Pattern.compile(",").split(REMOTEHOST);
      if (addrs.length>0) REMOTEHOST=addrs[addrs.length-1];
    }
    else
    {
      REMOTEHOST=request.getRemoteAddr(); // client (may be proxy)
    }
    rb=ResourceBundle.getBundle("LocalStrings",request.getLocale());

    MultipartRequest mrequest=null;
    if (request.getMethod().equalsIgnoreCase("POST"))
    {
      try { mrequest=new MultipartRequest(request,UPLOADDIR,10*1024*1024,"ISO-8859-1",
                                    new DefaultFileRenamePolicy()); }
      catch (IOException lEx) {
        this.getServletContext().log("not a valid MultipartRequest",lEx); }
    }

    // main logic:
    ArrayList<String> cssincludes = new
ArrayList<String>(Arrays.asList(PROXY_PREFIX+CONTEXTPATH+"/css/biocomp.css"));
    ArrayList<String> jsincludes = new ArrayList<String>(Arrays.asList(PROXY_PREFIX+CONTEXTPATH+"/js/biocomp.js",PROXY_PREFIX+CONTEXTPATH+"/js/ddtip.js"));
    boolean ok=initialize(request,mrequest);
    if (!ok)
    {
      response.setContentType("text/html");
      out=response.getWriter();
      out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
      out.println(HtmUtils.FooterHtm(errors,true));
      return;
    }
    else if (mrequest!=null)	//method=POST, normal operation
    {
      if (mrequest.getParameter("hscaf").equals("TRUE"))
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(formHtm(mrequest,response));
        ArrayList<ScaffoldTree> scaftrees = HScafsGenerate(mols);

        if (params.getVal("runmode").equals("1xN"))
        {
          ArrayList<Scaffold> mcscafs = MCScafs_1xN(scaftrees);
          HScafsResults_1xN(response,scaftrees,mcscafs);
        }
        else if (params.getVal("runmode").equals("NxN"))
        {
          ArrayList<ArrayList<Scaffold> > mcscaf_matrix = MCScafs_NxN(scaftrees);
          HScafsResults_NxN(response,scaftrees,mcscaf_matrix);
        }
        else	//runmode "analyze"
        {
          HScafsResults(response,scaftrees);
        }
        out.println(HtmUtils.OutputHtm(outputs));
        out.println(HtmUtils.FooterHtm(errors,true));
        HtmUtils.PurgeScratchDirs(Arrays.asList(SCRATCHDIR_PATH),scratch_retire_sec,false,".",(HttpServlet) this);
      }
    }
    else
    {
      String help=request.getParameter("help");	// GET param
      String downloadtxt=request.getParameter("downloadtxt"); // POST param
      String downloadfile=request.getParameter("downloadfile"); // POST param
      if (help!=null)	// GET method, help=TRUE
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(HelpHtm());
        out.println(HtmUtils.FooterHtm(errors,true));
      }
      else if (downloadtxt!=null && downloadtxt.length()>0) // POST param
      {
        ServletOutputStream ostream=response.getOutputStream();
        HtmUtils.DownloadString(response,ostream,downloadtxt,
          request.getParameter("fname"));
      }
      else if (downloadfile!=null && downloadfile.length()>0) // POST param
      {
        ServletOutputStream ostream=response.getOutputStream();
        HtmUtils.DownloadFile(response,ostream,downloadfile,
          request.getParameter("fname"));
      }
      else	// GET method, initial invocation of servlet w/ no params
      {
        response.setContentType("text/html");
        out=response.getWriter();
        out.println(HtmUtils.HeaderHtm(APPNAME, jsincludes, cssincludes, JavaScript(), "", color1, request));
        out.println(formHtm(mrequest,response));
        out.println("<SCRIPT>go_init(window.document.mainform)</SCRIPT>");
        out.println(HtmUtils.FooterHtm(errors,true));
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  private boolean initialize(HttpServletRequest request,MultipartRequest mrequest)
      throws IOException,ServletException
  {
    SERVLETNAME=this.getServletName();
    outputs=new ArrayList<String>();
    errors=new ArrayList<String>();
    mols=new ArrayList<Molecule>();
    params=new HttpParams();
    sizes_h=new LinkedHashMap<String,Integer>();
    sizes_w=new LinkedHashMap<String,Integer>();
    sizes_h.put("s",160); sizes_w.put("s",160);
    sizes_h.put("m",200); sizes_w.put("m",220);
    sizes_h.put("l",240); sizes_w.put("l",280);
    MOL2IMG_SERVLETURL=(PROXY_PREFIX+CONTEXTPATH+"/mol2img");

    String logo_htm="<TABLE CELLSPACING=5 CELLPADDING=5><TR><TD>";
    String imghtm=("<IMG BORDER=0 SRC=\""+PROXY_PREFIX+CONTEXTPATH+"/images/biocomp_logo_only.gif\">");
    String tiphtm=(APPNAME+" web app from UNM Translational Informatics.");
    String href=("http://medicine.unm.edu/informatics/");
    logo_htm+=(HtmUtils.HtmTipper(imghtm,tiphtm,href,200,"white"));
    logo_htm+="</TD><TD>";
    imghtm=("<IMG BORDER=0 SRC=\""+PROXY_PREFIX+CONTEXTPATH+"/images/chemaxon_powered_100px.png\">");
    tiphtm=("JChem from ChemAxon Ltd.");
    href=("http://www.chemaxon.com");
    logo_htm+=(HtmUtils.HtmTipper(imghtm,tiphtm,href,200,"white"));
    logo_htm+="</TD></TR></TABLE>";
    errors.add(logo_htm);
    //errors.add("DEBUG: SERVLETNAME="+SERVLETNAME);

    Calendar calendar=Calendar.getInstance();
    calendar.setTime(new Date());
    DATESTR = String.format("%04d%02d%02d%02d%02d", calendar.get(Calendar.YEAR), calendar.get(Calendar.MONTH)+1, calendar.get(Calendar.DAY_OF_MONTH), calendar.get(Calendar.HOUR_OF_DAY), calendar.get(Calendar.MINUTE));
    Random rand = new Random();
    PREFIX = SERVLETNAME+"."+DATESTR+"."+String.format("%03d",rand.nextInt(1000));

    try {
      LicenseManager.setLicenseFile(CONTEXT.getRealPath("")+"/.chemaxon/license.cxl");
    } catch (Exception e) {
      errors.add("ERROR: ChemAxon LicenseManager error: "+e.getMessage());
    }
    LicenseManager.refresh();
    if (!LicenseManager.isLicensed(LicenseManager.JCHEM))
    {
      errors.add("ERROR: ChemAxon license error; JCHEM is required.");
      return false;
    }

    if (mrequest==null) return true;

    /// Stuff for a run:

    for (Enumeration e=mrequest.getParameterNames(); e.hasMoreElements(); )
    {
      String key=(String)e.nextElement();
      if (mrequest.getParameter(key)!=null)
        params.setVal(key,mrequest.getParameter(key));
    }

    if (params.isChecked("verbose"))
    {
      errors.add("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());
      errors.add("server: "+CONTEXT.getServerInfo());
    }

    String fname="infile";
    File ifile=mrequest.getFile(fname);
    String intxt=params.getVal("intxt").replaceFirst("[\\s]+$","");
    byte[] inbytes=new byte[1024];
    if (ifile!=null)
    {
      FileInputStream fis=new FileInputStream(ifile);
      int asize=inbytes.length;
      int size=0;
      int b;
      while ((b=fis.read())>=0)
      {
        if (size+1>asize)
        {
          asize*=2;
          byte[] tmp=new byte[asize];
          System.arraycopy(inbytes,0,tmp,0,size);
          inbytes=tmp;
        }
        inbytes[size]=(byte)b;
        ++size; 
      }
      byte[] tmp=new byte[size];
      System.arraycopy(inbytes,0,tmp,0,size);
      inbytes=tmp;
    }
    else if (intxt.length()>0)
    {
      if (params.getVal("molfmt").equals("cdx"))
        inbytes=Base64Decoder.decodeToBytes(intxt);
      else
        inbytes=intxt.getBytes("utf-8");
    }
    else
    {
      errors.add("No input data.");
      return false;
    }

    MolImporter molReader=null;
    if (params.getVal("molfmt").equals("automatic"))
    {
      String orig_fname=mrequest.getOriginalFileName(fname);
      String molfmt_auto=MFileFormatUtil.getMostLikelyMolFormat(orig_fname);
      if (orig_fname!=null && molfmt_auto!=null)
        molReader=new MolImporter(new ByteArrayInputStream(inbytes),molfmt_auto);
      else
        molReader=new MolImporter(new ByteArrayInputStream(inbytes));
    }
    else
    {
      String ifmt=params.getVal("molfmt");
      molReader=new MolImporter(new ByteArrayInputStream(inbytes),ifmt);
    }
    String fmt=molReader.getFormat();
    params.setVal("molfmt_auto",fmt);
    if (ifile!=null) ifile.delete();
    MFileFormat mffmt=MFileFormatUtil.getFormat(fmt);
    if (params.isChecked("file2txt"))
    {
      if (mffmt==MFileFormat.CDX) //binary
      {
        intxt=Base64Encoder.encode(inbytes);
        if (params.getVal("molfmt").equals("automatic"))
          params.setVal("molfmt","cdx");
      }
      else
      {
        intxt=new String(inbytes,"utf-8");
      }
      params.setVal("intxt",intxt);
    }
    Molecule m;
    int n_failed=0;
    while (true)
    {
      try { m=molReader.read(); }
      catch (MolFormatException e)
      {
        errors.add("ERROR: MolImporter failed: "+e.getMessage());
        ++n_failed;
        continue;
      }
      if (m==null) break;
      mols.add(m);
    }
    molReader.close();
    if (n_failed>0) errors.add("Errors (unable to read mol): "+n_failed);

    SCRATCHDIR_PATH = CONTEXT.getRealPath(SCRATCHDIR);
    File scratchdir = new File(SCRATCHDIR_PATH);
    if (!(scratchdir.exists() && scratchdir.isDirectory()))
    {
      try {
        scratchdir.mkdir();
        scratchdir.setWritable(true,true);
      }
      catch (Exception e) {
        errors.add("ERROR: SCRATCHDIR creation failed ("+scratchdir+"): "+e.getMessage());
        return false;
      }
    }
    fout=File.createTempFile(PREFIX,"_out.smi",scratchdir);
    fout_scaf=File.createTempFile(PREFIX,"_scaf.smi",scratchdir);

    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String formHtm(MultipartRequest mrequest,HttpServletResponse response)
      throws IOException
  {
    String molfmt_menu="<SELECT NAME=\"molfmt\">\n";
    molfmt_menu+=("<OPTION VALUE=\"automatic\">automatic\n");
    for (String fmt: MFileFormatUtil.getMolfileFormats())
      molfmt_menu+=("<OPTION VALUE=\""+fmt+"\">"+MFileFormatUtil.getFormat(fmt).getDescription()+"\n");
    molfmt_menu+=("</SELECT>\n");
    molfmt_menu=molfmt_menu.replace(params.getVal("molfmt")+"\">",params.getVal("molfmt")+"\" SELECTED>\n");

    String size_menu="<SELECT NAME=\"size\">\n";
    for (String key: sizes_h.keySet())
      size_menu+=("<OPTION VALUE=\""+key+"\">"+key+" - "+sizes_h.get(key)+"x"+sizes_w.get(key)+"\n");
    size_menu+="</SELECT>\n";
    size_menu=size_menu.replace("\""+params.getVal("size")+"\">","\""+params.getVal("size")+"\" SELECTED>\n");

    String runmode_analyze=""; String runmode_1xN=""; String runmode_NxN="";
    if (params.getVal("runmode").equals("1xN")) runmode_1xN="CHECKED";
    else if (params.getVal("runmode").equals("NxN")) runmode_NxN="CHECKED";
    else runmode_analyze="CHECKED";

    String htm=
     ("<FORM NAME=\"mainform\" METHOD=POST ACTION=\""+response.encodeURL(SERVLETNAME)+"\" ENCTYPE=\"multipart/form-data\">\n")
    +("<INPUT TYPE=HIDDEN NAME=\"hscaf\">\n")
    +("<TABLE WIDTH=\"100%\"><TR><TD><H1>"+APPNAME+"</H1></TD>\n")
    +("<TD ALIGN=RIGHT>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"void window.open('"+response.encodeURL(SERVLETNAME)+"?help=TRUE','helpwin','width=600,height=400,scrollbars=1,resizable=1')\"><B>Help</B></BUTTON>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"go_demo(this.form)\"><B>Demo</B></BUTTON>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"window.location.replace('"+response.encodeURL(SERVLETNAME)+"')\"><B>Reset</B></BUTTON>\n")
    +("</TD></TR></TABLE>\n")
    +("<HR>\n")
    +("<TABLE WIDTH=\"100%\" CELLPADDING=5 CELLSPACING=5>\n")
    +("<TR BGCOLOR=\"#CCCCCC\"><TD VALIGN=TOP WIDTH=\"60%\">\n")
    +("format:"+molfmt_menu)
    +("<INPUT TYPE=CHECKBOX NAME=\"file2txt\" VALUE=\"CHECKED\" "+params.getVal("file2txt")+">file2txt<BR>\n")
    +("upload: <INPUT TYPE=\"FILE\" NAME=\"infile\"> ...or paste:")
    +("<BR><TEXTAREA NAME=\"intxt\" WRAP=OFF ROWS=12 COLS=60>"+params.getVal("intxt")+"</TEXTAREA>\n")
    +("</TD>\n")
    +("<TD VALIGN=TOP>\n")
    +("<B>task:</B><BR>\n")
    +("&nbsp;<INPUT TYPE=RADIO NAME=\"runmode\" VALUE=\"analyze\""+runmode_analyze+">analyze mols &amp; group<BR>")
    +("&nbsp;<INPUT TYPE=RADIO NAME=\"runmode\" VALUE=\"1xN\""+runmode_1xN+">max common scaffolds (1xN)<BR>\n")
    +("&nbsp;<INPUT TYPE=RADIO NAME=\"runmode\" VALUE=\"NxN\""+runmode_NxN+">max common scaffolds (NxN)<BR>\n")
    +("<B>output:</B><BR>\n")
    +("&nbsp;<INPUT TYPE=CHECKBOX NAME=\"out_view\" VALUE=\"CHECKED\" "+params.getVal("out_view")+">view table<BR>\n")
    +("&nbsp;<INPUT TYPE=CHECKBOX NAME=\"out_batch\" VALUE=\"CHECKED\" "+params.getVal("out_batch")+">batch<BR>\n")
    +("<B>depictions:</B><BR>\n")
    +("&nbsp;"+size_menu)
    +("&nbsp;<INPUT TYPE=CHECKBOX NAME=\"showarom\" VALUE=\"CHECKED\" "+params.getVal("showarom")+">showarom<BR>\n")
    +("&nbsp; <INPUT TYPE=CHECKBOX NAME=\"show_js\" VALUE=\"CHECKED\""+params.getVal("show_js")+">show attachments<BR>\n")
    +("<BR>\n")
    +("<B>misc:</B><BR>\n")
    +("max mol size (atoms):<INPUT TYPE=TEXT SIZE=4 NAME=\"maxmol\" VALUE=\""+params.getVal("maxmol")+"\"><BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"keep_nitro_attachments\" VALUE=\"CHECKED\""+params.getVal("keep_nitro_attachments")+">keep nitro attachments<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"stereo\" VALUE=\"CHECKED\""+params.getVal("stereo")+">stereo scaffolds<BR>\n")
    +("<INPUT TYPE=CHECKBOX NAME=\"verbose\" VALUE=\"CHECKED\" "+params.getVal("verbose")+">verbose<BR>\n")
    +("</TD></TR></TABLE>\n")
    +("<P>\n")
    +("<CENTER>\n")
    +("<BUTTON TYPE=BUTTON onClick=\"go_hscaf(this.form)\">\n")
    +("<B>Go "+APPNAME+"</B></BUTTON>\n")
    +("</CENTER>\n")
    +("</FORM>\n");
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate scaffold trees for each input molecule.
  */
  private static ArrayList<ScaffoldTree> HScafsGenerate(ArrayList<Molecule> mols)
      throws MolFormatException,MolExportException,IOException
  {
    ArrayList<ScaffoldTree> scaftrees = new ArrayList<ScaffoldTree>();
    ScaffoldSet scafset = new ScaffoldSet("scaffold set");
    MolExporter molWriter=new MolExporter(new FileOutputStream(fout),"smiles:+n");
    int n_mol=0;
    int n_err=0;
    int n_mol_toobig=0;
    int n_total_scaf=0;
    int n_mol_frag=0;
    int maxmol=Integer.parseInt(params.getVal("maxmol"));
    for (n_mol=0;n_mol<mols.size();++n_mol)
    {
      Molecule mol=mols.get(n_mol);
      String molname=mol.getName();
      String verbosetxt=(""+n_mol+". "+molname+"\n");
      verbosetxt+=("\t"+mol.exportToFormat(smifmt)+"\n");
      if (mol.getAtomCount()>maxmol)
      {
        errors.add("Warning: skipping molecule; natoms "+mol.getAtomCount() +">"+ +maxmol);
        ++n_mol_toobig;
        scaftrees.add(null);
        continue;
      }
      if (mol.getFragCount()>1)
      {
        errors.add("Warning: multi-frag mol; analyzing largest frag only: ["+n_mol+"] "+molname);
        ++n_mol_frag;
        mol=hier_scaffolds_utils.LargestPart(mol);
      }

      int ring_count=hier_scaffolds_utils.RawRingsystemCount(mol);
      verbosetxt+=("\traw_ringsys_count: "+ring_count);

      Molecule outmol = new Molecule();
      outmol.fuse(mol.cloneMolecule(),true);

      ScaffoldTree scaftree=null;
      try {
        scaftree = new ScaffoldTree(mol,
		params.isChecked("stereo"),
		params.isChecked("keep_nitro_attachments"),
		scafset);
      }
      catch (Exception e) {
        ++n_err;
        errors.add(e.getMessage());
        scaftrees.add(null);
        continue;
      }

      verbosetxt+=("\tn_scaf="+scaftree.getScaffoldCount());
      verbosetxt+=("\tn_link="+scaftree.getLinkerCount());
      verbosetxt+=("\tn_chain="+scaftree.getSidechainCount()+"\n");

      ArrayList<Long> scaflist = new ArrayList<Long>(); //for this mol
      int n_scaf=0;
      for (Scaffold scaf: scaftree.getScaffolds())
      {
        ++n_scaf;
        scaflist.add(scaf.getID());
        verbosetxt+=("\tscaf: "+n_scaf+". "+scaf.getJsmi()+"\n");
        verbosetxt+=("\t\tID="+(scaf.getID())
          +(" cIDs="+scaf.getChildIDs())
          +(scaf.isRoot()?" (root)":"")
          +(scaf.isLeaf()?" (leaf)":"")
          +"\n");
      }
      n_total_scaf+=n_scaf;
      verbosetxt+=("\t"+scaftree.toString()+"\n");

      Collections.sort(scaflist);
      String scaflist_str="S:";
      for (long id: scaflist) { scaflist_str+=""+id+","; }
      scaflist_str=scaflist_str.replaceFirst(",$","");

      outmol.setName(mol.getName()+" "+scaflist_str);
      outmol.setProperty("HSCAF_SCAFLIST",scaflist_str);

      int n_link=0;
      for (Linker link: scaftree.getLinkers())
      {
        ++n_link;
        verbosetxt+=("\tlinker: "+n_link+". "+link.getJsmi()+"\n");
      }
      int n_chain=0;
      for (Sidechain chain: scaftree.getSidechains())
      {
        ++n_chain;
        verbosetxt+=("\tsidechain: "+n_chain+". "+chain.getJsmi()+"\n");
      }
      molWriter.write(outmol);
      scaftrees.add(scaftree);
      if (params.isChecked("verbose"))
        errors.add("<PRE>"+verbosetxt+"</PRE>");
    }
    molWriter.close();
    MolExporter molWriter_scaf=new MolExporter(new FileOutputStream(fout_scaf),"smiles:+n");

    int n_scaf_unique=0;
    for (int scaf_id=1;scaf_id<=scafset.getCount();++scaf_id)
    {
      Scaffold scaf=scafset.getScaffoldByID(scaf_id);
      scaf.decompress();
      scaf.setName(""+scaf_id+" "+scaf.subTreeAsString());
      molWriter_scaf.write(scaf);
      ++n_scaf_unique;
    }
    molWriter_scaf.close();
    errors.add("total mols: "+n_mol);
    errors.add("errors: "+n_err);
    errors.add("total scaffolds found: "+n_total_scaf);
    errors.add("total unique scaffolds found: "+n_scaf_unique);
    errors.add("oversized mols rejected: "+n_mol_toobig);
    errors.add("Multi-fragment mols (largest part analyzed only): "+n_mol_frag);
    return scaftrees;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static ArrayList<Scaffold> MCScafs_1xN(ArrayList<ScaffoldTree> scaftrees)
  {
    ArrayList<Scaffold> mcscafs = new ArrayList<Scaffold>();
    ScaffoldTree scaftreeQ = scaftrees.get(0);
    int n_mcscaf=0;
    for (int i=1;i<scaftrees.size();++i)
    {
      Scaffold mcscaf=null;
      if (scaftrees.get(i)!=null) mcscaf=hier_scaffolds_utils.MaxCommonScaffold(scaftreeQ,scaftrees.get(i));
      mcscafs.add(mcscaf);
      if (mcscaf!=null) ++n_mcscaf;
    }
    errors.add("max common scafs found: "+n_mcscaf);
    return mcscafs;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static ArrayList<ArrayList<Scaffold> > MCScafs_NxN(ArrayList<ScaffoldTree> scaftrees)
  {
    ArrayList<ArrayList<Scaffold> > mcscaf_matrix = new ArrayList<ArrayList<Scaffold> >();
    int n_mcscaf=0;
    for (int i=0;i<scaftrees.size();++i)
    {
      ArrayList<Scaffold> mcscafs = new ArrayList<Scaffold>();
      for (int j=0;j<i;++j) mcscafs.add(null);
      for (int j=i;j<scaftrees.size();++j)
      {
        Scaffold mcscaf=null;
        if (scaftrees.get(i)!=null && scaftrees.get(j)!=null)
          mcscaf=hier_scaffolds_utils.MaxCommonScaffold(scaftrees.get(i),scaftrees.get(j));
        mcscafs.add(mcscaf);
        if (mcscaf!=null) ++n_mcscaf;
      }
      mcscaf_matrix.add(mcscafs);
    }
    errors.add("max common scafs found: "+n_mcscaf);
    return mcscaf_matrix;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate output for runmode "analyze".
  */
  private static void HScafsResults(HttpServletResponse response,ArrayList<ScaffoldTree> scaftrees)
      throws IOException,FileNotFoundException
  {
    int n_mol=0;
    int w=sizes_w.get(params.getVal("size"));
    int h=sizes_h.get(params.getVal("size"));
    String depopts=("mode=cow&imgfmt=png");
    depopts+=("&smi_is_cx=true");
    if (params.isChecked("showarom")) depopts+=("&arom_gen=true");
    else depopts+=("&kekule=true");
    String thtm=("<TABLE BORDER>\n");
    thtm+=("<TR><TH></TH><TH>mol</TH><TH>scaffolds</TH><TH>linkers</TH><TH>sidechains</TH></TR>\n");
    int n_total_scaf=0;
    for (ScaffoldTree scaftree: scaftrees)
    {
      if (scaftree==null) continue;
      n_total_scaf+=scaftree.getScaffoldCount();
      ++n_mol;
      String molname=scaftree.inmol.getName();
      if (params.isChecked("verbose"))
        errors.add("\tn_scaf:"+scaftree.getScaffoldCount()+" n_link:"+scaftree.getLinkerCount()+" n_chain:"+scaftree.getSidechainCount()); 
      thtm+="<TR>\n";
      thtm+=("<TD VALIGN=TOP ALIGN=RIGHT>"+n_mol+".</TD>\n");
      String smiles=scaftree.inmol.exportToFormat(smifmt_dep);
      String imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts,h-30,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
      thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm+"<BR>\n");
      thtm+=(molname.length()<35?molname:(molname.substring(0,31)+"...")+"</TD>\n");
      if (scaftree.getScaffoldCount()>0) {
        smiles=scaftree.getSmiForScaffoldsGroupmol(params.isChecked("show_js"));
        imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts,h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER>"+imghtm+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER>~</TD>\n");
      if (scaftree.getLinkerCount()>0) {
        smiles=scaftree.getSmiForLinkersGroupmol(params.isChecked("show_js"));
        imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts,h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER>"+imghtm+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER>~</TD>\n");
      if (scaftree.getSidechainCount()>0) {
        smiles=scaftree.getSmiForSidechainsGroupmol(params.isChecked("show_js"));
        imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts,h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER>"+imghtm+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER>~</TD>\n");
      thtm+="</TR>\n";
      if (n_mol==N_MAX) break;
    }
    thtm+=("</TABLE>");
    outputs.add("<H1>output:</H1>");
    outputs.add("total mols: "+scaftrees.size());
    outputs.add("total mols analyzed: "+n_mol);
    outputs.add("total scaffolds found: "+n_total_scaf);
    if (params.isChecked("out_view"))
      outputs.add("<CENTER>"+thtm+"</CENTER>");
    if (params.isChecked("out_batch"))
    {
      String fname=(SERVLETNAME+"_scaf.smi");
      outputs.add("&nbsp;"+
        "<FORM METHOD=\"POST\" ACTION=\""+response.encodeURL(SERVLETNAME)+"\">\n"+
        "<INPUT TYPE=HIDDEN NAME=\"downloadfile\" VALUE=\""+fout_scaf.getAbsolutePath()+"\">\n"+
        "<INPUT TYPE=HIDDEN NAME=\"fname\" VALUE=\""+fname+"\">\n"+
        "<BUTTON TYPE=BUTTON onClick=\"this.form.submit()\">"+
        "download "+fname+" ("+file_utils.NiceBytes(fout_scaf.length())+")</BUTTON></FORM>\n");
      fname=(SERVLETNAME+"_out.smi");
      outputs.add("&nbsp;"+
        "<FORM METHOD=\"POST\" ACTION=\""+response.encodeURL(SERVLETNAME)+"\">\n"+
        "<INPUT TYPE=HIDDEN NAME=\"downloadfile\" VALUE=\""+fout.getAbsolutePath()+"\">\n"+
        "<INPUT TYPE=HIDDEN NAME=\"fname\" VALUE=\""+fname+"\">\n"+
        "<BUTTON TYPE=BUTTON onClick=\"this.form.submit()\">"+
        "download "+fname+" ("+file_utils.NiceBytes(fout.length())+")</BUTTON></FORM>\n");
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate output for runmode "1xN".
  */
  private static void HScafsResults_1xN(HttpServletResponse response,ArrayList<ScaffoldTree> scaftrees,ArrayList<Scaffold> mcscafs)
      throws IOException,FileNotFoundException
  {
    int n_mol=0;
    int w=sizes_w.get(params.getVal("size"));
    int h=sizes_h.get(params.getVal("size"));
    String depopts=("mode=cow&imgfmt=png");
    depopts+=("&smi_is_cx=true");
    if (params.isChecked("showarom")) depopts+=("&arom_gen=true");
    else depopts+=("&kekule=true");

    String thtm=("<TABLE BORDER>\n");
    thtm+=("<TR><TH></TH><TH>query</TH><TH>max common scaffold</TH><TH>mol</TH></TR>\n");
    String molnameQ=scaftrees.get(0).inmol.getName();
    String smilesQ=scaftrees.get(0).inmol.exportToFormat(smifmt_dep);
    String imghtmQ=HtmUtils.Smi2ImgHtm(smilesQ,depopts,h-30,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
    int n_total_scaf=0;
    for (int i=1;i<scaftrees.size();++i)
    {
      if (scaftrees.get(i)==null) continue;
      n_total_scaf+=scaftrees.get(i).getScaffoldCount();
      ++n_mol;
      float sim=hier_scaffolds_utils.CommonScaffoldTanimoto(scaftrees.get(0),scaftrees.get(i),mcscafs.get(i-1));
      if (params.isChecked("verbose"))
        errors.add(String.format("scaffold Tanimoto (%d) = %.2f",i,sim));
      thtm+="<TR>\n";
      thtm+=("<TD VALIGN=TOP ALIGN=RIGHT>"+i+".</TD>\n");
      thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtmQ+"<BR>\n");
      thtm+=(molnameQ.length()<35?molnameQ:(molnameQ.substring(0,31)+"...")+"</TD>\n");
      if (mcscafs.size()>(i-1) && mcscafs.get(i-1)!=null)
      {
        String smiles_mcscaf=mcscafs.get(i-1).getSmi();
        String imghtm_mcscaf=HtmUtils.Smi2ImgHtm(smiles_mcscaf,depopts,h-30,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm_mcscaf+"<BR>\n");
        thtm+=(String.format("hscaf_Tanimoto = %.2f",sim)+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER VALIGN=MIDDLE>~</TD>\n");
      String smiles=scaftrees.get(i).inmol.exportToFormat(smifmt_dep);
      String molname=scaftrees.get(i).inmol.getName();
      String imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts,h-30,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
      thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm+"<BR>\n");
      thtm+=(molname.length()<35?molname:(molname.substring(0,31)+"...")+"</TD>\n");
      thtm+="</TR>\n";
      if (n_mol==N_MAX) break;
    }
    thtm+=("</TABLE>");
    outputs.add("<H1>output:</H1>");
    outputs.add("total mols: "+scaftrees.size());
    outputs.add("total mols analyzed: "+n_mol);
    outputs.add("total scaffolds found: "+n_total_scaf);
    if (params.isChecked("out_view"))
      outputs.add("<CENTER>"+thtm+"</CENTER>");
    if (params.isChecked("out_batch"))
      outputs.add("NOTE: 1xN batch output not yet supported.");
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate output for runmode "NxN".
  */
  private static void HScafsResults_NxN(HttpServletResponse response,ArrayList<ScaffoldTree> scaftrees,ArrayList<ArrayList<Scaffold> > mcscaf_matrix)
      throws IOException,FileNotFoundException
  {
    // This is our convention; Apache proxies the 8080 port via /tomcat.
    int n_mol=0;
    int w=sizes_w.get(params.getVal("size"))/2;
    int h=sizes_h.get(params.getVal("size"))/2;
    String depopts=("mode=cow&imgfmt=png");
    depopts+=("&smi_is_cx=true");
    if (params.isChecked("showarom")) depopts+=("&arom_gen=true");
    else depopts+=("&kekule=true");

    String thtm=("<TABLE BORDER>\n<TR><TD></TD>\n");
    String smiles=null;
    String imghtm=null;
    int n_total_scaf=0;
    for (int i=1;i<scaftrees.size();++i)
    {
      if (scaftrees.get(i)!=null)
      {
        n_total_scaf+=scaftrees.get(i).getScaffoldCount();
        ++n_mol;
        smiles=scaftrees.get(i).inmol.exportToFormat(smifmt_dep);
        imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts+"&bgcolor=%23DDDDDD",h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER VALIGN=MIDDLE WIDTH="+w+">~</TD>\n");
    }
    thtm+=("</TR>\n");
    for (int i=1;i<scaftrees.size();++i)
    {
      thtm+=("<TR>\n");
      if (scaftrees.get(i)!=null)
      {
        smiles=scaftrees.get(i).inmol.exportToFormat(smifmt_dep);
        imghtm=HtmUtils.Smi2ImgHtm(smiles,depopts+"&bgcolor=%23DDDDDD",h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
        thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm+"</TD>\n");
      }
      else
        thtm+=("<TD ALIGN=CENTER VALIGN=MIDDLE HEIGHT="+h+">~</TD>\n");
      for (int j=1;j<scaftrees.size();++j)
      {
        Scaffold mcscaf=null;
        if (scaftrees.get(i)!=null && scaftrees.get(j)!=null)
          mcscaf=((i<j)?mcscaf_matrix.get(i).get(j):mcscaf_matrix.get(j).get(i));
        if (mcscaf!=null)
        {
          imghtm=HtmUtils.Smi2ImgHtm(mcscaf.getSmi(),depopts,h,w,MOL2IMG_SERVLETURL,true,4,"go_zoom_smi2img");
          thtm+=("<TD ALIGN=CENTER VALIGN=TOP>"+imghtm+"</TD>\n");
        }
        else
          thtm+=("<TD ALIGN=CENTER VALIGN=MIDDLE>~</TD>\n");
      }
      thtm+=("</TR>\n");
    }
    thtm+=("</TABLE>");
    outputs.add("<H1>output:</H1>");
    outputs.add("total mols: "+scaftrees.size());
    outputs.add("total mols analyzed: "+n_mol);
    outputs.add("total scaffolds found: "+n_total_scaf);
    if (params.isChecked("out_view"))
      outputs.add("<CENTER>"+thtm+"</CENTER>");
    if (params.isChecked("out_batch"))
      outputs.add("NOTE: NxN batch output not yet supported.");
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String JavaScript()
      throws IOException
  {
    String js="var demotxt='';";
    if (DEMOSMIFILE!=null) {
      String absoluteDiskPath = CONTEXT.getRealPath(DEMOSMIFILE);
      BufferedReader buff=new BufferedReader(new FileReader(absoluteDiskPath ));
      String line=null;
      String startdate=null;
      while ((line=buff.readLine())!=null)
        js+=("demotxt+='"+line+"\\n';\n");
      buff.close();
    }
    js+=(
"function go_init(form)"+
"{\n"+
"  form.file2txt.checked=true;\n"+
"  form.intxt.value='';\n"+
"  var i;\n"+
"  for (i=0;i<form.molfmt.length;++i)\n"+
"    if (form.molfmt.options[i].value=='automatic')\n"+
"      form.molfmt.options[i].selected=true;\n"+
"  for (i=0;i<form.size.length;++i)\n"+
"    if (form.size.options[i].value=='s')\n"+
"      form.size.options[i].selected=true;\n"+
"  for (i=0;i<form.runmode.length;++i)\n"+ //radio
"    if (form.runmode[i].value=='analyze')\n"+
"      form.runmode[i].checked=true;\n"+
"  form.show_js.checked=false;\n"+
"  form.keep_nitro_attachments.checked=false;\n"+
"  form.stereo.checked=false;\n"+
"  form.showarom.checked=false;\n"+
"  form.maxmol.value='50';\n"+
"  form.out_view.checked=true;\n"+
"  form.out_batch.checked=true;\n"+
"  form.verbose.checked=false;\n"+
"}\n"+
"function checkform(form)\n"+
"{\n"+
"  if (!form.intxt.value && !form.infile.value) {\n"+
"    alert('ERROR: No input specified');\n"+
"    return 0;\n"+
"  }\n"+
"  return 1;\n"+
"}\n"+
"function go_demo(form) {\n"+
"  go_init(form);\n"+
"  form.intxt.value=demotxt;\n"+
"  form.hscaf.value='TRUE'\n"+
"  form.submit()\n"+
"}\n"+
"function go_hscaf(form)\n"+
"{\n"+
"  if (!checkform(form)) return;\n"+
"  form.hscaf.value='TRUE'\n"+
"  form.submit()\n"+
"}\n");
    return(js);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static String HelpHtm()
  {
    return (
    "<H2>"+APPNAME+" Help</H2>\n"+
    "Based on the HierS hierarchical scaffold algorithm (Wilkens, 2005).\n"+
    "Built with package edu.unm.health.biocomp.hscaf\n"+
    "(GitHub repo: <a href=\"https://github.com/unmtransinfo/unm_biocomp_hscaf\">unm_biocomp_hscaf</a>).\n"+
    "<P>\n"+
    "N_MAX = "+N_MAX+"\n"+
    "<P>\n"+
    "Thanks to <a href=\"http://www.chemaxon.com\">ChemAxon</a> for the use of JChem in this application.\n"+
    "<P>\n"+
    "Author/support: Jeremy Yang\n"
    );
  }
  /////////////////////////////////////////////////////////////////////////////
  public void init(ServletConfig conf) throws ServletException
  {
    super.init(conf);
    CONTEXT=getServletContext();
    CONTEXTPATH=CONTEXT.getContextPath();
    //CONFIG=conf;
    // read servlet parameters (from web.xml):
    try { APPNAME=conf.getInitParameter("APPNAME"); }
    catch (Exception e) { APPNAME=this.getServletName(); }
    UPLOADDIR=conf.getInitParameter("UPLOADDIR");
    if (UPLOADDIR==null)
      throw new ServletException("Please supply UPLOADDIR parameter");
    SCRATCHDIR=conf.getInitParameter("SCRATCHDIR");
    if (SCRATCHDIR==null) SCRATCHDIR="/tmp";
    DEMOSMIFILE=conf.getInitParameter("DEMOSMIFILE");
    try { N_MAX=Integer.parseInt(conf.getInitParameter("N_MAX")); }
    catch (Exception e) { N_MAX=100; }
    PROXY_PREFIX=((conf.getInitParameter("PROXY_PREFIX")!=null)?conf.getInitParameter("PROXY_PREFIX"):"");
  }
  /////////////////////////////////////////////////////////////////////////////
  public void doGet(HttpServletRequest request,HttpServletResponse response)
      throws IOException, ServletException
  {
    doPost(request,response);
  }
}

