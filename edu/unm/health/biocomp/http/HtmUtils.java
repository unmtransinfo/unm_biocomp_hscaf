package edu.unm.health.biocomp.http;

import java.io.*;
import java.awt.Color;
import java.util.*;
import java.util.zip.*;
import java.util.regex.*;
import java.text.*;
import java.net.URLEncoder;
import javax.servlet.*;
import javax.servlet.http.*;
import org.apache.commons.codec.binary.Base64;

/**	Static methods for web apps.
	Molfiles etc. are gzipped and base64ed by calling program.
	@author Jeremy J Yang
*/
public class HtmUtils
{
  private HtmUtils() {} //disable default constructor
  /////////////////////////////////////////////////////////////////////////////
  /**	HTML for a image from smiles.
  */
  public static String smi2ImgHtm(String smi,String opts,int h,int w,
    String mol2img_url,boolean zoomable,int zoomfactor,String jsfunc)
	throws UnsupportedEncodingException
  {
    if (opts.length()==0) opts="mode=cow";
    String smicode=URLEncoder.encode(smi,"UTF-8");
    String imgurl=(mol2img_url+"?"+opts+"&h="+h+"&w="+w+"&smiles="+smicode);
    String imghtm=("<IMG SRC='"+imgurl+"' BORDER=0>");
    if (zoomable)
    {
      imghtm=String.format("<A HREF=\"javascript:void(0)\" onClick=\"javascript:%s('%s','%s','%s',%d,%d)\">%s</A>",jsfunc,mol2img_url,smicode,opts+"&maxscale=0",w*zoomfactor,h*zoomfactor,imghtm);
    }
    return imghtm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	HTML for a image from molfile.
  */
  public static String mdl2ImgHtm(String mdl_gz_b64,String opts,int h,int w,
    String mol2img_url,boolean zoomable,int zoomfactor,String jsfunc)
	throws UnsupportedEncodingException,IOException
  {
    if (opts.length()==0) opts="mode=cow";
    mdl_gz_b64=URLEncoder.encode(mdl_gz_b64,"UTF-8");
    String imgurl=(mol2img_url+"?"+opts+"&h="+h+"&w="+w+"&mdlcode="+mdl_gz_b64);
    String imghtm=("<IMG SRC='"+imgurl+"' BORDER=0>");
    if (zoomable)
    {
      imghtm=String.format("<A HREF=\"javascript:void(0)\" onClick=\"javascript:%s('%s','%s','%s',%d,%d)\">%s</A>",jsfunc,mol2img_url,mdl_gz_b64,opts+"&maxscale=0",w*zoomfactor,h*zoomfactor,imghtm);
    }
    return imghtm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	HTML for a image from mrvfile (ChemAxon mrv format).
  */
  public static String mrv2ImgHtm(String mrv_gz_b64,ArrayList<Color> atomcolors,
	String opts,int h,int w,
	String mol2img_url,boolean zoomable,int zoomfactor,String jsfunc)
	throws UnsupportedEncodingException
  {
    if (opts.length()==0) opts="mode=cow";
    mrv_gz_b64=URLEncoder.encode(mrv_gz_b64,"UTF-8");
    if (atomcolors!=null)
    {
      opts+=("&atomcolors=true");
      for (int i=0;i<atomcolors.size();++i)
      {
        Color c=atomcolors.get(i);
        opts+=("&color"+i+"="+String.format("%02X%02X%02X,",c.getRed(),c.getGreen(),c.getBlue()));
      }
    }
    String imgurl=(mol2img_url+"?"+opts+"&h="+h+"&w="+w+"&mrvcode="+mrv_gz_b64);
    String imghtm=("<IMG SRC='"+imgurl+"' BORDER=0>");
    if (zoomable)
    {
      imghtm=String.format("<A HREF=\"javascript:void(0)\" onClick=\"javascript:%s('%s','%s','%s',%d,%d)\">%s</A>",jsfunc,mol2img_url,mrv_gz_b64,opts+"&maxscale=0",w*zoomfactor,h*zoomfactor,imghtm);
    }
    return imghtm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Header HTML with included JavaScript.
  	(For use with JS + CSS via Apache.)
  */
  public static String headerHtm(String title,String js,String color)
  {
    String htm=(
    "<HTML>\n"+
    "<HEAD><TITLE>"+title+"</TITLE>\n"+
    "<SCRIPT>\n"+js+"\n</SCRIPT>\n"+
    "<SCRIPT SRC=\"/js/Mol2Img.js\"></SCRIPT>\n"+
    "<SCRIPT SRC=\"/marvin/marvin.js\"></SCRIPT>\n"+
    "<SCRIPT SRC=\"/js/ddtip.js\"></SCRIPT>\n"+
    "<LINK REL=\"stylesheet\" type=\"text/css\" HREF=\"/css/biocomp.css\" />\n"+
    "</HEAD>\n"+
    "<BODY BGCOLOR=\""+color+"\">\n"+
    "<DIV ID=\"ddtooltip\"></DIV>\n"+
    "<SCRIPT TYPE=\"text/javascript\">ddtip_init()</SCRIPT>\n");
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Header HTML with included JavaScript.

	New: No more hard-coding of specific JS or CSS files in 
	this function.  Use jsincludes and cssincludes for all.
	Also no JS code here.

	If jsinclude|cssinclude bare filename (no '/'),
	prepend with conventional URL prefix, e.g. "/tomcat/biocomp/js".
	If path specified, then include unmodified (could be via
	Tomcat, httpd, or external URL).
  */
  public static String headerHtm(String title,
	List<String> jsincludes,
	List<String> cssincludes,
	String js,
	String color,
	HttpServletRequest request)
  {
    String htm=(
    "<HTML>\n"+
    "<HEAD><TITLE>"+title+"</TITLE>\n");
    //htm+=("<SCRIPT SRC=\"/tomcat"+request.getContextPath()+"/js/Mol2Img.js\"></SCRIPT>\n");
    //htm+=("<SCRIPT SRC=\"/tomcat"+request.getContextPath()+"/js/ddtip.js\"></SCRIPT>\n");
    //htm+=("<SCRIPT SRC=\"/marvin/marvin.js\"></SCRIPT>\n");
    if (cssincludes!=null)
    {
      for (String cssinclude: cssincludes)
      {
        if (cssinclude.contains("/"))
          htm+=("<LINK REL=\"stylesheet\" HREF=\""+cssinclude+"\" />\n");
        else
          htm+=("<LINK REL=\"stylesheet\" type=\"text/css\" HREF=\"/tomcat"+request.getContextPath()+"/css/"+cssinclude+"\" />\n");
      }
    }
    if (jsincludes!=null)
    {
      for (String jsinclude: jsincludes)
      {
        if (jsinclude.contains("/"))
          htm+=("<SCRIPT SRC=\""+jsinclude+"\"></SCRIPT>\n");
        else
          htm+=("<SCRIPT SRC=\"/tomcat"+request.getContextPath()+"/js/"+jsinclude+"\"></SCRIPT>\n");
      }
    }
    htm+=("<SCRIPT>\n"+js+"\n</SCRIPT>\n");
    htm+=("<LINK REL=\"stylesheet\" type=\"text/css\" HREF=\"/tomcat"+request.getContextPath()+"/css/biocomp.css\" />\n");
    htm+=("</HEAD>\n"+"<BODY BGCOLOR=\""+color+"\">\n");
    htm+=("<DIV ID=\"ddtooltip\"></DIV>\n"+"<SCRIPT TYPE=\"text/javascript\">ddtip_init()</SCRIPT>\n");
    return htm;
  }
  /**	Legacy entry point.
  */
  public static String headerHtm(String title,String js,String color,
	HttpServletRequest request)
  {
    return headerHtm(title,null,null,js,color,request);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Output lines to HTML delimited by line breaks.
  */
  public static String outputHtm(ArrayList<String> outputs)
  {
    String htm="";
    htm+=("<HR>\n");
    if (outputs.size()>0)
    {
      for (String output:outputs)
      {
        htm+=(output+"<BR>\n");
      }
    }
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Output lines to HTML delimited by line breaks; terminate HTML.
  */
  public static String footerHtm(ArrayList<String> errors, boolean showErrors)
  {
    String htm="";
    htm+=("<HR>\n");
    if (showErrors && errors.size()>0)
    {
      for (String error:errors)
      {
        htm+=(error+"<BR>\n");
      }
    }
    htm+=("</BODY>\n");
    htm+=("</HTML>\n");
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Human readable file size.
  */
  public static String niceBytes(long bytes)
  {
    if (bytes<1.0e3)
      return String.format("%d bytes",bytes);
    else if (bytes<1.0e6)
      return String.format("%.1f KB",((float)bytes/1.0e3));
    else if (bytes<1.0e9)
      return String.format("%.1f MB",((float)bytes/1.0e6));
    else if (bytes<1.0e12)
      return String.format("%.1f GB",((double)bytes/1.0e9));
    else
      return String.format("%.1f TB",((double)bytes/1.0e12));
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Human readable time interval.
  */
  public static String timeDeltaStr(Date t_i,Date t_f)
  {
    long t_d=t_f.getTime()-t_i.getTime();
    int t_d_m = (int)(t_d/60000L);
    int t_d_s = (int)((t_d/1000L)%60L);
    if (t_d_m==0)
      return (""+t_d_s+"s");
    else
      return (""+t_d_m+"m "+t_d_s+"s");
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Download specified string to browser with specified file name.
  */
  public static boolean downloadString(HttpServletResponse response,
        ServletOutputStream ostream,String str,String fname)
        throws IOException
  {
    //byte[] outbytes=Base64Decoder.decodeToBytes(str) ;
    byte[] outbytes=Base64.decodeBase64(str);
    response.setContentType("application/x-savefile"); //which is better?
    //response.setContentType("application/octet-stream"); //which is better?
    response.setContentLength(outbytes.length);
    response.setDateHeader("Expires",0);
    response.setHeader("Pragma","No-cache");
    response.setHeader("Cache-Control","no-cache");
    response.setHeader("Content-Disposition","attachment; filename="+fname);
    ostream.write(outbytes);
    ostream.flush();
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Download specified file to browser with specified file name.
  */
  public static boolean downloadFile(HttpServletResponse response,
        ServletOutputStream ostream,String fpath,String fname)
        throws IOException
  {
    File f = new File(fpath);
    if (!f.exists()) return false;
    long size = f.length();
    response.setContentType("application/x-savefile"); //which is better?
    //response.setContentType("application/octet-stream"); //which is better?
    response.setContentLength((int)size);
    response.setDateHeader("Expires",0);
    response.setHeader("Pragma","No-cache");
    response.setHeader("Cache-Control","no-cache");
    response.setHeader("Content-Disposition","attachment; filename="+fname);
    FileInputStream fis=new FileInputStream(f);
    int b;
    while ((b=fis.read())>=0) ostream.write(b);
    ostream.flush();
    return true;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate HTML with tooltip. Requires corresponding JavaScript.
  */
  public static String htmTipper(String htm_in,String tiphtm,String href,
	Integer width,String color)
  {
    int w=160;
    String c="yellow";
    if (width!=null) w=width;
    if (color!=null) c=color;
    tiphtm=tiphtm.replace("'","\\\\\\'");
    tiphtm=tiphtm.replace("\"","\\\"");
    String htm=("<A onMouseOver=\"dd_tip('"+tiphtm+"','"+c+"',"+w+")\"");
    if (href!=null && !href.isEmpty())
    {
      htm+=(" HREF=\""+href+"\"");
    }
    htm+=(" onMouseOut=\"dd_hidetip()\">"+htm_in+"</A>");
    return htm;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Generate HTML with tooltip. Requires corresponding JavaScript.
  */
  public static String htmTipper(String htm_in,String tiphtm,Integer width,
	String color)
  {
    return htmTipper(htm_in,tiphtm,null,width,color);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Delete temporary files; if non-existent no worries..
  */
  public static boolean DelFiles(ArrayList<String> files)
  {
    boolean ok=true;
    for (String file: files)
    {
      File f = new File(file);
      if (f.exists())
      {
        if (!f.delete()) {
          ok=false;
          //System.err.println("DEBUG: failed to delete: "+file);
        }
      }
    }
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Deprecated .
  */
  public static boolean PurgeScratchDirs(String [] dirs,int retire_sec,boolean verbose,
	String delim,HttpServlet servlet)
  {
    return PurgeScratchDirs(Arrays.asList(dirs),retire_sec,verbose,delim,servlet);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	For purging temporary files created by
	webapps (HttpServlet servlets), where the names indicate creation time.
	Note that timestamp is expected to be in format YYYYMMDDHHMMSS.
	@param dirs Directories to purge.
	@param retire_sec How old files must be to purge.
	@param verbose If true, log result to HttpServlet context.
	@param delim Filename delimiter defining timestamp.
	@param servlet HttpServlet calling object.
  */
  public static boolean PurgeScratchDirs(List<String> dirs,int retire_sec,boolean verbose,
	String delim,HttpServlet servlet)
  {
    ServletContext context=servlet.getServletContext();
    String servletname=servlet.getServletName();
    boolean ok=true;
    Calendar calendar=Calendar.getInstance();
    int ndeleted=0;
    int ntotal=0;
    for (String dir: dirs)
    {
      File d = new File(dir);
      File[] files = d.listFiles();
      if (files==null) continue;
      for (File f: files)
      {
        ++ntotal;
        String fname=f.getName();
        String[] fields=Pattern.compile("\\.").split(fname);
        if (fields.length<3) continue;
        String datestr=fields[1];
        Date date_file=null;
        if (!Pattern.matches("\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d$",datestr))
        {
          if (verbose)
            context.log(String.format("%s: ERROR: bad date format: \"%s\"",servletname,datestr));
          ok=false;
          date_file=new Date(f.lastModified());
        }
        else
        {
          calendar.set(Integer.parseInt(datestr.substring(0,4)),
               Integer.parseInt(datestr.substring(4,6))-1,
               Integer.parseInt(datestr.substring(6,8)),
               Integer.parseInt(datestr.substring(8,10)),
               Integer.parseInt(datestr.substring(10,12)),
               Integer.parseInt(datestr.substring(12,14)));
          date_file=calendar.getTime();
        }
        Date date_now = new Date();
        long delta_s = (date_now.getTime()-date_file.getTime())/1000;
        if (delta_s>retire_sec)
        {
          if (f.canWrite())
          {
            if (f.delete())
              ++ndeleted;
            else
              context.log(String.format("%s: ERROR: failed to delete \"%s\"",servletname,fname));
          }
          else
          {
            if (verbose)
              context.log(String.format("%s: ERROR: cannot delete \"%s\"",servletname,fname));
          }
        }
      }
    }
    if (verbose)
      context.log(String.format("%s: %d/%d scratch files deleted",servletname,ndeleted,ntotal));
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  public static String htmlEscape(String aText) { 
    if (aText == null) { 
      return ""; 
    } 
    StringBuffer escapedText = new StringBuffer(); 
    for (int i=0; i<aText.length(); i++) { 
      char ch = aText.charAt(i); 
    if (ch == '\'') 
      escapedText.append("&#39;"); 
    else if (ch == '\"') 
      escapedText.append("&#34;"); 
    else if (ch == '<') 
      escapedText.append("&lt;"); 
    else if (ch == '>') 
      escapedText.append("&gt;"); 
    else if (ch == '&') 
      escapedText.append("&amp;"); 
    else 
      escapedText.append(ch); 
    } 
    String result = escapedText.toString(); 
    return result; 
  } 
}
