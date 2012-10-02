package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
import java.sql.*;
import javax.sql.*;
import javax.naming.*;

import com.mysql.jdbc.Driver;


/**	Static utility methods for MySQL databases&#46;
	<br />
	Uses MySQL JDBC driver (com.mysql.jdbc.Driver)&#46;
	<br />
	@author Jeremy J Yang
	@see java.sql.DriverManager
	@see com.mysql.jdbc.Driver
*/
public class mysql_utils
{
  /////////////////////////////////////////////////////////////////////////////
  /**	Return text with server status information&#46;
  */
  public static String serverStatusTxt(Connection dbcon)
  {
    String txt="";
    try {
      ResultSet rset=executeSql(dbcon,"SHOW SERVER_VERSION");
      if (rset.next())
        txt+=("MySQL server version: "+rset.getString(1));
      rset.getStatement().close();
    }
    catch (Exception e)
    {
      txt+=("error: "+e.getMessage());
    }
    return txt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement&#46; Normally use this function&#46;
  */
  public static ResultSet executeSql(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement(ResultSet.TYPE_FORWARD_ONLY,ResultSet.CONCUR_READ_ONLY);
    ResultSet rset=stmt.executeQuery(sql);
    return rset;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement&#46; Use this function for scrollable ResultSet&#46;
  */
  public static ResultSet executeSqlScrollable(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,ResultSet.CONCUR_READ_ONLY);
    ResultSet rset=stmt.executeQuery(sql);
    return rset;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Return MySQL connection&#46;
  */
  public static Connection dbConnect(String dbhost,Integer dbport,String dbid,String dbusr,String dbpw)
    throws SQLException
  {
    DriverManager.registerDriver(new com.mysql.jdbc.Driver());
    Connection dbcon=DriverManager.getConnection("jdbc:mysql://"+dbhost+":"+dbport+"/"+dbid,dbusr,dbpw);
    return dbcon;
  }
  public static Connection dbConnectBard()
    throws SQLException,javax.naming.NamingException
  {
    javax.naming.Context ctx_init = new javax.naming.InitialContext();
    javax.naming.Context ctx = (javax.naming.Context) ctx_init.lookup("java:comp/env");
    javax.sql.DataSource ds = (javax.sql.DataSource) ctx.lookup("jdbc/bard");
    Connection dbcon=ds.getConnection();

    return dbcon;
  }
  /////////////////////////////////////////////////////////////////////////////
  private static void help(String msg)
  {
    System.err.println(msg+"\n"
      +"mysql_utils - mysql utilities\n"
      +"usage: mysql_utils [options]\n"
      +"  required:\n"
      +"    -dbhost <HOST> ... \n"
      +"    -dbname <NAME> ... \n"
      +"    -dbusr <USER> ... \n"
      +"    -dbpw <PW> ... \n"
      +"  mode (one of):\n"
      +"    -describe .............. describe database\n"
      +"    -v ..................... verbose\n"
      +"    -vv .................... very verbose\n"
      +"    -h ..................... this help\n");
    System.exit(1);
  }
  private static int verbose=0;
  private static String dbhost="localhost";
  private static String dbname="bard";
  private static String dbusr="bard";
  private static String dbpw="foobar";
  private static Boolean describe=false;
  private static Boolean bard=false;

  /////////////////////////////////////////////////////////////////////////////
  private static void parseCommand(String args[])
  {
    if (args.length==0) help("");
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-dbhost")) dbhost=args[++i];
      else if (args[i].equals("-dbname")) dbname=args[++i];
      else if (args[i].equals("-dbusr")) dbusr=args[++i];
      else if (args[i].equals("-dbpw")) dbpw=args[++i];
      else if (args[i].equals("-describe")) describe=true;
      else if (args[i].equals("-bard")) bard=true;
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-h")) help("");
      else help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Simple test connection&#46;
  */
  public static void main(String[] args)
	throws IOException,SQLException
  {
    parseCommand(args);

    Connection dbcon=null;

    if (bard)
    {
      try { dbcon=dbConnectBard(); }
      catch (Exception e) { help("MySQL BARD connection failed:"+e.getMessage()); }
      System.exit(0);
    }

    try { dbcon=dbConnect(dbhost,3306,dbname,dbusr,dbpw); }
    catch (SQLException e) { help("MySQL connection failed:"+e.getMessage()); }

    if (dbcon!=null)
      System.err.println("MySQL connection ok: jdbc:mysql://"+dbhost+":/"+dbname);

    if (describe)
    {
      DatabaseMetaData meta=dbcon.getMetaData();
      String txt="";
      txt+=(meta.getDatabaseProductName()+" "+meta.getDatabaseMajorVersion()+"."+meta.getDatabaseMinorVersion()+"\n");
      txt+=(meta.getDriverName()+" "+meta.getDriverVersion()+"\n");
      ArrayList<String> tables = new ArrayList<String>(Arrays.asList("assay","compound","substance","protein_target","source"));
      for (String table: tables)
      {
        ResultSet rset=executeSql(dbcon,"SELECT count(*) FROM "+table);
        if (rset.next()) txt+=("count("+table+"): "+rset.getString(1)+"\n");
        rset.getStatement().close();
      }
      System.err.println(txt);
    }
  }
}