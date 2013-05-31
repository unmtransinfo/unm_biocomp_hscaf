package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
//import java.util.regex.*;
import java.sql.*;

import org.postgresql.Driver;

//import edu.unm.health.biocomp.http.*;


/**	Static utility methods for PostgreSQL databases.
	<br />
	Uses PostgreSQL JDBC driver (org.postgresql.Driver).
	<br />
	@author Jeremy J Yang
	@see java.sql.DriverManager
	@see org.postgresql.Driver
*/
public class pg_utils
{
  /////////////////////////////////////////////////////////////////////////////
  /**	Return text with server status information.
  */
  public static String serverStatusTxt(Connection dbcon)
  {
    String txt="";
    try {
      ResultSet rset=executeSql(dbcon,"SHOW SERVER_VERSION");
      if (rset.next())
        txt+=("PostgreSQL server version: "+rset.getString(1));
      rset.getStatement().close();
    }
    catch (Exception e)
    {
      txt+=("error: "+e.getMessage());
    }
    return txt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement. Normally use this function for SELECTs.
  */
  public static ResultSet executeSql(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement(ResultSet.TYPE_FORWARD_ONLY,ResultSet.CONCUR_READ_ONLY);
    ResultSet rset=stmt.executeQuery(sql);
    return rset;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement. Use this function for scrollable ResultSet.
  */
  public static ResultSet executeSqlScrollable(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,ResultSet.CONCUR_READ_ONLY);
    ResultSet rset=stmt.executeQuery(sql);
    return rset;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**   Executes SQL statement. Use this function for non-queries,
	UPDATE, INSERT, DELETE, CREATE, REINDEX, etc.
	Note that autocommit is normally true so there is no
	need to call commit() directly.
  */
  public static boolean execute(Connection dbcon,String sql)
      throws SQLException
  {
    Statement stmt=dbcon.createStatement();
    boolean ok=stmt.execute(sql);
    stmt.close();
    return ok;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Return PostgreSQL connection.
  */
  public static Connection dbConnect(String dbhost,Integer dbport,String dbid,String dbusr,String dbpw)
    throws SQLException
  {
    DriverManager.registerDriver(new org.postgresql.Driver());
    Connection dbcon=DriverManager.getConnection("jdbc:postgresql://"+dbhost+":"+dbport+"/"+dbid,dbusr,dbpw);
    return dbcon;
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Encodes special chars as follows:
	\ 	---to---	'||E'\\'||'
  */
  public static String QuoteString(String str)
  {
    // This is ugly in Java because backslashes are parsed by javac and then by
    // the regex engine, so 4 become 1 and 8 become 2!
    String str2=str.replaceAll("\\\\","'||E'\\\\\\\\'||'");
    //if (!str2.equals(str)) System.err.println("DEBUG: pg_utils.QuoteString: \""+str+"\" -> \""+str2+"\"");
    return str2;
  }
}
