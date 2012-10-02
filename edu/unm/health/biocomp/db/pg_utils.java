package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
//import java.util.regex.*;
import java.sql.*;

import org.postgresql.Driver;

//import edu.unm.health.biocomp.http.*;


/**	Static utility methods for PostgreSQL databases&#46;
	<br />
	Uses PostgreSQL JDBC driver (org.postgresql.Driver)&#46;
	<br />
	@author Jeremy J Yang
	@see java.sql.DriverManager
	@see org.postgresql.Driver
*/
public class pg_utils
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
  /**   Executes SQL statement&#46; Use this function for non-queries,
	UPDATE, INSERT, DELETE, CREATE&#46;
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
  /**	Return PostgreSQL connection&#46;
  */
  public static Connection dbConnect(String dbhost,Integer dbport,String dbid,String dbusr,String dbpw)
    throws SQLException
  {
    DriverManager.registerDriver(new org.postgresql.Driver());
    Connection dbcon=DriverManager.getConnection("jdbc:postgresql://"+dbhost+":"+dbport+"/"+dbid,dbusr,dbpw);
    return dbcon;
  }
  /////////////////////////////////////////////////////////////////////////////
}
