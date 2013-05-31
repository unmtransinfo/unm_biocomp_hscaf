package edu.unm.health.biocomp.db;

import java.io.*;
import java.util.*;
import java.sql.*; //Connection

//import org.postgresql.Driver; //Driver
//import com.mysql.jdbc.Driver; //Driver
//import oracle.jdbc.driver.OracleDriver; //OracleDriver

/**	Abstraction class to support multiple db drivers: PostgreSQL, Oracle, MySql.
*/
public class DBCon
{
  private DBCon() {} //disallow default constructor
  private String dbtype;
  private java.sql.Connection con;

  /**	Normal constructor.
	Allowed values for database type: "postgres", "mysql", or "oracle".
  */
  public DBCon(String t,String dbhost,Integer dbport,String dbname,String dbusr,String dbpw)
	throws Exception //DEBUG
  {
    setDBType(t);
    try {
      if (t.equalsIgnoreCase("postgres"))
      {
        //System.err.println("DEBUG: Attempting PostgreSQL connection: jdbc:postgresQL://"+dbhost+":/"+dbname);
        DriverManager.registerDriver(new org.postgresql.Driver());
        java.sql.Connection dbcon=DriverManager.getConnection("jdbc:postgresql://"+dbhost+":"+dbport+"/"+dbname,dbusr,dbpw);
        //if (dbcon!=null) System.err.println("DEBUG: PostgreSQL connection ok: jdbc:postgresQL://"+dbhost+":/"+dbname);
        setConnection(dbcon);
      }
      else if (t.equalsIgnoreCase("mysql"))
      {
        //System.err.println("DEBUG: Attempting MySQL connection: jdbc:mysql://"+dbhost+":/"+dbname);
        DriverManager.registerDriver(new com.mysql.jdbc.Driver());
        java.sql.Connection dbcon=DriverManager.getConnection("jdbc:mysql://"+dbhost+":"+dbport+"/"+dbname,dbusr,dbpw);
        //if (dbcon!=null) System.err.println("DEBUG: MySQL connection ok: jdbc:mysql://"+dbhost+":/"+dbname);
        setConnection(dbcon);
      }
      else if (t.equalsIgnoreCase("oracle"))
      {
        //System.err.println("DEBUG: Attempting Oracle connection: jdbc:oracle:thin:@"+dbhost+":/"+dbname);
        DriverManager.registerDriver(new oracle.jdbc.driver.OracleDriver());
        java.sql.Connection dbcon=DriverManager.getConnection("jdbc:oracle:thin:@"+dbhost+":"+dbport+":"+dbname,dbusr,dbpw);
        //if (dbcon!=null) System.err.println("DEBUG: Oracle connection ok: jdbc:oracle:thin:@"+dbhost+":/"+dbname);
        setConnection(dbcon);
      }
      else
      {
        //System.err.println("DEBUG: Unknown dbtype: \""+t+"\"");
        throw new Exception("Unknown dbtype: \""+t+"\"");
      }
    }
    catch (SQLException e) {
      System.err.println("ERROR: connection failed:"+e.getMessage());
    }

  }

  public String getDBType() { return this.dbtype; }
  public void setDBType(String t) { this.dbtype=t; }
  public java.sql.Connection getConnection() { return this.con; }
  public void setConnection(java.sql.Connection c) { this.con=c; }

  public String serverStatusTxt()
	throws SQLException
  {
    if (this.getDBType().equalsIgnoreCase("postgres")) { return pg_utils.serverStatusTxt(this.con); }
    else if (this.getDBType().equalsIgnoreCase("mysql")) { return mysql_utils.serverStatusTxt(this.con); }
    else if (this.getDBType().equalsIgnoreCase("oracle")) { return ora_utils.serverStatusTxt(this.con); }
    return null;
  }

  public ResultSet executeSql(String sql)
	throws SQLException
  {
    if (this.getDBType().equalsIgnoreCase("postgres")) { return pg_utils.executeSql(this.con,sql); }
    else if (this.getDBType().equalsIgnoreCase("mysql")) { return mysql_utils.executeSql(this.con,sql); }
    else if (this.getDBType().equalsIgnoreCase("oracle")) { return ora_utils.executeSql(this.con,sql); }
    return null;
  }

  public Boolean execute(String sql)
	throws SQLException
  {
    if (this.getDBType().equalsIgnoreCase("postgres")) { return pg_utils.execute(this.con,sql); }
    else if (this.getDBType().equalsIgnoreCase("mysql")) { return mysql_utils.execute(this.con,sql); }
    else if (this.getDBType().equalsIgnoreCase("oracle")) { return ora_utils.execute(this.con,sql); }
    return false;
  }
  public void close()
	throws SQLException
  {
    this.getConnection().close();
  }
}
