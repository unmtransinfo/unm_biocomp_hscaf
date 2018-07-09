var jme_win=false;

// Starts JME using canned jme_window.html.  In future may need to
// use our own html.
function StartJME(loadmol)
{
  jme_win=window.open('/jme/jme_window.html','JME',
        'width=500,height=450,scrollbars=no,resizable=yes');

  // This is partly to allow time for window to initialize:
  if (loadmol && document.mainform.intxt.value)
  {
    if (jme_win.confirm("Load molecule?"))
    {
      jme_win.document.JME.readMolFile(document.mainform.intxt.value);
    }
  }
  jme_win.focus();
}

// This requires that the input molecule format name is "ifmt"
// and the input textarea field name is "intxt".
function fromEditor(smiles,jmefile)
{
  // this function is called from jme_window
  // editor fills variable smiles & jmefile
  if (smiles=="")
  {
    alert("no molecule submitted");
    return;
  }
  var form=document.mainform;
  var molfmt;
  var i;
  for (i=0;i<form.ifmt.length;++i)
  {
    if (form.ifmt.options[i].selected)
      molfmt=form.ifmt.options[i].value;
  }
  if ((molfmt==2 || molfmt==9) && jme_win)      //mdl,sdf
  {
    var molfile=jme_win.document.JME.molFile();
    document.mainform.intxt.value=molfile;
  }
  else
  {
    for (i=0;i<form.ifmt.length;++i)
    {
      if (form.ifmt.options[i].value==1)  //smi
        form.ifmt.options[i].selected=true;
    }
    form.intxt.value=smiles;
  }
}

