// This does not work; the window zooms but not the structure.
// Need to use api to "Fit to window" -- how?
// Also supported by ChemDraw:
//  'DATAURL="data:chemical/x-cdx;base64,'+cdx_base64
//
function zzoom_chemdraw(smi,tmpcdxurl,edit,h,w)
{
  if (!smi && !tmpcdxurl) return;
  var win=window.open('','zoomwin','width='+w+',height='+h+',resizable=1');
  win.focus();
  var cdhtml='<EMBED TYPE="chemical/x-cdx"';

  // More supported ChemDraw plugin parameters:
  // cdhtml+=' SHOWTOOLSWHENVISIBLE="1"'
  // cdhtml+=' NAME="foobar"'

  if (!edit) cdhtml+=' VIEWONLY="true"';
  cdhtml+=(' WIDTH='+w+' HEIGHT='+h);
  if (tmpcdxurl)
    cdhtml+=(' SRC="'+tmpcdxurl+'">');
  else
    cdhtml+=(' DATAURL="data:chemical/smiles,'+smi+'">');

  cdhtml=('<SCRIPT>cd_insertObjectStr(\''+cdhtml+'\');</SCRIPT>');

  win.document.open('text/html');

  win.document.writeln('<HTML><HEAD>');
  win.document.writeln('<SCRIPT SRC="/js/chemdraw.js"></SCRIPT>');
  win.document.writeln('<SCRIPT>cd_includeWrapperFile(\'/js/\',false);</SCRIPT>');
  win.document.writeln('</HEAD>');
  win.document.writeln('<BODY>');
  win.document.writeln(cdhtml);
  win.document.writeln('</BODY></HTML>');

  win.document.close()
}

function zoom_chemdraw(smi,tmpcdxurl,edit,h,w,cdcgi)
{
  var url=cdcgi+'?';
  if (!smi && !tmpcdxurl) return;
  url+=('h='+h+'&w='+w);
  if (edit) url+=('&edit=true');
  if (tmpcdxurl)
  {
    url+=('&cdxurl='+encodeURIComponent(tmpcdxurl));
  }
  else
  {
    url+=('&smi='+encodeURIComponent(smi));
  }
  var win=window.open(url,'zoomwin','width='+w+',height='+h+',resizable=1');
  win.focus();

}
