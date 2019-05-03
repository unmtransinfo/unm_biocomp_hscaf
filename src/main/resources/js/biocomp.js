/////////////////////////////////////////////////////////////////////////////
/// biocomp.js -- 
/// 
/// Jeremy Yang
///  19 Aug 2013
/////////////////////////////////////////////////////////////////////////////
function go_zoom_smi2img(smi2img,smicode,opts,w,h)
{
  wwin=w+35;
  hwin=h+35;
  depwin=window.open('','zoomwin','width='+wwin+',height='+hwin+',location=no,resizable=1')
  depwin.resizeTo(wwin,hwin);  // in case window exists
  depwin.focus();
  depwin.location=smi2img+'?'+opts+'&w='+w+'&h='+h+'&smicode='+smicode;
}
function go_zoom_mdl2img(mdl2img,mdlcode,opts,w,h)
{
  wwin=w+35;
  hwin=h+35;
  depwin=window.open('','zoomwin','width='+wwin+',height='+hwin+',location=no,resizable=1')
  depwin.resizeTo(wwin,hwin);  // in case window exists
  depwin.focus();
  depwin.location=mdl2img+'?'+opts+'&w='+w+'&h='+h+'&mdlcode='+mdlcode;
}
function go_zoom_mrv2img(mrv2img,mrvcode,opts,w,h)
{
  wwin=w+35;
  hwin=h+35;
  depwin=window.open('','zoomwin','width='+wwin+',height='+hwin+',location=no,resizable=1')
  depwin.resizeTo(wwin,hwin);  // in case window exists
  depwin.focus();
  depwin.location=mrv2img+'?'+opts+'&w='+w+'&h='+h+'&mrvcode='+mrvcode;
}
function go_zoom_mol2img(mol2img,f,opts,w,h)
{
  wwin=w+35;
  hwin=h+35;
  depwin=window.open('','zoomwin','width='+wwin+',height='+hwin+',location=no,resizable=1')
  depwin.resizeTo(wwin,hwin);  // in case window exists
  depwin.focus();
  depwin.location=mol2img+'?'+opts+'&w='+w+'&h='+h+'&fcode='+f;
}
function go_view3d(view3d,f,size,title,mode)
{
  var wwin=size+25;
  var hwin=size+50;     // extra for buttons
  if (!title) title='view3d';
  var x,y;
  if (!navigator.appName.match('Explorer')) { x=window.screenX+100; y=window.screenY+100; }
  else { x=window.screenLeft+100; y=window.screenTop+100; }
  var viewwin=window.open('',title,'width='+wwin+',height='+hwin+',left='+x+',top='+y+',location=no,resizable=1');
  viewwin.resizeTo(wwin,hwin);  // in case window exists
  viewwin.focus();
  viewwin.location=view3d+'?mode='+mode+'&size='+size+'&title='+title+'&file='+f;
}
function go_view3dAB(view3d,fA,fB,size,title,mode)
{
  var wwin=size+25;
  var hwin=size+50;     // extra for buttons
  if (!title) title='view3d';
  var x,y;
  if (!navigator.appName.match('Explorer')) { x=window.screenX+100; y=window.screenY+100; }
  else { x=window.screenLeft+100; y=window.screenTop+100; }
  var viewwin=window.open('',title,'width='+wwin+',height='+hwin+',left='+x+',top='+y+',location=no,resizable=1');
  viewwin.resizeTo(wwin,hwin);  // in case window exists
  viewwin.focus();
  viewwin.location=view3d+'?mode='+mode+'&size='+size+'&title='+title+'&fileA='+fA+'&fileB='+fB;
}
function go_zoom_chartimg(appurl,opts,w,h,title)
{
  var wwin=w+50;
  var hwin=h+50;
  var win=window.open('',title,'width='+wwin+',height='+hwin+',resizable=1');
  win.resizeTo(wwin,hwin);  // in case window exists
  win.focus();
  win.location=appurl+'?'+'w='+w+'&h='+h+'&'+opts;
}
function go_zoom_img(appurl,opts,w,h,title)
{
  var wwin=w*2;
  var hwin=h*2;
  var win=window.open('',title,'width='+wwin+',height='+hwin+',resizable=1');
  win.resizeTo(wwin,hwin);  // in case window exists
  win.focus();
  win.location=appurl+'?'+'w='+w+'&h='+h+'&'+opts;
}
