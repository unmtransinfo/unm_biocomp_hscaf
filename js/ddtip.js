/***********************************************
* Cool DHTML tooltip script- © Dynamic Drive DHTML code library (www.dynamicdrive.com)
* This notice MUST stay intact for legal use
* Visit Dynamic Drive at http://www.dynamicdrive.com/ for full source code
***********************************************/

var offsetxpoint=-60; //Customize x offset of tooltip
var offsetypoint=20; //Customize y offset of tooltip
var ie=null;
var ns6=null;
var enabletip=false;
var tipobj=null;

function ddtip_init()
{
  ie=document.all;
  ns6=document.getElementById && !document.all;
  if (ie||ns6)
  {
    tipobj=(document.all?(document.all["ddtooltip"]):(document.getElementById?document.getElementById("ddtooltip"):""));
  }
  document.onmousemove=dd_positiontip;
}

function dd_ietruebody()
{
  return ((document.compatMode&&document.compatMode!="BackCompat")?document.documentElement:document.body);
}

function dd_tip(htm,color,width)
{
  if (ns6||ie)
  {
    if (typeof width!="undefined") tipobj.style.width=width+"px";
    if (typeof color!="undefined" && color!="") tipobj.style.backgroundColor=color;
    tipobj.innerHTML=htm;
    enabletip=true;
    return false;
  }
}

function dd_positiontip(e)
{
  if (enabletip)
  {
    var curX=(ns6)?e.pageX : event.clientX+dd_ietruebody().scrollLeft;
    var curY=(ns6)?e.pageY : event.clientY+dd_ietruebody().scrollTop;
    //Find out how close the mouse is to the corner of the window
    var rightedge=(ie&&!window.opera)? dd_ietruebody().clientWidth-event.clientX-offsetxpoint : window.innerWidth-e.clientX-offsetxpoint-20;
    var bottomedge=(ie&&!window.opera)? dd_ietruebody().clientHeight-event.clientY-offsetypoint : window.innerHeight-e.clientY-offsetypoint-20;
    var leftedge=(offsetxpoint<0)? offsetxpoint*(-1) : -1000;

    //if the horizontal distance isn't enough to accomodate the width of the context menu
    if (rightedge<tipobj.offsetWidth)
    {
      //move the horizontal position of the menu to the left by it's width
      tipobj.style.left=ie? dd_ietruebody().scrollLeft+event.clientX-tipobj.offsetWidth+"px" : window.pageXOffset+e.clientX-tipobj.offsetWidth+"px"
    }
    else if (curX<leftedge)
    {
      tipobj.style.left="5px"
    }
    else
    {
      //position the horizontal position of the menu where the mouse is positioned
      tipobj.style.left=curX+offsetxpoint+"px"
    }

    //same concept with the vertical position
    if (bottomedge<tipobj.offsetHeight)
      tipobj.style.top=ie?dd_ietruebody().scrollTop+event.clientY-tipobj.offsetHeight-offsetypoint+"px":window.pageYOffset+e.clientY-tipobj.offsetHeight-offsetypoint+"px";
    else
      tipobj.style.top=curY+offsetypoint+"px";
    tipobj.style.visibility="visible";
  }
}

function dd_hidetip()
{
  if (ns6||ie)
  {
    enabletip=false;
    tipobj.style.visibility="hidden";
    tipobj.style.left="-1000px";
    tipobj.style.backgroundColor='';
    tipobj.style.width='';
  }
}
