;
var F0POCK=[],F12POCK=[];
function sussColor(p){var t=p/100,r=Math.round(150+105*t),g=Math.round(180-120*t),b=Math.round(180-120*t);return"rgb("+r+","+g+","+b+")";}
function fnum(x,d){return(x==null||isNaN(x))?"\u2013":Number(x).toFixed(d==null?2:d);}
function nodeData(id){return NET.nodes.find(function(n){return n.id===id;});}
var nodes=NET.nodes.map(function(d){var an=ANN[d.id];var lab=d.id+(an?("\n"+shortLab(an.label)):"");return{id:d.id,label:lab,value:d.n,color:{background:sussColor(d.suss),border:"#555"},title:d.id+": "+d.n+" members, "+d.suss+"% core SUSS"+(an?(" | "+an.label+" ("+an.pct_domain+"% w/domain, "+an.pct_eff+"% effector)"):"")};});
function shortLab(l){if(!l)return"";var t=l.trim();var map=[[/Auxiliary Activity family 9/i,"AA9/GH61"],[/Glycosyl hydrolases? family 16/i,"GH16"],[/Glycosyl hydrolase family 10/i,"GH10"],[/Glycoside hydrolase 131/i,"GH131"],[/Pectate lyase/i,"Pectate lyase"],[/Cutinase/i,"Cutinase"],[/LysM/i,"LysM"],[/CFEM/i,"CFEM"],[/CVNH/i,"CVNH"],[/GDSL-?like Lipase/i,"GDSL lipase"],[/Hydrophobic surface binding/i,"HsbA"],[/Metallopeptidase/i,"Metallopeptidase"],[/Alternaria alternata allerg/i,"Alt a1 allergen"],[/Domain of unknown function/i,"DUF"],[/Kre9\/KNH/i,"Kre9/KNH"],[/carbonic anh/i,"Carbonic anhydr."],[/Necrosis inducing/i,"NLP/necrosis"],[/Deuterolysin/i,"Deuterolysin"],[/Pregnancy-associated plasma/i,"PAPP/metallopept."],[/Metallo-beta-lactamase/i,"Metallo-\u03b2-lact."],[/novel fold/i,"novel"],[/mixed/i,"mixed"]];for(var k=0;k<map.length;k++){if(map[k][0].test(t))return map[k][1];}t=t.replace(/ domain.*/i,"").replace(/,.*/,"").replace(/family/i,"fam").trim();return t.length>16?t.slice(0,15)+"\u2026":t;}
var edges=NET.edges.map(function(e,i){return{id:"edge-"+i,from:e.from,to:e.to,value:e.tm,title:"mean structural TM "+e.tm.toFixed(2)+" ("+e.n+" cross-family pairs)",color:{color:"#ccc"}};});
var nodeColors={};nodes.forEach(function(n){nodeColors[n.id]=n.color.background;});
var networkNodes=new vis.DataSet(nodes),networkEdges=new vis.DataSet(edges);
var network=new vis.Network(document.getElementById("net"),{nodes:networkNodes,edges:networkEdges},
 {nodes:{shape:"dot",scaling:{min:8,max:40,label:{min:11,max:22}},font:{size:14}},edges:{smooth:false,scaling:{min:1,max:6}},
  physics:{barnesHut:{gravitationalConstant:-3200,springLength:130},stabilization:{iterations:220}},interaction:{hover:true}});
var atlasMode="clusters",singletonSort="acc",singletonSortDir=1,singletonPage=0,singletonPageSize=50,singletonFiltered=SINGLETONS.slice();
var searchMatches=NET.nodes.map(function(n){return n.id;});
function searchNorm(v){return v==null?"":String(v).toLowerCase().trim();}
function esc(v){return String(v==null?"":v).replace(/[&<>"']/g,function(c){return{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c];});}
function memberValues(m){return[m.acc,m.gene,m.eff,m.pfam,m.ipr,m.pdb,m.afdb,m.afdb_hit,m.tm+" tmr",m.tm+" tm",m.novel?"novel":"known"];}
function annotationValues(an){return[an.label,an.top_pfam,an.top_pdb,an.top_ipr].concat((an.members||[]).reduce(function(a,m){return a.concat(memberValues(m));},[]));}
function fieldMatch(id,field,value){
 var an=ANN[id]||{},d=nodeData(id)||{},mem=an.members||[],v=searchNorm(value);
 function any(vals){return vals.some(function(x){return searchNorm(x).indexOf(v)>=0;});}
 if(field==="family"||field==="cluster")return searchNorm(id).indexOf(v)>=0;
 if(field==="gene"||field==="acc"||field==="accession")return mem.some(function(m){return any([m.acc,m.gene]);});
 if(field==="annotation"||field==="anno"||field==="domain"||field==="pfam"||field==="interpro"||field==="pdb"||field==="afdb")return any(annotationValues(an));
 if(field==="effector"||field==="effectorp")return mem.some(function(m){var e=searchNorm(m.eff);if(v==="effector")return e.indexOf("effector")>=0&&e.indexOf("non")<0;if(v==="non-effector"||v==="noneffector")return e.indexOf("non")>=0&&e.indexOf("effector")>=0;return e.indexOf(v)>=0;});
 if(field==="tmr"||field==="deeptmhmm")return mem.some(function(m){return searchNorm(m.tm)===v||searchNorm(m.tm+" tmr").indexOf(v)>=0;});
 if(field==="tm")return fieldMatch(id,"tmr",v)||searchNorm(d.tm)===v;
 if(field==="structtm"||field==="structural-tm")return searchNorm(d.tm).indexOf(v)>=0;
 if(field==="novel")return mem.some(function(m){return searchNorm(m.novel).indexOf(v)>=0||(v==="novel"&&m.novel);});
 if(field==="suss")return searchNorm(d.suss).indexOf(v)>=0;
 var all=[id,d.n,d.tm,d.id_pct,d.suss,d.plddt,d.len,d.maxid,"structural tm "+d.tm,"suss "+d.suss].concat(annotationValues(an));
 if(v==="effector"||v==="non-effector"||v==="noneffector")return fieldMatch(id,"effectorp",v);
 return any(all);
}
function familyMatches(id,query){var terms=searchNorm(query).split(/\s+/).filter(Boolean);return terms.every(function(term){var p=term.indexOf(":");return p>0?fieldMatch(id,term.slice(0,p),term.slice(p+1)):fieldMatch(id,"all",term);});}
function applyNetworkSearch(query){
 var q=searchNorm(query),ids=NET.nodes.map(function(n){return n.id;});searchMatches=q?ids.filter(function(id){return familyMatches(id,q);}):ids;
 var matched={};searchMatches.forEach(function(id){matched[id]=true;});
 networkNodes.update(ids.map(function(id){var hit=!!matched[id]||!q;return{id:id,color:{background:hit?nodeColors[id]:"#e4e9ec",border:hit&&q?"#e67e22":(hit?"#555":"#c7d0d5")},borderWidth:hit&&q?4:1,shadow:hit&&q?{enabled:true,color:"rgba(230,126,34,0.35)",size:12,x:0,y:0}:false,font:{color:hit?"#222":"#a0a9ae"}};}));
 networkEdges.update(edges.map(function(e){var active=!q||(matched[e.from]&&matched[e.to]);return{id:e.id,color:{color:active?"#ccc":"#e7ebed",opacity:active?1:0.25}};}));
 var st=document.getElementById("searchstatus");if(st)st.textContent=q?(searchMatches.length+" cluster"+(searchMatches.length===1?"":"s")):"";
 var clear=document.getElementById("clearsearch");if(clear)clear.style.visibility=q?"visible":"hidden";
 return searchMatches;
}
function singletonMatches(s,query){
 var terms=searchNorm(query).split(/\s+/).filter(Boolean);
 function values(field){
  if(field==="gene"||field==="acc"||field==="accession")return[s.acc,s.gene];
  if(field==="annotation"||field==="anno"||field==="domain"||field==="pfam"||field==="interpro")return[s.label,s.pfam,s.ipr];
  if(field==="pdb")return[s.pdb,s.pdb_tm];
  if(field==="afdb"||field==="foldseek")return[s.afdb,s.afdb_hit,s.afdb_tm,s.pdb,s.pdb_tm];
  if(field==="effector"||field==="effectorp")return[s.eff];
  if(field==="tmr"||field==="deeptmhmm")return[s.tmr,s.tmr+" tmr"];
  if(field==="novel")return[s.novel?"novel":"known"];
  if(field==="pocket")return[s.pocket?"pocket":"no pocket",s.pocket_method,s.pocket_score];
  if(field==="rna"||field==="rnaseq")return[s.rna_condition,s.rna_peak].concat(Object.keys(s.rna||{}));
  return[s.acc,s.gene,s.label,s.eff,s.tmr,s.pfam,s.ipr,s.pdb,s.pdb_tm,s.afdb,s.afdb_hit,s.afdb_tm,s.novel?"novel":"known",s.pocket_method,s.pocket_score,s.plddt,s.length,s.rna_condition,s.rna_peak];
 }
 return terms.every(function(term){var p=term.indexOf(":"),field=p>0?term.slice(0,p):"all",value=p>0?term.slice(p+1):term;return values(field).some(function(x){return searchNorm(x).indexOf(value)>=0;});});
}
function singletonFilterOn(id){var el=document.getElementById(id);return !!(el&&el.checked);}
function singletonSortValue(s,key){
 if(key==="eff")return searchNorm(s.eff);
 if(key==="pocket")return s.pocket_score==null?-Infinity:Number(s.pocket_score);
 if(key==="rna")return s.rna_peak==null?-Infinity:Number(s.rna_peak);
 var value=s[key];return value==null?"":value;
}
function setSingletonSort(key){if(singletonSort===key)singletonSortDir*=-1;else{singletonSort=key;singletonSortDir=1;}singletonPage=0;renderSingletonTable();}
function applySingletonSearch(query){
 singletonFiltered=SINGLETONS.filter(function(s){
  if(query&&!singletonMatches(s,query))return false;
  if(singletonFilterOn("sf-eff")&&!(searchNorm(s.eff).indexOf("effector")>=0&&searchNorm(s.eff).indexOf("non")<0))return false;
  if(singletonFilterOn("sf-novel")&&!s.novel)return false;
  if(singletonFilterOn("sf-pocket")&&!s.pocket)return false;
  if(singletonFilterOn("sf-tmr")&&!(Number(s.tmr)>0))return false;
  return true;
 });
 singletonFiltered.sort(function(a,b){var av=singletonSortValue(a,singletonSort),bv=singletonSortValue(b,singletonSort);if(typeof av==="number"||typeof bv==="number"){av=Number(av);bv=Number(bv);return((isNaN(av)?-Infinity:av)-(isNaN(bv)?-Infinity:bv))*singletonSortDir;}return String(av).localeCompare(String(bv))*singletonSortDir;});
 var maxPage=Math.max(0,Math.ceil(singletonFiltered.length/singletonPageSize)-1);if(singletonPage>maxPage)singletonPage=maxPage;
 renderSingletonRows();
 var st=document.getElementById("searchstatus");if(st)st.textContent=singletonFiltered.length+" singleton"+(singletonFiltered.length===1?"":"s");
 return singletonFiltered;
}
function singletonHeader(label,key){return'<th onclick="setSingletonSort(\''+key+'\')" title="Sort by '+esc(label)+'">'+esc(label)+(singletonSort===key?(singletonSortDir>0?" \u25b2":" \u25bc"):"")+'</th>';}
function singletonHit(name,tm){if(!name)return"\u2013";return'<span title="'+esc(name)+'">'+esc(name)+'</span>'+(tm!=null?'<br><span class="hint">TM '+fnum(tm,3)+'</span>':"");}
function renderSingletonTable(){
 var box=document.getElementById("singletons");if(!box)return;
 box.innerHTML='<div class="singleton-head"><div><h2>Singleton proteins</h2><div class="hint">Independent proteins without a within-dataset structural family. Database Foldseek annotation remains available.</div></div><div class="singleton-actions"><button onclick="dlFilteredSingletons()">Download filtered CSV</button></div></div>'+
 '<div class="singleton-filters"><b>Filter:</b><label><input id="sf-eff" type="checkbox" onchange="singletonFiltersChanged()"> Effector</label><label><input id="sf-novel" type="checkbox" onchange="singletonFiltersChanged()"> Novel</label><label><input id="sf-pocket" type="checkbox" onchange="singletonFiltersChanged()"> Has pocket</label><label><input id="sf-tmr" type="checkbox" onchange="singletonFiltersChanged()"> Has TM helix</label><label>Rows <select id="singletonPageSize" onchange="setSingletonPageSize(this.value)"><option>25</option><option selected>50</option><option>100</option></select></label></div>'+
 '<div class="singleton-table-wrap"><table class="singleton-table"><thead><tr>'+singletonHeader("Protein","acc")+singletonHeader("Annotation","label")+singletonHeader("PDB100 hit","pdb_tm")+singletonHeader("AFDB / Swiss-Prot","afdb_tm")+singletonHeader("Effector / TMR","eff")+singletonHeader("Pocket","pocket")+singletonHeader("pLDDT / length","plddt")+singletonHeader("RNA-seq peak","rna")+'</tr></thead><tbody id="singletonRows"></tbody></table></div><div id="singletonPager" class="pager"></div>';
 applySingletonSearch((document.getElementById("searchinput")||{}).value||"");
}
function renderSingletonRows(){
 var body=document.getElementById("singletonRows");if(!body)return;
 var start=singletonPage*singletonPageSize,rows=singletonFiltered.slice(start,start+singletonPageSize);
 body.innerHTML=rows.map(function(s){
  var eff=esc(s.eff||"\u2013")+(Number(s.tmr)>0?'<br><span class="status-pill warn">'+s.tmr+' TMR</span>':"");
  var pocket=s.pocket?'<span class="status-pill good">'+esc(s.pocket_method||"pocket")+'</span>'+(s.pocket_score!=null?'<br><span class="hint">'+fnum(s.pocket_score,3)+'</span>':""):'<span class="hint">none</span>';
  var tags=(s.novel?'<span class="status-pill novel">novel</span> ':"")+(s.gene?'<span class="hint">'+esc(s.gene)+'</span>':"");
  var rna=s.rna_condition?esc(s.rna_condition)+'<br><span class="hint">'+fnum(s.rna_peak,2)+'</span>':'\u2013';
  return'<tr data-singleton="'+esc(s.id)+'" onclick="showSingleton(\''+String(s.id).replace(/'/g,"\\'")+'\')"><td><b>'+esc(s.acc)+'</b><br>'+tags+'</td><td title="'+esc(s.label)+'">'+esc(s.label)+'</td><td>'+singletonHit(s.pdb,s.pdb_tm)+'</td><td>'+singletonHit(s.afdb||s.afdb_hit,s.afdb_tm)+'</td><td>'+eff+'</td><td>'+pocket+'</td><td>'+fnum(s.plddt,1)+'<br><span class="hint">'+(s.length==null?"\u2013":s.length+' aa')+'</span></td><td>'+rna+'</td></tr>';
 }).join("");
 var pages=Math.max(1,Math.ceil(singletonFiltered.length/singletonPageSize)),pager=document.getElementById("singletonPager");
 if(pager)pager.innerHTML='<button onclick="changeSingletonPage(-1)" '+(singletonPage===0?"disabled":"")+' title="Previous page">\u2039</button><span>'+(singletonFiltered.length?start+1:0)+'\u2013'+Math.min(start+singletonPageSize,singletonFiltered.length)+' of '+singletonFiltered.length+'</span><button onclick="changeSingletonPage(1)" '+(singletonPage>=pages-1?"disabled":"")+' title="Next page">\u203a</button>';
}
function singletonFiltersChanged(){singletonPage=0;applySingletonSearch((document.getElementById("searchinput")||{}).value||"");}
function changeSingletonPage(delta){var pages=Math.max(1,Math.ceil(singletonFiltered.length/singletonPageSize));singletonPage=Math.max(0,Math.min(pages-1,singletonPage+delta));renderSingletonRows();}
function setSingletonPageSize(value){singletonPageSize=Number(value)||50;singletonPage=0;renderSingletonRows();}
function csvCell(value){var text=value==null?"":String(value);return'"'+text.replace(/"/g,'""')+'"';}
function dlFilteredSingletons(){var columns=["acc","gene","label","pdb","pdb_tm","afdb_hit","afdb","afdb_tm","eff","tmr","novel","pfam","ipr","pocket_method","pocket_score","plddt","length","rna_condition","rna_peak"],lines=[columns.join(",")];singletonFiltered.forEach(function(s){lines.push(columns.map(function(c){return csvCell(s[c]);}).join(","));});var blob=new Blob([lines.join("\n")+"\n"],{type:"text/csv"}),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download="singletons_filtered.csv";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function applyActiveSearch(query){return atlasMode==="singletons"?applySingletonSearch(query):applyNetworkSearch(query);}
function clearAtlasSearch(){var input=document.getElementById("searchinput");if(input){input.value="";input.focus();}applyActiveSearch("");}
function setAtlasMode(mode){
 atlasMode=mode;var single=mode==="singletons";document.getElementById("net").classList.toggle("mode-hidden",single);document.getElementById("singletons").classList.toggle("mode-hidden",!single);document.getElementById("modeclusters").className=single?"":"on";document.getElementById("modesingletons").className=single?"on":"";
 var input=document.getElementById("searchinput");if(input)input.placeholder=single?"Search singleton gene, annotation, Foldseek hit":"Search clusters";
 if(single){renderSingletonTable();document.getElementById("side").innerHTML='<p class="hint">Select a singleton row to inspect its structure and evidence.</p>';}else{applyNetworkSearch(input?input.value:"");setTimeout(function(){network.redraw();network.fit();},50);document.getElementById("side").innerHTML='<p class="hint">Click a family node to load data, structure, tree and downloads.</p>';}
}
var searchInput=document.getElementById("searchinput"),searchTimer=null;
if(searchInput){searchInput.addEventListener("input",function(){var q=this.value;clearTimeout(searchTimer);searchTimer=setTimeout(function(){applyActiveSearch(q);},80);});searchInput.addEventListener("keydown",function(e){if(e.key==="Escape"){clearAtlasSearch();e.preventDefault();}else if(e.key==="Enter"&&atlasMode==="clusters"){var found=applyNetworkSearch(this.value);if(found.length===1){network.selectNodes(found);network.focus(found[0],{scale:1.15,animation:true});showFamily(found[0]);}else if(found.length>1){network.fit({nodes:found,animation:true});}e.preventDefault();}});}
var clearSearchButton=document.getElementById("clearsearch");if(clearSearchButton)clearSearchButton.addEventListener("click",clearAtlasSearch);document.getElementById("modesingletons").textContent="Singletons ("+SINGLETONS.length+")";applyNetworkSearch("");
var curFam=null,glviewer=null,structMode="cons",repMode="cartoon",selMembers={},curTree=null,pockMethod="fpocket";

// ---------- Newick parser ----------
function parseNewick(s){
 var i=0; s=s.trim(); if(s[s.length-1]===";")s=s.slice(0,-1);
 function node(){var n={children:[],name:null,length:0};
   if(s[i]==="("){i++;
     do{ if(s[i]===",")i++; n.children.push(node()); }while(s[i]===",");
     i++; // skip )
   }
   // read name
   var nm=""; while(i<s.length&&":,()".indexOf(s[i])<0){nm+=s[i++];}
   if(nm)n.name=nm;
   if(s[i]===":"){i++; var l=""; while(i<s.length&&",()".indexOf(s[i])<0){l+=s[i++];} n.length=parseFloat(l)||0;}
   return n;
 }
 return node();
}
// assign leaf order + depths
function layoutTree(root){
 var leaves=[]; var maxDepth=0;
 (function collect(n,depth){ n.depth=depth;
   if(n.children.length===0){ n.y=leaves.length; leaves.push(n); if(depth>maxDepth)maxDepth=depth; }
   else { n.children.forEach(function(c){collect(c,depth+n.length+0.0);}); }
 })(root,0);
 // cumulative branch length as x
 (function setx(n,x){ n.x=x; n.children.forEach(function(c){setx(c,x+c.length);}); })(root,0);
 // internal y = mean of children y
 (function sety(n){ if(n.children.length){ n.children.forEach(sety); n.y=(n.children[0].y+n.children[n.children.length-1].y)/2; } })(root);
 var maxX=0; (function mx(n){ if(n.x>maxX)maxX=n.x; n.children.forEach(mx); })(root);
 return {leaves:leaves,maxX:maxX};
}
function leavesOf(n){ if(n.children.length===0)return[n.name]; var r=[]; n.children.forEach(function(c){r=r.concat(leavesOf(c));}); return r; }

function renderTree(root, box){
 var lay=layoutTree(root), leaves=lay.leaves, n=leaves.length;
 var rowH=18, padT=10, padB=10, padL=8, labelW=95, plotW=230;
 var H=padT+padB+n*rowH, W=padL+plotW+labelW+10;
 var sx=function(x){return padL+(lay.maxX>0? x/lay.maxX*plotW : 0);};
 var sy=function(y){return padT+rowH/2+y*rowH;};
 var svg='<svg width="'+W+'" height="'+H+'" xmlns="http://www.w3.org/2000/svg">';
 var maxid=PAY[curFam].maxid;
 // draw branches (rectangular)
 (function draw(nd){
   nd.children.forEach(function(c){
     // vertical connector at parent x between children handled once; here draw horizontal + vertical
     draw(c);
   });
   if(nd.children.length){
     // vertical line spanning children y at nd.x
     var y1=sy(nd.children[0].y), y2=sy(nd.children[nd.children.length-1].y);
     svg+='<line x1="'+sx(nd.x)+'" y1="'+y1+'" x2="'+sx(nd.x)+'" y2="'+y2+'" stroke="#333" stroke-width="1.2"/>';
   }
 })(root);
 // horizontal branches + clickable
 var nodeId=0;
 (function draw2(nd,parentX){
   var x0=sx(parentX), x1=sx(nd.x), y=sy(nd.y);
   if(nd!==root){
     svg+='<line class="branch" x1="'+x0+'" y1="'+y+'" x2="'+x1+'" y2="'+y+'" stroke="#333" stroke-width="1.2" '+
          'data-leaves="'+leavesOf(nd).join(",")+'"/>';
   }
   if(nd.children.length===0){
     var mid=(maxid[nd.name]!=null?maxid[nd.name]:1);
     var col=(mid<0.3)?"#c0392b":"#2a6b8a";
     var isHub=(EXTRA[curFam]&&EXTRA[curFam].hub===nd.name);
     svg+='<circle cx="'+x1+'" cy="'+y+'" r="'+(isHub?6:4)+'" fill="#2a6b8a" stroke="'+(isHub?"#e8a90c":"#fff")+'" stroke-width="'+(isHub?3:1)+'" data-leaf="'+nd.name+'" style="cursor:pointer"/>';
     if(isHub){svg+='<text x="'+(x1)+'" y="'+(y-7)+'" font-size="13" fill="#e8a90c" text-anchor="middle">\u2605</text>';}
     svg+='<text x="'+(x1+9)+'" y="'+(y+3.5)+'" font-size="10" fill="'+col+'" font-weight="'+(isHub?"700":"400")+'" data-leaf="'+nd.name+'">'+nd.name+(isHub?'  \u2605 hub':'')+'</text>';
   } else {
     svg+='<circle cx="'+x1+'" cy="'+y+'" r="3.5" fill="#e67e22" stroke="#fff" data-leaves="'+leavesOf(nd).join(",")+'" style="cursor:pointer"><title>click: toggle this clade ('+leavesOf(nd).length+' members)</title></circle>';
   }
   nd.children.forEach(function(c){draw2(c,nd.x);});
 })(root,0);
 svg+='</svg>';
 box.innerHTML=svg;
 curTree={root:root};
 // wire clicks
 box.querySelectorAll("[data-leaf]").forEach(function(el){
   el.addEventListener("click",function(){toggleLeaf(el.getAttribute("data-leaf"));});
 });
 box.querySelectorAll("[data-leaves]").forEach(function(el){
   el.addEventListener("click",function(){toggleClade(el.getAttribute("data-leaves").split(","));});
 });
 paintTree();
}
function paintTree(){
 var box=document.getElementById("treebox"); if(!box)return;
 box.querySelectorAll("circle[data-leaf]").forEach(function(el){
   var m=el.getAttribute("data-leaf"); el.setAttribute("fill",selMembers[m]?"#2a6b8a":"#ccc");
 });
 box.querySelectorAll("text[data-leaf]").forEach(function(el){
   var m=el.getAttribute("data-leaf"); el.setAttribute("opacity",selMembers[m]?"1":"0.35");
 });
}
function toggleLeaf(m){selMembers[m]=!selMembers[m];paintTree();if(structMode==="super")drawStruct();}
function toggleClade(arr){
 // if all selected -> deselect all; else select all
 var allSel=arr.every(function(m){return selMembers[m];});
 arr.forEach(function(m){selMembers[m]=!allSel;});
 paintTree();if(structMode==="super")drawStruct();
}

function row(k,v){return"<tr><td>"+k+"</td><td><b>"+v+"</b></td></tr>";}
function showFamily(id){
 curFam=id;var d=nodeData(id),hasS=!!PAY[id];
 var h='<h2>'+id+' <span class="badge" style="background:'+sussColor(d.suss)+'">'+d.suss+'% core SUSS</span></h2>';
 h+='<table>'+row("Members",d.n)+row("Mean structural sim (TM)",fnum(d.tm))+row("Mean sequence identity",fnum(d.id_pct*100,1)+"%")+row("Max identity",fnum(d.maxid*100,1)+"%")+row("Mean pLDDT",fnum(d.plddt,1))+'</table>';
 var an=ANN[id];
 if(hasS){
  h+='<div class="tabs"><div class="tab on" onclick="tab(0)">Structure + Tree</div><div class="tab" onclick="tab(1)">FoldTree (figure)</div><div class="tab" onclick="tab(2)">Struct sim (TM)</div><div class="tab" onclick="tab(3)">Sequence + MSA</div><div class="tab" onclick="tab(4)">RNAseq</div><div class="tab" onclick="tab(5)">Annotation</div></div>';
  h+='<div id="p0" class="pane on"></div><div id="p1" class="pane"></div><div id="p2" class="pane"></div><div id="p3" class="pane"></div><div id="p4" class="pane"></div><div id="p5" class="pane"></div>';
 } else {
  h+='<div class="tabs"><div class="tab on" onclick="tab(5)">Annotation</div></div><div id="p5" class="pane on"></div>';
  h+='<p class="hint" style="margin-top:6px">Structure/tree/matrix layers computed for demo families F0 &amp; F12; annotation below covers all 39 families.</p>';
 }
 document.getElementById("side").innerHTML=h;
 if(hasS){ buildStructPane(id); setTimeout(function(){initViewer();renderTree(parseNewick(PAY[id].newick),document.getElementById("treebox"));},50); } else { document.getElementById("p5").innerHTML=annHTML(id); }
}
function singletonById(id){return SINGLETONS.find(function(s){return s.id===id;});}
function showSingleton(id){
 curFam=id;var s=singletonById(id),p=PAY[id],ex=EXTRA[id]||{};if(!s||!p)return;
 document.querySelectorAll(".singleton-table tr.selected").forEach(function(row){row.classList.remove("selected");});
 var selected=null;document.querySelectorAll(".singleton-table tr[data-singleton]").forEach(function(row){if(row.getAttribute("data-singleton")===id)selected=row;});if(selected)selected.classList.add("selected");
 var flags=(s.novel?'<span class="status-pill novel">novel</span> ':"")+
   (searchNorm(s.eff).indexOf("effector")>=0&&searchNorm(s.eff).indexOf("non")<0?'<span class="status-pill good">EffectorP</span> ':"")+
   (Number(s.tmr)>0?'<span class="status-pill warn">'+esc(s.tmr)+' TMR</span> ':"");
 var h='<h2>'+esc(s.acc)+'</h2><div style="margin-bottom:6px">'+flags+'</div>'+
   '<table>'+row("Annotation",esc(s.label||"annotation unavailable"))+
   row("Mean pLDDT",fnum(s.plddt,1))+row("Length",s.length==null?"\u2013":esc(s.length)+" aa")+
   row("Pocket",s.pocket?esc(s.pocket_method||"detected")+" \u00b7 score "+fnum(s.pocket_score,3):"not detected")+
   row("RNA-seq peak",s.rna_condition?esc(s.rna_condition)+" \u00b7 "+fnum(s.rna_peak,2):"not available")+'</table>'+
   '<div class="tabs"><div class="tab on" onclick="singletonTab(0)">Structure</div><div class="tab" onclick="singletonTab(1)">RNA-seq</div><div class="tab" onclick="singletonTab(2)">Annotation</div><div class="tab" onclick="singletonTab(3)">Sequence</div></div>'+
   '<div id="sp0" class="pane on"></div><div id="sp1" class="pane"></div><div id="sp2" class="pane"></div><div id="sp3" class="pane"></div>';
 document.getElementById("side").innerHTML=h;
 buildSingletonStructPane(id);
 setTimeout(initViewer,50);
}
function singletonTab(i){
 var tabs=document.querySelectorAll("#side .tab");for(var k=0;k<tabs.length;k++)tabs[k].className="tab"+(k===i?" on":"");
 for(var j=0;j<4;j++){var p=document.getElementById("sp"+j);if(p)p.className="pane"+(j===i?" on":"");}
 var assets=PAY[curFam].assets||{};
 if(i===1&&!document.getElementById("sp1").innerHTML){
  document.getElementById("sp1").innerHTML=assets.rna_svg?'<h3>RNA-seq expression</h3><img src="'+assets.rna_svg+'">'+singletonDlbtn():'<p class="hint">RNA-seq data were not available for this protein.</p>'+singletonDlbtn();
 }
 if(i===2&&!document.getElementById("sp2").innerHTML)document.getElementById("sp2").innerHTML=annHTML(curFam)+singletonDlbtn();
 if(i===3&&!document.getElementById("sp3").innerHTML)buildSequencePane(curFam,false,"sp3");
}
function tab(i){var tabs=document.querySelectorAll(".tab");for(var k=0;k<tabs.length;k++)tabs[k].className="tab"+(tabs[k].getAttribute("onclick").indexOf("tab("+i+")")>=0?" on":"");
  for(var j=0;j<6;j++){var p=document.getElementById("p"+j);if(p)p.className="pane"+(j===i?" on":"");}
  if(i===5){if(!document.getElementById("p5").innerHTML)document.getElementById("p5").innerHTML=annHTML(curFam);return;}
  var a=PAY[curFam]?PAY[curFam].assets:null;if(!a)return;
  if(i===1&&!document.getElementById("p1").innerHTML)document.getElementById("p1").innerHTML='<h3>FoldTree structural phylogeny</h3><img src="'+a.tree_svg+'"><div class="hint">Foldseek \u2192 FoldTree distance \u2192 QuickTree NJ \u2192 '+((EXTRA[curFam]&&EXTRA[curFam].foldtree_rooting_label)||"rooting status unavailable")+'.</div>'+dlbtn();
  if(i===2&&!document.getElementById("p2").innerHTML){var ex2=EXTRA[curFam]||{};var h2='<h3>Structural similarity (Foldseek TM)</h3><img src="'+a.tm_svg+'">'+dlbtn();if(a.tmus_svg){var cr=(ex2.tm_cons_r!=null)?(' \\u00b7 Foldseek\\u2194US-align r='+fnum(ex2.tm_cons_r,3)+', max\\u0394='+fnum(ex2.tm_cons_maxdiff,2)+', '+(ex2.tm_disagree||0)+' pair(s) disagree >0.1'):'';h2+='<h3 style="margin-top:14px">Structural similarity (US-align TM \\u2014 independent algorithm)</h3><div class="hint" style="margin:2px 0 6px">Foldseek builds the families; US-align (TM-align successor) recomputes true TM within the family as an algorithm-independent cross-check.'+cr+'</div><img src="'+a.tmus_svg+'">';}document.getElementById("p2").innerHTML=h2;}
  if(i===3&&!document.getElementById("p3").innerHTML){document.getElementById("p3").innerHTML=(a.id_svg?'<h3>Sequence identity (BLASTp)</h3><img src="'+a.id_svg+'">':'<p class="hint">BLASTp identity matrix was not available for this cluster.</p>')+'<div id="familySequencePane"></div>'+dlbtn();buildSequencePane(curFam,true,"familySequencePane");}
  if(i===4&&!document.getElementById("p4").innerHTML)document.getElementById("p4").innerHTML='<h3>RNA-seq expression (GSE178879)</h3><img src="'+a.rna_svg+'">'+dlbtn();
}
function annHTML(fam){
  var an=ANN[fam];if(!an)return '<p class="hint">No annotation for this family.</p>';
  if(PAY[fam]&&PAY[fam].kind==="singleton"){
    var m=(an.members||[])[0]||{},novel=m.novel===true?"novel":m.novel===false?"not novel":"indeterminate";
    var status=function(v){return v?esc(v):"\u2013";};
    var hit=function(name,tm){return name?esc(name)+(tm!=null?' <span class="hint">(TM '+fnum(tm,3)+')</span>':""):"\u2013";};
    var sh='<h3>Direct protein annotation</h3><table>';
    sh+=row("Label",esc(an.label||"annotation unavailable"));
    sh+=row("Pfam",status(m.pfam));
    sh+=row("InterPro",status(m.ipr));
    sh+=row("Foldseek PDB100",hit(m.pdb,m.pdb_tm));
    sh+=row("Foldseek AFDB / Swiss-Prot",hit(m.afdb||m.afdb_hit,m.afdb_tm));
    sh+=row("EffectorP",status(m.eff));
    sh+=row("DeepTMHMM",Number(m.tm)>0?esc(m.tm)+" predicted transmembrane region(s)":"no transmembrane region");
    sh+=row("Novel status",novel);
    sh+='</table><h3>Analysis status</h3><table>';
    sh+=row("Annotation",status(m.annotation_status));
    sh+=row("PDB100 search",status(m.foldseek_pdb_status));
    sh+=row("AFDB search",status(m.foldseek_afdb_status));
    sh+='</table><div class="hint" style="margin-top:5px">Singleton novelty is reported only when the required domain and structural searches completed. Foldseek database hits remain valid even though no within-dataset family was formed.</div>';
    return sh;
  }
  var h='<h3>Cluster consensus</h3><table>';
  h+='<tr><td>Consensus label</td><td><b>'+an.label+'</b></td></tr>';
  h+='<tr><td>Members with a domain</td><td><b>'+an.pct_domain+'%</b></td></tr>';
  h+='<tr><td>Novel (no domain, no known fold)</td><td><b>'+an.pct_novel+'%</b></td></tr>';
  h+='<tr><td>Predicted effector (EffectorP)</td><td><b>'+an.pct_eff+'%</b></td></tr>';
  if(an.top_pfam!=="\u2014"&&an.top_pfam!=="—")h+='<tr><td>Top Pfam ('+Math.round(an.top_pfam_frac*100)+'% of members)</td><td>'+an.top_pfam+'</td></tr>';
  if(an.top_pdb!=="\u2014"&&an.top_pdb!=="—")h+='<tr><td>Top PDB fold ('+Math.round(an.top_pdb_frac*100)+'%)</td><td>'+an.top_pdb+'</td></tr>';
  if(an.fusion)h+='<tr><td>Multi-domain / fusion</td><td><b style="color:#c0392b">'+an.n_multi+' members</b></td></tr>';
  h+='</table>';
  h+='<h3>Per-member annotation ('+an.n+')</h3>';
  h+='<div style="max-height:340px;overflow:auto"><table style="font-size:11px">';
  h+='<tr style="position:sticky;top:0;background:#eef"><td style="width:auto"><b>Protein</b></td><td><b>EffectorP</b></td><td><b>TM (DeepTMHMM)</b></td><td><b>Pfam domain(s)</b></td><td><b>PDB fold</b></td><td><b>AFDB-SwissProt (Foldseek)</b></td></tr>';
  an.members.forEach(function(m){
    var tag=m.novel?' <span style="color:#c0392b;font-weight:600">novel</span>':'';
    var tm=m.tm>0?' <span style="color:#e67e22">['+m.tm+'TM]</span>':'';
    h+='<tr><td>'+m.acc+tag+'</td><td>'+m.eff.replace(" effector","")+'</td><td>'+(m.tm>0?'<span style="color:#e67e22">'+m.tm+' TMR</span>':'\u2014')+'</td><td>'+m.pfam+'</td><td>'+m.pdb+'</td><td>'+(m.afdb||'\u2014')+'</td></tr>';
  });
  h+='</table></div>';
  h+='<div class="hint" style="margin-top:5px">InterProScan (Pfam+CDD+Gene3D) · Foldseek vs PDB100/AFDB-SwissProt · EffectorP 3.0 · DeepTMHMM. \u201cnovel\u201d = no domain &amp; no structural homolog. [nTM] = predicted transmembrane region.</div>';
  return h;
}
function b64toBlob(b64,mime){var bin=atob(b64),len=bin.length,arr=new Uint8Array(len);for(var i=0;i<len;i++)arr[i]=bin.charCodeAt(i);return new Blob([arr],{type:mime});}
function dlSummary(kind){var b64=(SUMMARY||{})[kind];if(!b64){alert("No "+kind+" summary available in this atlas.");return;}var blob=b64toBlob(b64,"text/csv");var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download="family_summary_"+kind+".csv";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function dlText(txt,fname,mime){var blob=new Blob([txt],{type:mime||"chemical/x-pdb"});var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fname;document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function basePdb(fam){var p=PAY[fam]||{},member=(p.members||[])[0];return REFPDB[fam+"_base"]||REFPDB[fam+"_cons"]||(p.struct||{})[fam]||(p.struct||{})[member]||"";}
function pdbWithValues(pdb,values,def){values=values||{};def=def==null?0:def;return pdb.split("\n").map(function(l){if((l.substring(0,4)==="ATOM"||l.substring(0,6)==="HETATM")&&l.length>=66){var ri=parseInt(l.substring(22,26)),v=Object.prototype.hasOwnProperty.call(values,ri)?Number(values[ri]):def;return l.substring(0,60)+v.toFixed(2).padStart(6," ")+l.substring(66);}return l;}).join("\n");}
function esmPdb(fam){var d=EXTRA[fam]||{};return REFPDB[fam+"_esm"]||(d.esm_values?pdbWithValues(basePdb(fam),d.esm_values,0):"");}
function alignedPdb(fam,m){var pdb=(PAY[fam].struct||{})[m],tr=(PAY[fam].transforms||{})[m];if(!pdb||!tr)return null;var r=tr.rotation,t=tr.translation;return pdb.split("\n").map(function(l){if((l.substring(0,4)==="ATOM"||l.substring(0,6)==="HETATM")&&l.length>=54){var x=parseFloat(l.substring(30,38)),y=parseFloat(l.substring(38,46)),z=parseFloat(l.substring(46,54));if(isFinite(x)&&isFinite(y)&&isFinite(z)){var nx=x*r[0][0]+y*r[1][0]+z*r[2][0]+t[0],ny=x*r[0][1]+y*r[1][1]+z*r[2][1]+t[1],nz=x*r[0][2]+y*r[1][2]+z*r[2][2]+t[2];return l.substring(0,30)+nx.toFixed(3).padStart(8," ")+ny.toFixed(3).padStart(8," ")+nz.toFixed(3).padStart(8," ")+l.substring(54);}}return l;}).join("\n");}
function dlStruct(kind){var fam=curFam,d=EXTRA[fam];
 if(kind==="quality"){var base=basePdb(fam);if(!base){alert("No structure for "+fam);return;}dlText(base,fam+"_AlphaFold_structure.pdb");return;}
 if(kind==="cons"){dlText(REFPDB[fam+"_cons"],fam+"_conservation.pdb");return;}
 if(kind==="esm"){if(!d.has_esm){alert("No ESM scan for "+fam);return;}dlText(esmPdb(fam),fam+"_ESM_tolerance.pdb");return;}
 if(kind==="pocket"){var pock=(pockMethod==="p2rank")?(d.p2rank_resi||[]):(d.fpocket_resi||[]);var ps={};pock.forEach(function(r){ps[r]=999;});dlText(pdbWithValues(basePdb(fam),ps,0),fam+"_"+pockMethod+"_pocket.pdb");return;}
 if(kind==="super"){var mem=PAY[fam].members,sel=[],mdl=1;mem.forEach(function(m){var pdb=alignedPdb(fam,m);if(selMembers[m]&&pdb){sel.push("MODEL "+(mdl++)+"\n"+pdb+"\nENDMDL");}});if(!sel.length){alert("No aligned members selected. Tick members on the tree first.");return;}dlText(sel.join("\n"),fam+"_superposed_"+(mdl-1)+"members.pdb");return;}}
function dlPockResidues(){var fam=curFam,d=EXTRA[fam];
 // map residue number -> 3-letter AA from the reference PDB CA atoms
 var aa={},lines=basePdb(fam).split("\n");
 lines.forEach(function(l){if(l.substring(0,4)==="ATOM"&&l.substring(12,16).trim()==="CA"){aa[parseInt(l.substring(22,26))]=l.substring(17,20).trim();}});
 var rows=["method,residue_number,amino_acid"];
 (d.fpocket_resi||[]).forEach(function(r){rows.push("fpocket,"+r+","+(aa[r]||""));});
 (d.p2rank_resi||[]).forEach(function(r){rows.push("P2Rank,"+r+","+(aa[r]||""));});
 if(rows.length===1){alert("No pocket residues for "+fam);return;}
 var blob=new Blob([rows.join("\n")+"\n"],{type:"text/csv"});var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fam+"_pocket_residues.csv";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function dlXlsx(){var fam=curFam,b64=PAY[fam].assets.xlsx_b64;var blob=b64toBlob(b64,"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fam+"_data.xlsx";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function fastaText(records,members){var out=[],n=0;(members||[]).forEach(function(m){var s=records[m];if(s){out.push(">"+m);for(var i=0;i<s.length;i+=60)out.push(s.substring(i,i+60));n++;}});return{text:out.join("\n")+(out.length?"\n":""),count:n};}
function dlSeqs(){var fam=curFam,p=PAY[fam],fa=fastaText(p.seq||{},p.members);if(!fa.count){alert("No sequences available for "+fam);return;}dlText(fa.text,fam+"_members_"+fa.count+"seqs.fasta","text/plain");}
function dlMsa(){var fam=curFam,p=PAY[fam],fa=fastaText(p.msa||{},p.members);if(fa.count<2){alert("A multiple sequence alignment is not applicable or unavailable for "+fam);return;}dlText(fa.text,fam+"_FoldMason_MSA_"+fa.count+"seqs.fasta","text/plain");}
function dlAllStruct(){var fam=curFam,b64=PAY[fam].assets.structures_zip_b64;if(!b64){alert("No structures embedded for "+fam);return;}var blob=b64toBlob(b64,"application/zip");var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fam+"_member_structures.zip";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function dlMemberStruct(){var fam=curFam,sel=document.getElementById("memSel"),m=sel?sel.value:(PAY[fam].kind==="singleton"?fam:"");var st=(PAY[fam].struct||{})[m];if(!st){alert("No structure for "+(m||fam));return;}dlText(st,m+".pdb");}
function dlMemberSeq(){var fam=curFam,sel=document.getElementById("memSel"),m=sel?sel.value:(PAY[fam].kind==="singleton"?fam:"");var s=(PAY[fam].seq||{})[m];if(!s){alert("No sequence for "+(m||fam));return;}var fa=fastaText((function(){var x={};x[m]=s;return x;})(),[m]);dlText(fa.text,m+".fasta","text/plain");}
function memberDlBar(fam){var st=PAY[fam].struct||{},sq=PAY[fam].seq||{},msa=PAY[fam].msa||{},mem=PAY[fam].members;var ns=0,nq=0,nm=0;mem.forEach(function(m){if(st[m])ns++;if(sq[m])nq++;if(msa[m])nm++;});var opts=mem.map(function(m){return '<option value="'+m+'">'+m+'</option>';}).join("");return '<div class="dlbar" style="margin-top:8px;padding:8px;background:#f7f9fb;border:1px solid #e3e8ee;border-radius:5px">'+'<b>Download members</b> &middot; <span class="hint">'+nq+' sequences, '+ns+' structures</span><br>'+'<button class="dl" onclick="dlSeqs()">\u2b07 All sequences (FASTA)</button> '+(nm>1?'<button class="dl" onclick="dlMsa()">\u2b07 FoldMason MSA (FASTA)</button> ':"")+'<button class="dl" onclick="dlAllStruct()">\u2b07 All structures (ZIP)</button>'+'<br><span class="hint">single member:</span> <select id="memSel">'+opts+'</select> '+'<button class="dl" onclick="dlMemberSeq()">\u2b07 .fasta</button> '+'<button class="dl" onclick="dlMemberStruct()">\u2b07 .pdb</button>'+'</div>';}
function dlbtn(){return '<br><button class="dl" onclick="dlXlsx()">\u2b07 Download all '+curFam+' data (Excel: pockets / FoldTree / Foldseek / US-align / sequence / RNA-seq / per-site)</button>';}
function singletonDlbtn(){return '<div class="dlbar" style="margin-top:8px;padding:8px;background:#f7f9fb;border:1px solid #e3e8ee;border-radius:5px"><b>Download singleton</b><br><button class="dl" onclick="dlMemberSeq()">\u2b07 Sequence (FASTA)</button><button class="dl" onclick="dlMemberStruct()">\u2b07 Structure (PDB)</button><button class="dl" onclick="dlXlsx()">\u2b07 All evidence (Excel)</button></div>';}
var seqViewMode="member";
function sequenceMsaText(msa,members){var present=(members||[]).filter(function(m){return !!msa[m];});if(!present.length)return"No alignment available.";var nlen=Math.min(18,Math.max.apply(null,present.map(function(m){return m.length;}))),alen=Math.max.apply(null,present.map(function(m){return msa[m].length;})),width=60,out=[];for(var pos=0;pos<alen;pos+=width){out.push("Alignment columns "+(pos+1)+"-"+Math.min(pos+width,alen));present.forEach(function(m){var label=(m.length>nlen?m.slice(0,nlen):m).padEnd(nlen," ");out.push(label+"  "+msa[m].substring(pos,pos+width));});out.push("");}return out.join("\n");}
function setSequenceView(mode){seqViewMode=mode;["member","msa"].forEach(function(x){var b=document.getElementById("seqmode_"+x);if(b)b.className=(x===mode?"on":"");});var sel=document.getElementById("seqMember");if(sel)sel.style.display=mode==="member"?"inline-block":"none";renderSequenceViewer();}
function renderSequenceViewer(){var p=PAY[curFam]||{},box=document.getElementById("sequenceViewer");if(!box)return;if(seqViewMode==="msa"){box.textContent=sequenceMsaText(p.msa||{},p.members);return;}var sel=document.getElementById("seqMember"),m=sel?sel.value:(p.members||[])[0],s=(p.seq||{})[m]||"";box.textContent=s?(">"+m+" | "+s.length+" aa\n"+fastaText((function(){var x={};x[m]=s;return x;})(),[m]).text.split("\n").slice(1).join("\n")):"Sequence unavailable.";}
function buildSequencePane(fam,allowMsa,targetId){var p=PAY[fam]||{},mem=p.members||[],msa=p.msa||{},hasMsa=allowMsa&&mem.filter(function(m){return !!msa[m];}).length>1,target=document.getElementById(targetId);if(!target)return;seqViewMode="member";var opts=mem.map(function(m){return'<option value="'+esc(m)+'">'+esc(m)+'</option>';}).join(""),downloadFn=p.kind==="singleton"?"dlMemberSeq()":"dlSeqs()",missingMsa=allowMsa?"FoldMason MSA was not available for this cluster.":"A multiple sequence alignment requires at least two proteins and is not applicable to a singleton.";target.innerHTML='<h3>Sequence viewer</h3><div class="sequence-toolbar"><button id="seqmode_member" class="on" onclick="setSequenceView(\'member\')">Member</button>'+(hasMsa?'<button id="seqmode_msa" onclick="setSequenceView(\'msa\')">MSA</button>':"")+'<select id="seqMember" onchange="renderSequenceViewer()">'+opts+'</select><button onclick="'+downloadFn+'">\u2b07 '+(p.kind==="singleton"?"Sequence":"All sequences")+' (FASTA)</button>'+(hasMsa?'<button onclick="dlMsa()">\u2b07 MSA (FASTA)</button>':"")+'</div>'+(hasMsa?'<div class="hint">MSA is the FoldMason structure-guided alignment used for structural correspondence and superposition.</div>':'<div class="hint">'+missingMsa+'</div>')+'<pre id="sequenceViewer" class="sequence-view"></pre>';renderSequenceViewer();}
function buildStructPane(id){
 var ex=EXTRA[id],mem=PAY[id].members;pockMethod="fpocket";selMembers={};mem.forEach(function(m){selMembers[m]=true;});
 var h='<div><button id="bcons" class="on" onclick="setMode(\'cons\')">Conservation</button>'+
   '<button id="bpocket" onclick="setMode(\'pocket\')">Pocket</button>'+
   '<button id="besm" onclick="setMode(\'esm\')">ESM tolerance</button>'+
   '<button id="bsuper" onclick="setMode(\'super\')">Superpose selected</button></div>'+
   '<div id="pockrow" style="margin-top:5px;display:none"><span class="hint">Pocket method:</span> '+
   '<button id="pk_fpocket" class="on" onclick="setPock(\'fpocket\')">fpocket</button>'+
   '<button id="pk_p2rank" onclick="setPock(\'p2rank\')">P2Rank</button></div>'+
   '<div style="margin-top:5px"><span class="hint">Style:</span> '+
   '<button id="r_cartoon" class="on" onclick="setRep(\'cartoon\')">Cartoon</button>'+
   '<button id="r_surface" onclick="setRep(\'surface\')">Surface</button>'+
   '<button id="r_stick" onclick="setRep(\'stick\')">Stick</button>'+
   '<button id="r_sphere" onclick="setRep(\'sphere\')">Sphere</button>'+
   '<button id="r_line" onclick="setRep(\'line\')">Line</button></div>'+
   '<div id="v3d"></div><div id="leg" class="hint"></div>'+'<div style="margin-top:6px"><span class="hint">Download structure (ChimeraX/PyMOL):</span><br>'+'<button onclick="dlStruct(\'cons\')">Conservation PDB</button>'+'<button onclick="dlStruct(\'esm\')">ESM PDB</button>'+'<button onclick="dlStruct(\'pocket\')" title="pocket residues in B-factor column (999=pocket); see Help for ChimeraX/PyMOL">Pocket-annotated PDB</button>'+'<button onclick="dlStruct(\'super\')">Superpose selected PDB</button></div>';
 h+=memberDlBar(id);
 h+='<h3>FoldTree \u2014 click tips or internal nodes to pick members</h3>'+'<div class="hint" style="margin:2px 0 6px">\u2605 <b>gold star</b> = family hub (highest mean structural similarity to all other members, TM '+fnum(EXTRA[curFam]&&EXTRA[curFam].hub_meanTM,2)+'). '+'The hub is the most representative fold; pocket/conservation are currently projected on '+(EXTRA[curFam]?EXTRA[curFam].ref_used:"?")+' (see note).</div>'+
    '<div style="margin:3px 0"><button onclick="allMem(1)">All</button><button onclick="allMem(0)">None</button>'+
    '<span class="hint" style="margin-left:8px">tip = one member &middot; orange node = whole clade</span></div>'+
    '<div id="treebox"></div>';
 function reslist(a){return (a&&a.length)?a.join(", "):"—";}
 h+='<h3>Selection &amp; sites</h3><table>'+row("fpocket",(ex.fpocket_score!=null?"score "+fnum(ex.fpocket_score,3)+", ":"")+ex.fpocket_resi.length+" res")+row("P2Rank",ex.p2rank_resi.length?("prob "+fnum(ex.p2rank_prob,3)+", "+ex.p2rank_resi.length+" res"):"no pocket")+row("Cys anchors",ex.n_cys)+row("Conserved-buried r",fnum(ex.cons_sasa_r,2))+row("ESM vs conservation r",ex.esm_vs_cons_r!=null?fnum(ex.esm_vs_cons_r,2):"n/a")+row("ESM vs SASA r",ex.esm_vs_sasa_r!=null?fnum(ex.esm_vs_sasa_r,2):"n/a")+(ex.tm_us_mean!=null?row("US-align TM (mean)",fnum(ex.tm_us_mean,3))+row("Foldseek\\u2194US-align r",ex.tm_cons_r!=null?fnum(ex.tm_cons_r,3):"n/a"):"")+'</table>';
 h+='<div class="hint" style="margin-top:4px"><b>fpocket lining residues</b> (ref '+(ex.ref_used||"?")+'): <span style="color:#233">'+reslist(ex.fpocket_resi)+'</span></div>';
 if(ex.p2rank_resi&&ex.p2rank_resi.length)h+='<div class="hint"><b>P2Rank lining residues</b>: <span style="color:#233">'+reslist(ex.p2rank_resi)+'</span></div>';
 h+='<div style="margin-top:5px"><button onclick="dlPockResidues()">\u2b07 Pocket residues (CSV)</button></div>';
 document.getElementById("p0").innerHTML=h;
}
function buildSingletonStructPane(id){
 var ex=EXTRA[id]||{},p=PAY[id],hasStruct=!!((p.struct||{})[id]||REFPDB[id+"_base"]);
 pockMethod=(ex.p2rank_resi||[]).length?"p2rank":"fpocket";
 var h='<div><button id="bquality" class="on" onclick="setMode(\'quality\')">pLDDT</button>'+
   '<button id="bpocket" onclick="setMode(\'pocket\')">Pocket</button>'+
   '<button id="besm" onclick="setMode(\'esm\')">ESM tolerance</button></div>'+
   '<div id="pockrow" style="margin-top:5px;display:none"><span class="hint">Pocket method:</span> '+
   '<button id="pk_fpocket" class="'+(pockMethod==="fpocket"?"on":"")+'" onclick="setPock(\'fpocket\')">fpocket</button>'+
   '<button id="pk_p2rank" class="'+(pockMethod==="p2rank"?"on":"")+'" onclick="setPock(\'p2rank\')">P2Rank</button></div>'+
   '<div style="margin-top:5px"><span class="hint">Style:</span> '+
   '<button id="r_cartoon" class="on" onclick="setRep(\'cartoon\')">Cartoon</button>'+
   '<button id="r_surface" onclick="setRep(\'surface\')">Surface</button>'+
   '<button id="r_stick" onclick="setRep(\'stick\')">Stick</button>'+
   '<button id="r_sphere" onclick="setRep(\'sphere\')">Sphere</button>'+
   '<button id="r_line" onclick="setRep(\'line\')">Line</button></div>'+
   (hasStruct?'<div id="v3d"></div><div id="leg" class="hint"></div>':'<p class="hint">No structure was embedded for this protein.</p>')+
   '<div style="margin-top:6px"><span class="hint">Download structure:</span><br>'+
   '<button onclick="dlStruct(\'quality\')">Original PDB</button>'+
   '<button onclick="dlStruct(\'esm\')">ESM PDB</button>'+
   '<button onclick="dlStruct(\'pocket\')">Pocket-annotated PDB</button></div>';
 h+=singletonDlbtn();
 function reslist(a){return(a&&a.length)?a.join(", "):"\u2013";}
 h+='<h3>Pocket and sequence evidence</h3><table>'+
   row("fpocket",(ex.fpocket_resi||[]).length?("score "+fnum(ex.fpocket_score,3)+" \u00b7 "+ex.fpocket_resi.length+" residues"):"no pocket")+
   row("P2Rank",(ex.p2rank_resi||[]).length?("score "+fnum(ex.p2rank_prob,3)+" \u00b7 "+ex.p2rank_resi.length+" residues"):"no pocket")+
   row("Cysteines",ex.n_cys||0)+row("ESM scan",ex.has_esm?"available":"not available")+'</table>'+
   '<div class="hint"><b>fpocket lining residues:</b> '+reslist(ex.fpocket_resi)+'</div>'+
   '<div class="hint"><b>P2Rank lining residues:</b> '+reslist(ex.p2rank_resi)+'</div>'+
   '<div style="margin-top:5px"><button onclick="dlPockResidues()">\u2b07 Pocket residues (CSV)</button></div>';
 document.getElementById("sp0").innerHTML=h;
}
function allMem(v){for(var m in selMembers)selMembers[m]=!!v;paintTree();if(structMode==="super")drawStruct();}
function initViewer(){var el=document.getElementById("v3d");if(!el)return;el.innerHTML="";glviewer=$3Dmol.createViewer(el,{backgroundColor:"white"});var singleton=PAY[curFam]&&PAY[curFam].kind==="singleton";structMode=singleton?"quality":"cons";repMode="cartoon";setMode(structMode);}
function setPock(x){pockMethod=x;["fpocket","p2rank"].forEach(function(y){var b=document.getElementById("pk_"+y);if(b)b.className=(y===x?"on":"");});if(structMode==="pocket")drawStruct();}
function setMode(m){structMode=m;var pr=document.getElementById("pockrow");if(pr)pr.style.display=(m==="pocket")?"block":"none";["quality","cons","pocket","esm","super"].forEach(function(x){var b=document.getElementById("b"+x);if(b)b.className=(x===m?"on":"");});drawStruct();}
function setRep(r){repMode=r;["cartoon","surface","stick","sphere","line"].forEach(function(x){var b=document.getElementById("r_"+x);if(b)b.className=(x===r?"on":"");});drawStruct();}
function applyStyle(sel,cs){
 if(repMode==="surface"){glviewer.setStyle(sel,{cartoon:{color:(cs.color||"white"),opacity:0.0}});glviewer.addSurface($3Dmol.SurfaceType.VDW,Object.assign({opacity:0.9},cs),sel);}
 else{var st={};
  if(repMode==="cartoon")st.cartoon=cs;
  else if(repMode==="stick")st.stick=Object.assign({radius:0.18},cs);
  else if(repMode==="sphere")st.sphere=Object.assign({scale:0.28},cs);
  else if(repMode==="line")st.line=Object.assign({linewidth:2},cs);
  glviewer.setStyle(sel,st);}
}
function drawStruct(){
 if(!glviewer)return;glviewer.clear();_resLabels={};try{glviewer.removeAllSurfaces();}catch(e){}var d=EXTRA[curFam];
 if(structMode==="super"){
  var cols=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];
  var ci=0,shown=0,mem=PAY[curFam].members;
  mem.forEach(function(m){var pdb=alignedPdb(curFam,m);if(selMembers[m]&&pdb){glviewer.addModel(pdb,"pdb");applyStyle({model:ci},{color:cols[ci%cols.length]});ci++;shown++;}});
  glviewer.zoomTo();glviewer.render();
  document.getElementById("leg").innerHTML="Hub-referenced FoldMason/Kabsch superposition of "+shown+" selected members ("+repMode+"). Tight core = conserved scaffold; splayed loops = variable surface.";
  return;
 }
 var pdbtext=(structMode==="esm"&&d.has_esm)?esmPdb(curFam):(structMode==="quality"?basePdb(curFam):(REFPDB[curFam+"_cons"]||basePdb(curFam)));
 if(!pdbtext){document.getElementById("leg").innerHTML="Structure unavailable.";return;}
 glviewer.addModel(pdbtext,"pdb");
 if(structMode==="esm"){
  if(d.has_esm){applyStyle({},{colorscheme:{prop:"b",gradient:"rwb",min:d.esm_min,max:d.esm_max}});
   document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>constrained <span class="swatch" style="background:#b2182b"></span>tolerant &middot; ESM-1b '+fnum(d.esm_min,1)+"\u2026"+fnum(d.esm_max,1)+" &middot; red = mutation-tolerant (variable), blue = constrained";}
  else{applyStyle({},{color:"lightgrey"});document.getElementById("leg").innerHTML="ESM scan unavailable for this family.";}
 }else if(structMode==="quality"){
  var qmin=d.plddt_min==null?50:d.plddt_min,qmax=d.plddt_max==null?100:d.plddt_max;
  applyStyle({},{colorscheme:{prop:"b",gradient:"rwb",min:qmin,max:qmax}});
  document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>lower confidence <span class="swatch" style="background:#b2182b"></span>higher confidence &middot; AlphaFold pLDDT '+fnum(qmin,1)+"\u2026"+fnum(qmax,1);
 }else if(structMode==="cons"){
  applyStyle({},{colorscheme:{prop:"b",gradient:"rwb",min:d.cons_max,max:d.cons_min}});
  document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>variable <span class="swatch" style="background:#b2182b"></span>conserved &middot; Rate4Site '+fnum(d.cons_min,1)+"\u2026"+fnum(d.cons_max,1);
 } else {
  var pock=(pockMethod==="p2rank")?(d.p2rank_resi||[]):(d.fpocket_resi||[]);
  applyStyle({},{color:"lightgrey"});
  applyStyle({resi:pock},{color:"red"});
  var mlab=(pockMethod==="p2rank")?"P2Rank (ML)":"fpocket (geometric)";
  var sc=(pockMethod==="p2rank")?(d.p2rank_prob!=null?"prob "+fnum(d.p2rank_prob,3):"\u2014"):(d.fpocket_score!=null?"score "+fnum(d.fpocket_score,3):"\u2014");
  var legtxt=pock.length?('<span class="swatch" style="background:lightgrey"></span>scaffold <span class="swatch" style="background:red"></span>'+mlab+' pocket ('+pock.length+' res, '+sc+')'):('<span class="swatch" style="background:lightgrey"></span>'+mlab+' found no pocket for this family');
  document.getElementById("leg").innerHTML=legtxt;
 }
 addResidueClick();
 glviewer.zoomTo();glviewer.render();
}
// click any residue to label it with amino-acid name + position; click again to remove
var _resLabels={};
function addResidueClick(){
 glviewer.setClickable({},true,function(atom){
  var key=atom.chain+":"+atom.resi;
  if(_resLabels[key]){glviewer.removeLabel(_resLabels[key]);delete _resLabels[key];glviewer.render();return;}
  var txt=(atom.resn||"")+" "+atom.resi;
  _resLabels[key]=glviewer.addLabel(txt,{position:{x:atom.x,y:atom.y,z:atom.z},backgroundColor:"black",backgroundOpacity:0.72,fontColor:"white",fontSize:12,borderThickness:0});
  glviewer.render();
 });
 glviewer.setHoverable({},true,function(atom,vw){if(!atom.__hl){atom.__hl=1;}},function(atom){});
}
network.on("click",function(p){if(p.nodes.length)showFamily(p.nodes[0]);});
