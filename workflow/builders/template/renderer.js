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
var atlasMode="clusters",lastFullMode="clusters",lastDomainMode="domains",domainDiagnosticAutoCollapsed=false,singletonSort="acc",singletonSortDir=1,singletonPage=0,singletonPageSize=50,singletonFiltered=SINGLETONS.slice();
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
  if(field==="pocket")return[s.pocket?"pocket":"no pocket",s.pocket_method,s.pocket_metric,s.pocket_value,s.pocket_score];
  if(field==="rna"||field==="rnaseq")return[s.rna_condition,s.rna_peak].concat(Object.keys(s.rna||{}));
  return[s.acc,s.gene,s.label,s.eff,s.tmr,s.pfam,s.ipr,s.pdb,s.pdb_tm,s.afdb,s.afdb_hit,s.afdb_tm,s.novel?"novel":"known",s.pocket_method,s.pocket_metric,s.pocket_value,s.pocket_score,s.plddt,s.length,s.rna_condition,s.rna_peak];
 }
 return terms.every(function(term){var p=term.indexOf(":"),field=p>0?term.slice(0,p):"all",value=p>0?term.slice(p+1):term;return values(field).some(function(x){return searchNorm(x).indexOf(value)>=0;});});
}
function singletonFilterOn(id){var el=document.getElementById(id);return !!(el&&el.checked);}
function singletonSortValue(s,key){
 if(key==="eff")return searchNorm(s.eff);
 if(key==="pocket")return s.pocket_value==null?-Infinity:Number(s.pocket_value);
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
  var pocket=s.pocket?'<span class="status-pill good">'+esc(s.pocket_method||"pocket")+'</span>'+(s.pocket_value!=null?'<br><span class="hint">'+(s.pocket_metric==="probability"?"prob ":"score ")+fnum(s.pocket_value,3)+'</span>':""):'<span class="hint">none</span>';
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
function dlFilteredSingletons(){var columns=["acc","gene","label","pdb","pdb_tm","afdb_hit","afdb","afdb_tm","eff","tmr","novel","pfam","ipr","pocket_method","pocket_metric","pocket_value","pocket_score","plddt","length","rna_condition","rna_peak"],lines=[columns.join(",")];singletonFiltered.forEach(function(s){lines.push(columns.map(function(c){return csvCell(s[c]);}).join(","));});var blob=new Blob([lines.join("\n")+"\n"],{type:"text/csv"}),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download="singletons_filtered.csv";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
var domainFiltered=DOMAIN_FAMILIES.slice(),domainNetwork=null,domainNetworkNodes=null,domainNetworkEdges=null,domainSegmentNetwork=null,domainViewer=null,curDomain=null,curDomainSegment=null,domainSuperpose=false,domainSelected={},domainColorMode="domain",domainRepMode="cartoon",domainBackground="white",domainPocketMethod="fpocket";
function domainMembers(id){return DOMAIN_MEMBERS.filter(function(m){return String(m.domain_family)===String(id);});}
function domainEdges(id){return DOMAIN_EDGES.filter(function(e){return String(e.domain_family)===String(id);});}
function domainNodeColor(lddt){var x=Math.max(0,Math.min(1,Number(lddt)||0));return x>=.7?"#287c8e":(x>=.5?"#58a68d":"#d68b3c");}
function domainAnnotationText(m){var text=(m.overlap_annotations||[]).map(function(x){return x.source+" "+x.label;}).join(" ");if(/necrosis|npp1|nlp/i.test(text))text+=" nlp npp npp1 necrosis-inducing protein";return searchNorm(text);}
function domainEvidenceText(m){return searchNorm([m.acc,m.segment_id,m.parent_family,m.gene,m.eff,m.tmr,m.pfam,m.ipr,m.pdb,m.pdb_tm,m.afdb,m.afdb_hit,m.afdb_tm,m.novel?"novel":"known"].join(" "));}
function domainMemberMatches(m,query){
 var terms=searchNorm(query).split(/\s+/).filter(Boolean),annotation=domainAnnotationText(m),evidence=domainEvidenceText(m),general=searchNorm([m.domain_family,m.start,m.end,annotation,evidence].join(" "));
 return terms.every(function(term){var p=term.indexOf(":"),field=p>0?term.slice(0,p):"",value=p>0?term.slice(p+1):term;if(!value)return true;if(field==="domain"||field==="annotation")return annotation.indexOf(value)>=0;if(field==="gene"||field==="acc"||field==="accession"||field==="protein")return searchNorm([m.acc,m.segment_id,m.gene].join(" ")).indexOf(value)>=0;if(field==="evidence")return evidence.indexOf(value)>=0;if(field==="effectorp")return searchNorm(m.eff).indexOf(value)>=0;if(field==="tmr")return searchNorm(m.tmr).indexOf(value)>=0;if(field==="pdb")return searchNorm([m.pdb,m.pdb_tm].join(" ")).indexOf(value)>=0;if(field==="afdb")return searchNorm([m.afdb,m.afdb_hit,m.afdb_tm].join(" ")).indexOf(value)>=0;return general.indexOf(value)>=0;});
}
function domainMatchingMembers(id,query){var members=domainMembers(id);return query?members.filter(function(m){return domainMemberMatches(m,query);}):members;}
function domainMatches(d,query){return !query||domainMatchingMembers(d.domain_family,query).length>0||searchNorm([d.domain_family,d.top_annotation].join(" ")).indexOf(searchNorm(query))>=0;}
function ensureDomainNetwork(){
 var box=document.getElementById("domains");if(!box||domainNetwork)return;
 var nodes=(DNET.nodes||DOMAIN_FAMILIES).map(function(d){var label=d.domain_family+(d.top_annotation?("\n"+shortLab(d.top_annotation)):"");return{id:d.domain_family,label:label,value:d.n_segments,color:{background:domainNodeColor(d.mean_lddt),border:"#46616c"},title:esc(d.domain_family)+": "+d.n_segments+" segments / "+d.n_proteins+" proteins<br>mean local lDDT "+fnum(d.mean_lddt,3)+(d.top_annotation?"<br>"+esc(d.top_annotation):"")};});
 var edges=(DNET.edges||[]).map(function(e,i){return{id:"dbridge-"+i,from:e.from,to:e.to,value:Math.max(1,e.n),title:e.n+" local structural bridge"+(e.n===1?"":"s")+"<br>mean lDDT "+fnum(e.lddt,3)+"<br>mean probability "+fnum(e.prob,3),color:{color:"#a5b4ba",opacity:.65}};});
 var fixedSmallNetwork=edges.length===0&&nodes.length>1&&nodes.length<=6;if(fixedSmallNetwork)nodes.forEach(function(node,index){node.x=(index-(nodes.length-1)/2)*190;node.y=0;node.fixed={x:true,y:true};});
 domainNetworkNodes=new vis.DataSet(nodes);domainNetworkEdges=new vis.DataSet(edges);
 domainNetwork=new vis.Network(box,{nodes:domainNetworkNodes,edges:domainNetworkEdges},{nodes:{shape:"dot",scaling:{min:9,max:36,label:{min:10,max:18}},font:{size:12}},edges:{smooth:false,scaling:{min:1,max:5}},physics:fixedSmallNetwork?false:{barnesHut:{gravitationalConstant:-2800,springLength:120},stabilization:{iterations:180}},interaction:{hover:true}});
 domainNetwork.on("click",function(p){if(p.nodes.length)showDomain(p.nodes[0]);});
 domainNetwork.once("stabilizationIterationsDone",function(){domainNetwork.stopSimulation();});
}
function applyDomainSearch(query){
 ensureDomainNetwork();domainFiltered=DOMAIN_FAMILIES.filter(function(d){return !query||domainMatches(d,query);});var matched={};domainFiltered.forEach(function(d){matched[d.domain_family]=true;});
 if(domainNetworkNodes)domainNetworkNodes.update(DOMAIN_FAMILIES.map(function(d){var hit=!query||matched[d.domain_family];return{id:d.domain_family,color:{background:hit?domainNodeColor(d.mean_lddt):"#e2e7e9",border:hit&&query?"#c4492d":"#73858d"},borderWidth:hit&&query?4:1,font:{color:hit?"#1d2f37":"#a0a9ae"},shadow:hit&&query?{enabled:true,color:"rgba(196,73,45,.3)",size:10,x:0,y:0}:false};}));
 if(domainNetworkEdges)domainNetworkEdges.update((DNET.edges||[]).map(function(e,i){var active=!query||(matched[e.from]&&matched[e.to]);return{id:"dbridge-"+i,color:{color:active?"#a5b4ba":"#e4e9eb",opacity:active?.65:.18}};}));
 var st=document.getElementById("searchstatus");if(st)st.textContent=domainFiltered.length+" domain famil"+(domainFiltered.length===1?"y":"ies");return domainFiltered;
}
function renderDomainTable(){ensureDomainNetwork();applyDomainSearch((document.getElementById("searchinput")||{}).value||"");}
function domainDiagnosticLabel(status){return status==="borderline"?"Borderline local hit":(status==="filtered_hit"?"Filtered local hit":(status==="unclustered_retained_hit"?"Passed hit, no family":(status==="no_raw_hit"?"No raw local hit":"Diagnostic unavailable")));}
function domainDiagnosticClass(status){return status==="borderline"?"borderline":(status==="filtered_hit"||status==="unclustered_retained_hit"?"filtered":"none");}
function openDomainDiagnosticProtein(acc){var side=document.getElementById("side");if(side&&side.classList.contains("collapsed"))toggleView("side");domainDiagnosticAutoCollapsed=false;openProtein(acc);}
function dlDomainDiagnostics(){var cols=["acc","parent_family","domain_status","best_match","protein_start","protein_end","match_start","match_end","evalue","prob","alntmscore","lddt","fident","protein_coverage","match_coverage","shorter_coverage","alnlen","raw_hit_count","failed_filters","domain_summary","annotation"],lines=[cols.join(",")];DOMAIN_UNASSIGNED.forEach(function(item){var row=Object.assign({annotation:item.label},item);lines.push(cols.map(function(c){return csvCell(row[c]);}).join(","));});dlText(lines.join("\n")+"\n","domain_unassigned_match_diagnostics.csv","text/csv");}
function renderDomainUnassigned(){
 var box=document.getElementById("unassigned"),query=searchNorm((document.getElementById("searchinput")||{}).value||""),rows=DOMAIN_UNASSIGNED.filter(function(item){return!query||searchNorm([item.acc,item.parent_family,item.label,item.domain_status,item.best_match,item.failed_filters,item.domain_summary,item.lddt,item.alntmscore,item.prob,Object.keys(item.expression||{}).join(" ")].join(" ")).indexOf(query)>=0;});
 var counts={borderline:0,filtered_hit:0,unclustered_retained_hit:0,no_raw_hit:0,diagnostic_unavailable:0};DOMAIN_UNASSIGNED.forEach(function(item){counts[item.domain_status]=(counts[item.domain_status]||0)+1;});
 box.innerHTML='<div class="singleton-head"><div><h2>Unassigned domain evidence</h2><div class="hint">These proteins are not in a retained D family. Borderline means the best local hit passed every configured filter except local lDDT and was within the diagnostic margin. This label does not change family membership.</div><div class="domain-diagnostic-summary"><span><b>'+counts.borderline+'</b> borderline</span><span><b>'+(counts.filtered_hit+counts.unclustered_retained_hit)+'</b> filtered / unclustered</span><span><b>'+counts.no_raw_hit+'</b> no raw hit</span></div></div><div class="singleton-actions"><button onclick="dlDomainDiagnostics()">Download diagnostic CSV</button></div></div><div class="singleton-table-wrap"><table class="singleton-table"><thead><tr><th>Protein</th><th>Full-length context</th><th>Domain evidence</th><th>Best local match</th><th>Local scores</th><th>Why no D family</th><th>Annotation</th></tr></thead><tbody>'+rows.map(function(item){var hasHit=!!item.best_match,sourceRange=item.protein_start==null?"":item.protein_start+"\u2013"+item.protein_end,targetRange=item.match_start==null?"":item.match_start+"\u2013"+item.match_end;return'<tr onclick="openDomainDiagnosticProtein(\''+String(item.acc).replace(/'/g,"\\'")+'\')"><td><b>'+esc(item.acc)+'</b><br><span class="hint">'+fnum(item.plddt,1)+' pLDDT \u00b7 '+(item.length==null?"\u2013":item.length+" aa")+'</span></td><td>'+esc(item.parent_family==="singleton"?"protein singleton":item.parent_family)+'</td><td><span class="status-pill '+domainDiagnosticClass(item.domain_status)+'">'+esc(domainDiagnosticLabel(item.domain_status))+'</span><br><span class="hint">'+item.raw_hit_count+' directed hit records</span></td><td>'+(hasHit?'<button onclick="event.stopPropagation();openDomainDiagnosticProtein(\''+String(item.best_match).replace(/'/g,"\\'")+'\')">'+esc(item.best_match)+'</button><br><span class="hint">'+esc(sourceRange)+' matched '+esc(targetRange)+'</span>':"\u2013")+'</td><td class="domain-hit-scores">'+(hasHit?'lDDT '+fnum(item.lddt,3)+'<br>alnTM '+fnum(item.alntmscore,3)+'<br>prob '+fnum(item.prob,3)+'<br>coverage '+fnum(item.protein_coverage,3)+' / '+fnum(item.match_coverage,3):"\u2013")+'</td><td><span class="domain-filter-reason">'+esc(item.failed_filters||item.domain_summary)+'</span></td><td>'+esc(item.label)+'</td></tr>';}).join("")+'</tbody></table></div>';var status=document.getElementById("searchstatus");if(status)status.textContent=rows.length+" unassigned protein"+(rows.length===1?"":"s");
}
function dlDomainFamilies(){var cols=["domain_family","segment_id","acc","start","end","length","parent_family","gene","eff","tmr","pfam","ipr","pdb","pdb_tm","afdb","afdb_tm"],lines=[cols.join(",")];DOMAIN_MEMBERS.forEach(function(m){lines.push(cols.map(function(c){return csvCell(m[c]);}).join(","));});dlText(lines.join("\n")+"\n","domain_families.csv","text/csv");}
function domainWorkspace(m){return(PAY[m.workspace]||PAY[m.parent_family]||PAY[m.acc]||{});}
function domainPdb(m){return(domainWorkspace(m).struct||{})[m.acc]||"";}
function domainSeq(m){return(domainWorkspace(m).seq||{})[m.acc]||"";}
function dlDomainFamily(id){var cols=["segment_id","acc","start","end","length","parent_family","gene","eff","tmr","pfam","ipr","pdb","pdb_tm","afdb","afdb_tm"],lines=[cols.join(",")];domainMembers(id).forEach(function(m){lines.push(cols.map(function(c){return csvCell(m[c]);}).join(","));});dlText(lines.join("\n")+"\n",id+"_segments.csv","text/csv");}
function dlDomainEdges(id){var cols=["source","target","evalue","prob","bits","alntmscore","lddt","fident","alnlen","qcov","tcov","shorter_coverage"],lines=[cols.join(",")];domainEdges(id).forEach(function(e){lines.push(cols.map(function(c){return csvCell(e[c]);}).join(","));});dlText(lines.join("\n")+"\n",id+"_local_foldseek_edges.csv","text/csv");}
function dlDomainSeqs(id){var records={},members=[];domainMembers(id).forEach(function(m){var seq=domainSeq(m);if(seq){records[m.segment_id]=seq.substring(Math.max(0,m.start-1),m.end);members.push(m.segment_id);}});var fa=fastaText(records,members);if(!fa.count){alert("No sequences available for "+id);return;}dlText(fa.text,id+"_domain_segments.fasta","text/plain");}
function dlDomainMsa(id){var wb=(DOMAIN_WORKBENCH.families||{})[id]||{},records=wb.sequence_msa||{},fa=fastaText(records,Object.keys(records));if(fa.count<2){alert("Domain sequence MSA is unavailable for "+id);return;}dlText(fa.text,id+"_MAFFT_domain_sequence_MSA.fasta","text/plain");}
function dlDomainAlignment(id,kind){var wb=domainWorkbench(id),records=kind==="aa"?(wb.structural_msa||{}):(wb.three_di_msa||{}),suffix=kind==="aa"?"FoldMason_AA_MSA":"FoldMason_3Di_MSA",fa=fastaText(records,Object.keys(records));if(fa.count<2){alert(suffix+" is unavailable for "+id);return;}dlText(fa.text,id+"_"+suffix+".fasta","text/plain");}
function dlDomainParentSeqs(id){var records={},members=[];domainMembers(id).forEach(function(m){var seq=domainSeq(m);if(seq&&!records[m.acc]){records[m.acc]=seq;members.push(m.acc);}});var fa=fastaText(records,members);if(!fa.count){alert("No parent sequences available for "+id);return;}dlText(fa.text,id+"_parent_proteins.fasta","text/plain");}
function dlDomainAsset(id,kind,key,filename,mime){var wb=domainWorkbench(id),b64=(wb.assets||{})[key];if(!b64&&BACKEND.enabled){window.location.href=artifactUrl(kind,id);return;}if(!b64){alert("This download is unavailable for "+id);return;}var blob=b64toBlob(b64,mime||"application/zip"),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=filename;document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function dlDomainStructures(id){dlDomainAsset(id,"structures","domain_structures_zip_b64",id+"_domain_structures.zip");}
function dlDomainParents(id){dlDomainAsset(id,"domain_parents","parent_structures_zip_b64",id+"_parent_structures.zip");}
function dlDomainPackage(id){dlDomainAsset(id,"domain_package","package_zip_b64",id+"_domain_package.zip");}
function dlDomainXlsx(id){dlDomainAsset(id,"xlsx","xlsx_b64",id+"_data.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");}
function segmentPdbText(m){var pdb=domainPdb(m);if(!pdb)return"";var lines=pdb.split("\n").filter(function(l){if(l.substring(0,4)!=="ATOM"&&l.substring(0,6)!=="HETATM")return false;var r=parseInt(l.substring(22,26));return isFinite(r)&&r>=m.start&&r<=m.end;});return lines.join("\n")+(lines.length?"\nEND\n":"");}
function dlDomainSegment(){var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(!m)return;var pdb=segmentPdbText(m);if(!pdb&&BACKEND.enabled&&!structureFetchFailed[m.acc]){fetchStructure(m.acc,function(text){if(text)dlDomainSegment();else alert("Structure unavailable for "+m.acc);});return;}if(!pdb){alert("Structure unavailable for "+m.acc);return;}dlText(pdb,m.segment_id.replace(":","_")+".pdb");}
function domainWorkbench(id){return((DOMAIN_WORKBENCH.families||{})[id]||{});}
function showDomain(id){
 var query=arguments.length>1?arguments[1]:"",d=DOMAIN_FAMILIES.find(function(x){return String(x.domain_family)===String(id);}),members=domainMembers(id),wb=domainWorkbench(id);if(!d)return;curDomain=id;domainSuperpose=false;
 var support=d.top_annotation?"Annotation support: "+d.top_annotation_count+"/"+d.n_segments+" segments overlap this label":"Annotation support: no coordinate-overlapping annotation",parents=(d.parent_family_counts||[]).map(function(x){return'<button class="domain-link" onclick="openDomainParent(\''+String(x.family).replace(/'/g,"\\'")+'\')">'+esc(x.family)+' ('+x.n_segments+')</button>';}).join(" ");
 document.getElementById("side").innerHTML='<h2>'+esc(id)+(d.top_annotation?' \u00b7 '+esc(d.top_annotation):"")+'</h2><div class="domain-support"><strong>'+esc(support)+'</strong>. A D family is defined by retained local Foldseek 3Di+AA links; coordinate-overlapping annotation is evidence, not the clustering rule.</div><div class="domain-parent-links"><b>Parent full families:</b> '+(parents||'<span class="hint">not computed</span>')+'</div><div class="domain-summary"><div class="domain-stat"><b>'+d.n_proteins+'</b>proteins</div><div class="domain-stat"><b>'+d.n_segments+'</b>segments</div><div class="domain-stat"><b>'+fnum(d.mean_lddt,3)+'</b>retained-edge lDDT</div><div class="domain-stat"><b>'+fnum(d.mean_alntm,3)+'</b>alignment TM</div><div class="domain-stat"><b>'+fnum(d.mean_aligned_residues,0)+'</b>aligned residues</div><div class="domain-stat"><b>'+d.n_annotated_segments+'</b>annotated segments</div></div><div class="tabs domain-tabs"><div class="tab on" onclick="domainTab(0)">Structure + Network</div><div class="tab" onclick="domainTab(1)">Trees</div><div class="tab" onclick="domainTab(2)">Struct sim</div><div class="tab" onclick="domainTab(3)">Sequence + MSA</div><div class="tab" onclick="domainTab(4)">Conservation + Pockets</div><div class="tab" onclick="domainTab(5)">RNA-seq</div><div class="tab" onclick="domainTab(6)">Annotation</div></div><div id="dp0" class="pane on"></div><div id="dp1" class="pane"></div><div id="dp2" class="pane"></div><div id="dp3" class="pane"></div><div id="dp4" class="pane"></div><div id="dp5" class="pane"></div><div id="dp6" class="pane"></div>';
 buildDomainStructurePane(id,query);
}
function domainTab(i){for(var k=0;k<7;k++){var pane=document.getElementById("dp"+k);if(pane)pane.className="pane"+(k===i?" on":"");}document.querySelectorAll("#side .tab").forEach(function(tab,k){tab.className="tab"+(k===i?" on":"");});if(i===1&&!document.getElementById("dp1").innerHTML)buildDomainTreesPane(curDomain);if(i===2&&!document.getElementById("dp2").innerHTML)buildDomainSimilarityPane(curDomain);if(i===3&&!document.getElementById("dp3").innerHTML)buildDomainSequencePane(curDomain);if(i===4&&!document.getElementById("dp4").innerHTML)buildDomainConservationPane(curDomain);if(i===5&&!document.getElementById("dp5").innerHTML)buildDomainRnaPane(curDomain);if(i===6&&!document.getElementById("dp6").innerHTML)buildDomainAnnotationPane(curDomain);}
function openDomainParent(family){if(/^F\d+$/i.test(family)&&PAY[family]){setAtlasMode("clusters");network.selectNodes([family]);network.focus(family,{scale:1.1,animation:true});showFamily(family);}else if(singletonById(family)){setAtlasMode("singletons");showSingleton(family);}}
function buildDomainStructurePane(id,query){
 var members=domainMembers(id),edges=domainEdges(id),wb=domainWorkbench(id),matched=domainMatchingMembers(id,query||""),initial=((query&&matched[0])||members.find(function(m){return m.segment_id===wb.hub;})||matched[0]||members[0]);
 domainSelected={};if(initial)domainSelected[initial.segment_id]=true;domainColorMode="domain";
 var choices=members.map(function(m){var plen=(domainSeq(m)||"").length||Number(m.parent_length)||Number(m.end)||1,left=Math.max(0,Math.min(100,(Number(m.start)-1)/plen*100)),width=Math.max(1,Math.min(100-left,(Number(m.end)-Number(m.start)+1)/plen*100)),color=domainSegmentColor(m),hasEsm=Object.keys(m.esm_values||{}).length>0;return'<div class="domain-choice" data-domain-choice="'+esc(m.segment_id)+'"><input type="checkbox" id="dc-'+esc(m.segment_id).replace(/[^A-Za-z0-9_-]/g,"-")+'" data-domain-check="'+esc(m.segment_id)+'" onchange="toggleDomainSelection(\''+String(m.segment_id).replace(/'/g,"\\'")+'\',this.checked)"><label onclick="showDomainSegment(\''+String(m.segment_id).replace(/'/g,"\\'")+'\')" title="'+esc(m.segment_id)+'"><span class="swatch" style="background:'+color+'"></span>'+esc(m.acc)+(m.segment_id===wb.hub?" · D hub":"")+(hasEsm?" · ESM":"")+'</label><div class="domain-architecture" title="'+m.start+"–"+m.end+" / "+plen+" aa"+'"><span style="left:'+left+'%;width:'+width+'%;background:'+color+'"></span></div></div>';}).join("");
 var hasSequenceConservation=Object.keys(wb.sequence_conservation||{}).some(function(segment){return Object.keys((wb.sequence_conservation||{})[segment]||{}).length;});
 document.getElementById("dp0").innerHTML=
  '<div class="domain-toolbar"><button id="domainSuperBtn" onclick="setDomainSuperpose(true)">Superpose selected</button><button onclick="setDomainSuperpose(false)">Single member</button><button onclick="fitDomainSegments()">Fit matched domains</button><button onclick="domainViewer&&domainViewer.zoomTo()">Fit full proteins</button></div>'+
  '<div class="domain-toolbar"><span class="hint">Evidence:</span><span class="domain-color-modes"><button id="dcm-domain" class="on" onclick="setDomainColorMode(\'domain\')">Domain</button><button id="dcm-pocket" onclick="setDomainColorMode(\'pocket\')">Pocket</button><button id="dcm-esm" onclick="setDomainColorMode(\'esm\')">ESM</button><button id="dcm-cons" onclick="setDomainColorMode(\'cons\')">Structural conservation</button><button id="dcm-seqcons" onclick="setDomainColorMode(\'seqcons\')" title="'+(hasSequenceConservation?'Rate4Site scores available for at least one D sequence subgroup':'No D sequence subgroup met the Rate4Site requirements')+'">Sequence conservation'+(hasSequenceConservation?'':' · unavailable')+'</button></span></div>'+
  '<div id="domainPocketRow" class="domain-toolbar" style="display:none"><span class="hint">Pocket method:</span><button id="dpk-fpocket" onclick="setDomainPocket(\'fpocket\')">fpocket</button><button id="dpk-p2rank" onclick="setDomainPocket(\'p2rank\')">P2Rank</button></div>'+
  '<div class="domain-toolbar"><span class="hint">Style:</span><button id="dr-cartoon" onclick="setDomainRep(\'cartoon\')">Cartoon</button><button id="dr-surface" onclick="setDomainRep(\'surface\')">Surface</button><button id="dr-stick" onclick="setDomainRep(\'stick\')">Stick</button><button id="dr-sphere" onclick="setDomainRep(\'sphere\')">Sphere</button><button id="dr-line" onclick="setDomainRep(\'line\')">Line</button><span class="hint" style="margin-left:7px">Background:</span><button id="dbg-white" onclick="setDomainBackground(\'white\')">White</button><button id="dbg-black" onclick="setDomainBackground(\'black\')">Black</button></div>'+
  '<div class="domain-selection"><div class="domain-selection-head"><b>Members · one checked = single protein; multiple checked = superpose</b><span><button onclick="selectDomainMembers(\'all\')">All</button><button onclick="selectDomainMembers(\'none\')">None</button><button onclick="selectDomainMembers(\'neighbors\')">Neighbors</button></span></div><div class="domain-selection-list">'+choices+'</div></div><div id="domain3d"></div><div id="domainPosition" class="domain-position"></div><div id="domainLegend" class="domain-superpose-legend"></div><div class="domain-downloads"><b>Download:</b> <button onclick="dlDomainSegment()">Focused domain PDB</button><button onclick="dlDomainFullParent()">Focused full PDB</button><button onclick="dlDomainSelectedSuperposition(false)">Selected domains aligned</button><button onclick="dlDomainSelectedSuperposition(true)">Selected full proteins aligned</button><button onclick="dlDomainStructures(\''+esc(id)+'\')">All domain structures ZIP</button><button onclick="dlDomainParents(\''+esc(id)+'\')">All parent structures ZIP</button><button onclick="dlDomainPackage(\''+esc(id)+'\')">Complete D-family ZIP</button><button onclick="dlDomainXlsx(\''+esc(id)+'\')">All data Excel</button></div><div class="hint">Complete parent proteins are transformed using the matched domain coordinates. Saturated colour marks the aligned domain; pale colour shows the remaining parent context.</div><h3>Segment similarity network</h3><div id="domainSegmentNet"></div><div id="domainEvidence" class="domain-evidence"></div>';
 updateDomainSelectionUi();
 if(domainSegmentNetwork){domainSegmentNetwork.destroy();domainSegmentNetwork=null;}
 var matchedIds={};matched.forEach(function(m){matchedIds[m.segment_id]=true;});
 var snodes=members.map(function(m){var hit=!query||matchedIds[m.segment_id],ann=(m.overlap_annotations||[]).map(function(x){return x.label;}).join(", ");return{id:m.segment_id,label:m.acc+"\n"+m.start+"\u2013"+m.end,value:m.length,color:{background:hit?(m.segment_id===wb.hub?"#d39b2a":m.parent_family==="singleton"?"#5a9f68":"#4388a1"):"#dfe5e8",border:hit&&query?"#c4492d":"#355b69"},borderWidth:hit&&query?4:1,title:esc(m.segment_id)+(ann?"<br>"+esc(ann):"")+"<br>"+esc(m.parent_family==="singleton"?"singleton":m.parent_family)};});
 var sedges=edges.map(function(e,i){return{id:"dseg-"+i,from:e.source,to:e.target,value:Math.max(.1,Number(e.lddt)||.1),title:"local lDDT "+fnum(e.lddt,3)+"<br>alignment TM "+fnum(e.alntmscore,3)+"<br>probability "+fnum(e.prob,3)+"<br>aligned "+e.alnlen+" residues",color:{color:"#9aadb5"}};});
 domainSegmentNetwork=new vis.Network(document.getElementById("domainSegmentNet"),{nodes:new vis.DataSet(snodes),edges:new vis.DataSet(sedges)},{nodes:{shape:"dot",scaling:{min:8,max:24},font:{size:11}},edges:{smooth:false,scaling:{min:1,max:5}},physics:{barnesHut:{gravitationalConstant:-1800,springLength:95},stabilization:{iterations:140}},interaction:{hover:true}});
 domainSegmentNetwork.on("click",function(p){if(p.nodes.length)showDomainSegment(p.nodes[0]);});domainSegmentNetwork.once("stabilizationIterationsDone",function(){domainSegmentNetwork.stopSimulation();});if(initial){domainSegmentNetwork.selectNodes([initial.segment_id]);showDomainSegment(initial.segment_id);}
}
function showDomainSegment(segmentId){
 var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===segmentId;});if(!m)return;curDomainSegment=segmentId;if(!domainSuperpose){domainSelected={};domainSelected[segmentId]=true;}if(domainSegmentNetwork)domainSegmentNetwork.selectNodes([segmentId]);var anns=m.overlap_annotations||[],neighbors=domainEdges(m.domain_family).filter(function(e){return e.source===segmentId||e.target===segmentId;});
 var annHtml=anns.length?anns.map(function(a){return'<span class="domain-tag">'+esc(a.source)+': '+esc(a.label)+' ('+a.start+'\u2013'+a.end+')</span>';}).join(""):'<span class="hint">No Pfam/InterPro call overlaps these coordinates.</span>';
 var edgeRows=neighbors.slice(0,8).map(function(e){var other=e.source===segmentId?e.target:e.source;return'<tr><td>'+esc(other)+'</td><td>lDDT '+fnum(e.lddt,3)+' \u00b7 alnTM '+fnum(e.alntmscore,3)+' \u00b7 '+e.alnlen+' aa</td></tr>';}).join("");
 var plen=(domainSeq(m)||"").length||Number(m.parent_length)||Number(m.end),seglen=Number(m.end)-Number(m.start)+1,where=Number(m.start)<=Math.max(20,plen*.15)?"N-terminal":(Number(m.end)>=plen-Math.max(20,plen*.15)?"C-terminal":"internal"),position=document.getElementById("domainPosition");if(position)position.innerHTML='<b>'+esc(m.acc)+'</b>: domain '+m.start+'\u2013'+m.end+' / '+plen+' aa ('+fnum(seglen/plen*100,1)+'%, '+where+'). The colored interval is the local structure shared by this D family.';
 document.querySelectorAll("[data-domain-choice]").forEach(function(row){row.classList.toggle("focus",row.getAttribute("data-domain-choice")===segmentId);});
 var p2=m.p2rank||{},fp=m.fpocket||{},box=document.getElementById("domainEvidence"),hasEsm=Object.keys(m.esm_values||{}).length>0;if(box)box.innerHTML='<h3>'+esc(m.segment_id)+(m.segment_id===domainWorkbench(m.domain_family).hub?' · family hub':"")+'</h3><div class="domain-actions"><button onclick="openProtein(\''+String(m.acc).replace(/'/g,"\\'")+'\')">Open full protein workbench</button><span class="status-pill">'+esc(m.parent_family==="singleton"?"full singleton":m.parent_family)+'</span><span class="status-pill">ESM '+(hasEsm?"parent-context available":"not scanned for this member")+'</span><span class="status-pill">P2Rank '+(p2.top_score==null?"not available":("score "+fnum(p2.top_score,2)+" · "+(p2.domain_lining_residues||[]).length+" domain residues"))+'</span><span class="status-pill">fpocket '+(fp.top_score==null?"not available":("score "+fnum(fp.top_score,2)+" · "+(fp.domain_lining_residues||[]).length+" domain residues"))+'</span></div><div class="domain-tags">'+annHtml+'</div>'+(edgeRows?'<table>'+edgeRows+'</table>':'<p class="hint">No retained local edge.</p>');
 updateDomainSelectionUi();
 renderDomainStructure(m);
}
function transformDomainPdb(pdb,tr){if(!tr)return pdb;var r=tr.rotation,t=tr.translation;return pdb.split("\n").map(function(l){if((l.substring(0,4)==="ATOM"||l.substring(0,6)==="HETATM")&&l.length>=54){var x=parseFloat(l.substring(30,38)),y=parseFloat(l.substring(38,46)),z=parseFloat(l.substring(46,54));if(isFinite(x)&&isFinite(y)&&isFinite(z)){var nx=x*r[0][0]+y*r[1][0]+z*r[2][0]+t[0],ny=x*r[0][1]+y*r[1][1]+z*r[2][1]+t[1],nz=x*r[0][2]+y*r[1][2]+z*r[2][2]+t[2];return l.substring(0,30)+nx.toFixed(3).padStart(8," ")+ny.toFixed(3).padStart(8," ")+nz.toFixed(3).padStart(8," ")+l.substring(54);}}return l;}).join("\n");}
function domainSegmentColor(m){var wb=domainWorkbench(m.domain_family),palette=["#00879b","#d06b35","#5a8f4e","#8b65a5","#c15d73","#4776a8","#b14f87","#3d9184"],members=domainMembers(m.domain_family),idx=members.findIndex(function(x){return x.segment_id===m.segment_id;});return m.segment_id===wb.hub?"#d39b2a":palette[Math.max(0,idx)%palette.length];}
function domainSelectedMembers(id){return domainMembers(id).filter(function(m){return !!domainSelected[m.segment_id];});}
function updateDomainSelectionUi(){if(!curDomain)return;var selected=domainSelectedMembers(curDomain);document.querySelectorAll("[data-domain-check]").forEach(function(box){box.checked=!!domainSelected[box.getAttribute("data-domain-check")];});var superBtn=document.getElementById("domainSuperBtn");if(superBtn){superBtn.textContent="Superpose selected ("+selected.length+")";superBtn.disabled=selected.length<2;superBtn.className=domainSuperpose?"on":"";}document.querySelectorAll("[id^=dcm-]").forEach(function(button){button.className=button.id==="dcm-"+domainColorMode?"on":"";});document.querySelectorAll("[id^=dr-]").forEach(function(button){button.className=button.id==="dr-"+domainRepMode?"on":"";});document.querySelectorAll("[id^=dbg-]").forEach(function(button){button.className=button.id==="dbg-"+domainBackground?"on":"";});document.querySelectorAll("[id^=dpk-]").forEach(function(button){button.className=button.id==="dpk-"+domainPocketMethod?"on":"";});var pocketRow=document.getElementById("domainPocketRow");if(pocketRow)pocketRow.style.display=domainColorMode==="pocket"?"flex":"none";var shown=domainSuperpose?selected:[DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;})],visible=shown.filter(Boolean),legend=document.getElementById("domainLegend");if(legend)legend.innerHTML=visible.map(function(m){return'<span><i class="swatch" style="background:'+domainSegmentColor(m)+'"></i>'+esc(m.acc+" "+m.start+"–"+m.end)+(m.segment_id===domainWorkbench(m.domain_family).hub?" · hub":"")+'</span>';}).join("")+'<span><i class="swatch" style="background:#b9c1c5"></i>complete parent context</span>'+domainEvidenceLegend(visible);}
function domainEvidenceLegend(shown){if(domainColorMode==="pocket"){var count=shown.reduce(function(total,m){return total+(((domainPocketMethod==="p2rank"?m.p2rank:m.fpocket)||{}).lining_residues||[]).length;},0);return'<span><i class="swatch" style="background:#c83b35"></i>'+esc(domainPocketMethod)+' top-pocket lining residues ('+count+')</span>';}if(domainColorMode==="esm"){var esmAvailable=shown.some(function(m){return Object.keys(m.esm_values||{}).length>0;});return esmAvailable?'<span><i class="swatch" style="background:#2166ac"></i>constrained → <i class="swatch" style="background:#b2182b"></i>tolerant · parent-context ESM-1b</span>':'<span><i class="swatch" style="background:#9aa4a8"></i>parent-context ESM was not run for this selection</span>';}if(domainColorMode==="cons"){var available=shown.some(function(m){return Object.keys(domainStructuralValues(m,domainWorkbench(m.domain_family))).length>0;});return available?'<span><i class="swatch" style="background:#2166ac"></i>variable → <i class="swatch" style="background:#b2182b"></i>structurally conserved</span>':'<span><i class="swatch" style="background:#9aa4a8"></i>structural conservation unavailable for this selection</span>';}if(domainColorMode==="seqcons"){var sequenceAvailable=shown.some(function(m){return Object.keys(((domainWorkbench(m.domain_family).sequence_conservation||{})[m.segment_id]||{})).length>0;});return sequenceAvailable?'<span><i class="swatch" style="background:#2166ac"></i>variable → <i class="swatch" style="background:#b2182b"></i>sequence conserved (D-segment MAFFT + Rate4Site)</span>':'<span><i class="swatch" style="background:#9aa4a8"></i>Rate4Site unavailable for the selected D-segment subgroup</span>';}return"";}
function toggleDomainSelection(segmentId,checked){domainSelected[segmentId]=!!checked;var selected=domainSelectedMembers(curDomain);if(selected.length===1){curDomainSegment=selected[0].segment_id;domainSuperpose=false;}else if(selected.length<2){domainSuperpose=false;}updateDomainSelectionUi();var focus=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;})||selected[0];if(focus)renderDomainStructure(focus);else{var el=document.getElementById("domain3d");if(el)el.innerHTML='<p class="hint" style="padding:12px">Select one or more members.</p>';}}
function selectDomainMembers(mode){var members=domainMembers(curDomain),chosen={};if(mode==="all")members.forEach(function(m){chosen[m.segment_id]=true;});else if(mode==="neighbors"){if(curDomainSegment)chosen[curDomainSegment]=true;domainEdges(curDomain).forEach(function(e){if(e.source===curDomainSegment)chosen[e.target]=true;if(e.target===curDomainSegment)chosen[e.source]=true;});}domainSelected=chosen;if(domainSelectedMembers(curDomain).length<2)domainSuperpose=false;updateDomainSelectionUi();var focus=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(focus)renderDomainStructure(focus);}
function setDomainSuperpose(value){var selected=domainSelectedMembers(curDomain);if(value&&selected.length<2){alert("Select at least two domain members to superpose.");return;}if(!value&&selected.length>1){domainSelected={};if(curDomainSegment)domainSelected[curDomainSegment]=true;}domainSuperpose=value;updateDomainSelectionUi();var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;})||domainSelectedMembers(curDomain)[0];if(m)renderDomainStructure(m);}
function setDomainColorMode(mode){domainColorMode=mode;updateDomainSelectionUi();var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(m)renderDomainStructure(m);}
function setDomainRep(mode){domainRepMode=mode;updateDomainSelectionUi();var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(m)renderDomainStructure(m);}
function setDomainBackground(color){domainBackground=color;updateDomainSelectionUi();if(domainViewer){domainViewer.setBackgroundColor(color);domainViewer.render();}}
function setDomainPocket(method){domainPocketMethod=method;updateDomainSelectionUi();if(domainColorMode==="pocket"){var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(m)renderDomainStructure(m);}}
function dlDomainFullParent(){var m=DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;});if(!m)return;var pdb=domainPdb(m);if(!pdb&&BACKEND.enabled){window.location.href=artifactUrl("structure",m.acc);return;}if(!pdb){alert("Full parent structure unavailable.");return;}dlText(pdb,m.acc+".pdb");}
function dlDomainSelectedSuperposition(fullParents){var selected=domainSelectedMembers(curDomain),wb=domainWorkbench(curDomain);if(selected.length<2){alert("Select at least two domain members.");return;}var missing=selected.filter(function(seg){return !domainPdb(seg)&&BACKEND.enabled&&!structureFetchFailed[seg.acc];});if(missing.length){var pending=missing.length;missing.forEach(function(seg){fetchStructure(seg.acc,function(){pending--;if(!pending)dlDomainSelectedSuperposition(fullParents);});});return;}var models=[];selected.forEach(function(seg,i){var pdb=fullParents?domainPdb(seg):segmentPdbText(seg),tr=(wb.transforms||{})[seg.segment_id];if(pdb&&tr)models.push("MODEL     "+String(i+1).padStart(4," ")+"\nREMARK 900 "+seg.segment_id+" aligned by domain coordinates\n"+transformDomainPdb(pdb,tr).replace(/\nEND\s*\n?$/,"\n")+"ENDMDL");});if(models.length<2){alert("Aligned structures were unavailable for the selected members.");return;}dlText(models.join("\n")+"\nEND\n",curDomain+(fullParents?"_selected_full_parents_aligned.pdb":"_selected_domains_aligned.pdb"));}
function domainResidues(m){var residues=[];for(var residue=Number(m.start);residue<=Number(m.end);residue++)residues.push(residue);return residues;}
function domainStructuralValues(m,wb){var aligned=((wb.structural_msa||{})[m.segment_id]||""),scores=wb.structural_conservation||[],values={},residue=Number(m.start);if(!aligned||!scores.length)return values;for(var column=0;column<aligned.length&&column<scores.length;column++){var symbol=aligned[column];if(symbol==="-"||symbol===".")continue;var score=Number(scores[column]);if(isFinite(score))values[residue]=score;residue++;}return values;}
function fitDomainSegments(){if(!domainViewer)return;var shown=domainSuperpose?domainSelectedMembers(curDomain):[DOMAIN_MEMBERS.find(function(x){return x.segment_id===curDomainSegment;})],sels=[];shown.filter(Boolean).forEach(function(m,index){sels.push({model:index,resi:domainResidues(m)});});domainViewer.zoomTo(sels.length?{or:sels}:{});domainViewer.render();}
function applyDomainStyle(selection,style){if(domainRepMode==="surface"){domainViewer.setStyle(selection,{cartoon:{color:style.color||"white",opacity:0}});domainViewer.addSurface($3Dmol.SurfaceType.VDW,Object.assign({opacity:.9},style),selection);return;}var rep={};if(domainRepMode==="cartoon")rep.cartoon=style;else if(domainRepMode==="stick")rep.stick=Object.assign({radius:.18},style);else if(domainRepMode==="sphere")rep.sphere=Object.assign({scale:.28},style);else rep.line=Object.assign({linewidth:2},style);domainViewer.setStyle(selection,rep);}
function renderDomainStructure(m){
 var el=document.getElementById("domain3d");if(!el)return;
 var selected=domainSuperpose?domainSelectedMembers(m.domain_family):[m],missing=selected.filter(function(seg){return !domainPdb(seg)&&BACKEND.enabled&&!structureFetchFailed[seg.acc];});
 if(missing.length){el.innerHTML='<p class="hint" style="padding:12px">Loading selected full proteins…</p>';var pending=missing.length;missing.forEach(function(seg){fetchStructure(seg.acc,function(){pending--;if(!pending)renderDomainStructure(m);});});return;}
 domainViewer=$3Dmol.createViewer(el,{backgroundColor:domainBackground});
 var wb=domainWorkbench(m.domain_family),modelIndex=0;
 selected.forEach(function(seg){
  var pdb=domainPdb(seg),tr=domainSuperpose?(wb.transforms||{})[seg.segment_id]:{rotation:[[1,0,0],[0,1,0],[0,0,1]],translation:[0,0,0]};if(!pdb||!tr)return;
  var values=domainColorMode==="esm"?(seg.esm_values||{}):(domainColorMode==="cons"?domainStructuralValues(seg,wb):(domainColorMode==="seqcons"?((wb.sequence_conservation||{})[seg.segment_id]||{}):{}));
  var modelText=Object.keys(values).length?pdbWithValues(pdb,values,0):pdb;if(domainSuperpose)modelText=transformDomainPdb(modelText,tr);domainViewer.addModel(modelText,"pdb");
  applyDomainStyle({model:modelIndex},{color:"#b9c1c5",opacity:domainRepMode==="surface"?.52:.9});
  if(domainColorMode==="pocket"){
   var pocketResult=domainPocketMethod==="p2rank"?(seg.p2rank||{}):(seg.fpocket||{}),pocket=pocketResult.lining_residues||[];
   if(pocket.length)applyDomainStyle({model:modelIndex,resi:pocket},{color:"#c83b35",opacity:.98});
  }else{
   var domainStyle={color:domainSegmentColor(seg),opacity:.96};
   if(Object.keys(values).length){
    var vals=Object.keys(values).map(function(k){return Number(values[k]);}).filter(isFinite),mn=Math.min.apply(null,vals),mx=Math.max.apply(null,vals);
    if(domainColorMode==="cons")domainStyle={colorscheme:{prop:"b",gradient:"rwb",min:1,max:0},opacity:.98};
    else domainStyle={colorscheme:{prop:"b",gradient:"rwb",min:mx,max:mn},opacity:.98};
   }else if(domainColorMode!=="domain")domainStyle={color:"#9aa4a8",opacity:.8};
   applyDomainStyle({model:modelIndex,resi:domainResidues(seg)},domainStyle);
  }
  modelIndex++;
 });
 if(!modelIndex){el.innerHTML='<p class="hint" style="padding:12px">Structure unavailable.</p>';return;}
 fitDomainSegments();updateDomainSelectionUi();
}
function buildDomainTreesPane(id){var wb=domainWorkbench(id),foldtrees=wb.foldtree_tree_svgs||{},structural=Object.keys(foldtrees).map(function(metric){return'<section class="tree-block"><h3>FoldTree structural relationship · '+esc(metric)+'</h3><img src="'+foldtrees[metric]+'"><div class="hint">Cropped domain structures; structural relationship evidence, not a sequence phylogeny.</div></section>';}).join(""),guide=wb.foldmason_guide_tree_svg?'<section class="tree-block"><h3>FoldMason structural guide tree</h3><img src="'+wb.foldmason_guide_tree_svg+'"><div class="hint">Guide tree used during multiple-structure alignment; kept separate from FoldTree.</div></section>':"",sequence=(wb.sequence_subgroups||[]).filter(function(group){return group.tree_svg;}).map(function(group){return'<section class="tree-block"><h3>Sequence relationship · '+esc(group.id)+'</h3><img src="'+group.tree_svg+'"><div class="hint">Independent D-segment BLASTp subgroup → D-segment MAFFT → FastTree WAG.</div></section>';}).join("");document.getElementById("dp1").innerHTML=(structural||'<p class="hint">FoldTree was not available or fewer than three domain structures passed.</p>')+guide+(sequence||'<p class="hint">No D-segment sequence-homologous subgroup was large enough for a sequence tree.</p>');}
function buildDomainSimilarityPane(id){var wb=domainWorkbench(id),rows=domainEdges(id).slice().sort(function(a,b){return Number(b.lddt)-Number(a.lddt);}).map(function(e){return'<tr><td>'+esc(e.source)+'</td><td>'+esc(e.target)+'</td><td>'+fnum(e.lddt,3)+'</td><td>'+fnum(e.alntmscore,3)+'</td><td>'+fnum(e.prob,3)+'</td><td>'+e.alnlen+'</td><td>'+fnum(e.qcov,2)+' / '+fnum(e.tcov,2)+'</td></tr>';}).join("");document.getElementById("dp2").innerHTML='<h3>Independent US-align matrix</h3>'+(wb.usalign_matrix_svg?'<img src="'+wb.usalign_matrix_svg+'">':'<p class="hint">US-align matrix unavailable.</p>')+'<h3>Retained local Foldseek links</h3><div class="hint">Foldseek 3Di+AA links define the D family. US-align independently validates cropped segment similarity.</div><div class="domain-edge-table"><table><thead><tr><th>Segment A</th><th>Segment B</th><th>lDDT</th><th>alnTM</th><th>Prob</th><th>Aligned</th><th>Coverage</th></tr></thead><tbody>'+rows+'</tbody></table></div><div class="domain-actions"><button onclick="dlDomainEdges(\''+esc(id)+'\')">Foldseek edges CSV</button><button onclick="dlDomainXlsx(\''+esc(id)+'\')">All data Excel</button></div>';}
function domainMsaBlock(title,records){return'<h3>'+esc(title)+'</h3>'+(Object.keys(records||{}).length?'<pre class="sequence-view">'+esc(sequenceMsaText(records,Object.keys(records)))+'</pre>':'<p class="hint">Not available.</p>');}
function buildDomainSequencePane(id){var wb=domainWorkbench(id),identity=wb.sequence_identity_matrix_svg?'<h3>Domain-segment BLASTp identity</h3><div class="hint">The heatmap is calculated from an independent all-vs-all search of the cropped D segments. Similarity elsewhere in the complete parent proteins cannot create a D sequence relationship.</div><img src="'+wb.sequence_identity_matrix_svg+'">':'<p class="hint">Domain BLASTp identity matrix was not available.</p>',structuralIdentity=wb.structural_alignment_identity_matrix_svg?'<h3>FoldMason-aligned amino-acid identity</h3><div class="hint">This second matrix measures amino-acid identity over structurally corresponding FoldMason columns and is kept separate from the BLASTp sequence search.</div><img src="'+wb.structural_alignment_identity_matrix_svg+'">':"";document.getElementById("dp3").innerHTML=identity+structuralIdentity+'<div class="hint">Local Foldseek structure links define the D family. Domain-segment BLASTp defines sequence-homologous subgroups; MAFFT and sequence trees are built independently inside those subgroups. FoldMason covers the complete structural family.</div>'+domainMsaBlock("MAFFT D-segment sequence MSA",wb.sequence_msa||{})+domainMsaBlock("FoldMason structural MSA · amino acids",wb.structural_msa||{})+domainMsaBlock("FoldMason structural MSA · 3Di",wb.three_di_msa||{})+'<div class="domain-actions"><button onclick="dlDomainSeqs(\''+esc(id)+'\')">All domain sequences</button><button onclick="dlDomainParentSeqs(\''+esc(id)+'\')">All parent sequences</button><button onclick="dlDomainMsa(\''+esc(id)+'\')">MAFFT MSA</button><button onclick="dlDomainAlignment(\''+esc(id)+'\',\'aa\')">FoldMason AA MSA</button><button onclick="dlDomainAlignment(\''+esc(id)+'\',\'3di\')">FoldMason 3Di MSA</button></div>';}
function buildDomainConservationPane(id){var wb=domainWorkbench(id),members=domainMembers(id),sequenceConservation=wb.sequence_conservation||{},rows=members.map(function(m){var p2=m.p2rank||{},fp=m.fpocket||{},esm=Object.keys(m.esm_values||{}).length,seq=Object.keys(sequenceConservation[m.segment_id]||{}).length;return'<tr><td><b>'+esc(m.segment_id)+'</b></td><td>'+((p2.domain_lining_residues||[]).join(", ")||"none")+'</td><td>'+fnum(p2.top_score,3)+' / '+fnum(p2.top_probability,3)+'</td><td>'+((fp.domain_lining_residues||[]).join(", ")||"none")+'</td><td>'+fnum(fp.top_score,3)+'</td><td>'+esm+'</td><td>'+seq+'</td></tr>';}).join(""),structuralSegments=members.filter(function(m){return Object.keys(domainStructuralValues(m,wb)).length;}).length,sequenceSegments=members.filter(function(m){return Object.keys(sequenceConservation[m.segment_id]||{}).length;}).length,subgroups=(wb.sequence_subgroups||[]).map(function(group){var required=Number(group.minimum_conservation_sequences)||4,status=group.sequence_conservation_status==="complete"?"complete":(group.n_sequences+" / "+required+" sequences · "+String(group.sequence_conservation_status||"not applicable").replace(/_/g," "));return'<tr><td><b>'+esc(group.id)+'</b></td><td>'+Number(group.n_sequences||0)+'</td><td>'+Number(group.n_relationships||0)+'</td><td>'+esc(status)+'</td></tr>';}).join("");document.getElementById("dp4").innerHTML='<h3>Structural and sequence conservation</h3><div class="hint">FoldMason structural conservation is calculated independently across the cropped structures in this D family. Rate4Site uses only independently searched D-segment BLASTp subgroups followed by D-segment MAFFT. P2Rank, fpocket, and parent-context ESM are single-protein evidence mapped onto the domain coordinates.</div><div class="domain-summary"><div class="domain-stat"><b>'+structuralSegments+'</b>segments with FoldMason lDDT</div><div class="domain-stat"><b>'+sequenceSegments+'</b>segments with Rate4Site scores</div><div class="domain-stat"><b>'+members.filter(function(m){return(m.pocket_residues||[]).length;}).length+'</b>segments overlapping a pocket</div><div class="domain-stat"><b>'+members.filter(function(m){return Object.keys(m.esm_values||{}).length;}).length+'</b>segments with parent-context ESM</div></div><h3>D-segment sequence subgroup status</h3><div class="domain-edge-table"><table><thead><tr><th>Subgroup</th><th>D sequences</th><th>BLASTp relationships</th><th>Rate4Site</th></tr></thead><tbody>'+subgroups+'</tbody></table></div><h3>Per-segment evidence</h3><div class="domain-edge-table"><table><thead><tr><th>Segment</th><th>P2Rank domain residues</th><th>score / probability</th><th>fpocket domain residues</th><th>score</th><th>ESM sites</th><th>Rate4Site sites</th></tr></thead><tbody>'+rows+'</tbody></table></div>';}
function buildDomainRnaPane(id){var wb=domainWorkbench(id),members=domainMembers(id),conditions={},unique={};members.forEach(function(m){unique[m.acc]=m;Object.keys(m.expression||{}).forEach(function(k){conditions[k]=true;});});var keys=Object.keys(conditions),head='<tr><th>Parent protein</th>'+keys.map(function(k){return'<th>'+esc(k)+'</th>';}).join("")+'</tr>',rows=Object.keys(unique).sort().map(function(acc){var m=unique[acc];return'<tr><td>'+esc(acc)+'</td>'+keys.map(function(k){return'<td>'+fnum((m.expression||{})[k],2)+'</td>';}).join("")+'</tr>';}).join("");document.getElementById("dp5").innerHTML=keys.length?'<h3>RNA-seq by unique parent protein</h3><div class="hint">Expression is a parent-protein measurement. A protein appearing in multiple domain segments is counted once.</div>'+(wb.rna_svg?'<img src="'+wb.rna_svg+'">':'')+'<div class="domain-edge-table"><table><thead>'+head+'</thead><tbody>'+rows+'</tbody></table></div>':'<p class="hint">RNA-seq data were not available for this run.</p>';}
function buildDomainAnnotationPane(id){var rows=domainMembers(id).map(function(m){var overlap=(m.overlap_annotations||[]).map(function(a){return a.source+": "+a.label+" ("+a.start+"–"+a.end+")";}).join("; ")||"none",protein=[m.gene,m.eff,m.pfam,m.ipr,m.pdb,m.afdb||m.afdb_hit].filter(Boolean).join("; ")||"none";return'<tr><td><b>'+esc(m.segment_id)+'</b><br><button onclick="openProtein(\''+String(m.acc).replace(/'/g,"\\'")+'\')">Open full protein</button></td><td>'+esc(m.parent_family==="singleton"?"protein singleton":m.parent_family)+'</td><td>'+esc(overlap)+'</td><td>'+esc(protein)+'</td></tr>';}).join("");document.getElementById("dp6").innerHTML='<h3>Coordinate-aware annotation</h3><div class="hint">Domain calls overlapping the segment are separated from parent-protein annotation and external Foldseek hits.</div><div class="domain-edge-table"><table><thead><tr><th>Segment</th><th>Parent F family</th><th>Overlapping domain call</th><th>Parent-protein evidence</th></tr></thead><tbody>'+rows+'</tbody></table></div><div class="domain-actions"><button onclick="dlDomainFamily(\''+esc(id)+'\')">Segments CSV</button><button onclick="dlDomainPackage(\''+esc(id)+'\')">Complete D-family ZIP</button></div>';}
function fullFamilyForAcc(acc){var ids=Object.keys(PAY);for(var i=0;i<ids.length;i++){var p=PAY[ids[i]]||{};if(p.kind!=="singleton"&&(p.members||[]).indexOf(acc)>=0)return ids[i];}return null;}
function openProtein(acc){var singleton=singletonById(acc);if(singleton){setAtlasMode("singletons");showSingleton(acc);return;}var family=fullFamilyForAcc(acc);if(family){setAtlasMode("clusters");network.selectNodes([family]);network.focus(family,{scale:1.1,animation:true});showFamily(family);}}
function applyActiveSearch(query){return atlasMode==="singletons"?applySingletonSearch(query):(atlasMode==="domains"?applyDomainSearch(query):(atlasMode==="unassigned"?(renderDomainUnassigned(),DOMAIN_UNASSIGNED):applyNetworkSearch(query)));}
function clearAtlasSearch(){var input=document.getElementById("searchinput");if(input){input.value="";input.focus();}applyActiveSearch("");}
function setAnalysisAxis(axis){if(axis==="full"){if(ANALYSIS_SCOPE==="domain")return;setAtlasMode(lastFullMode);}else{if(ANALYSIS_SCOPE==="full")return;setAtlasMode(lastDomainMode);}}
function setAtlasMode(mode){
 if((mode==="domains"||mode==="unassigned")&&ANALYSIS_SCOPE==="full")return;if((mode==="clusters"||mode==="singletons")&&ANALYSIS_SCOPE==="domain")return;
 var detailSide=document.getElementById("side");if(mode==="unassigned"){if(detailSide&&!detailSide.classList.contains("collapsed")){toggleView("side");domainDiagnosticAutoCollapsed=true;}}else if(domainDiagnosticAutoCollapsed){if(detailSide&&detailSide.classList.contains("collapsed"))toggleView("side");domainDiagnosticAutoCollapsed=false;}
 atlasMode=mode;var single=mode==="singletons",domains=mode==="domains",unassigned=mode==="unassigned",domainAxis=domains||unassigned;if(domainAxis)lastDomainMode=mode;else lastFullMode=mode;document.getElementById("net").classList.toggle("mode-hidden",single||domains||unassigned);document.getElementById("singletons").classList.toggle("mode-hidden",!single);document.getElementById("domains").classList.toggle("mode-hidden",!domains);document.getElementById("unassigned").classList.toggle("mode-hidden",!unassigned);document.getElementById("modeclusters").className=(!single&&!domains&&!unassigned)?"on":"";document.getElementById("modedomains").className=domains?"on":"";document.getElementById("modeunassigned").className=unassigned?"on":"";document.getElementById("modesingletons").className=single?"on":"";document.getElementById("scopefull").className=domainAxis?"":"on";document.getElementById("scopedomain").className=domainAxis?"on":"";document.getElementById("fullModeTabs").classList.toggle("mode-hidden",domainAxis);document.getElementById("domainModeTabs").classList.toggle("mode-hidden",!domainAxis);
 var input=document.getElementById("searchinput");if(input)input.placeholder=single?"Search protein singleton":(domains?"Search D family, protein, or segment":(unassigned?"Search unassigned protein, best hit, or filter":"Search full-length families"));
 if(single){renderSingletonTable();document.getElementById("side").innerHTML='<p class="hint">Select a singleton row to inspect its structure and evidence.</p>';}else if(domains){renderDomainTable();document.getElementById("side").innerHTML='<p class="hint">The overview network connects D families when retained local Foldseek hits bridge their segment communities. Select a node for the segment-level workbench.</p>';setTimeout(function(){if(domainNetwork){domainNetwork.redraw();domainNetwork.fit();}},50);}else if(unassigned){renderDomainUnassigned();document.getElementById("side").innerHTML='<p class="hint">Select a protein to open its full-protein evidence workspace.</p>';}else{applyNetworkSearch(input?input.value:"");setTimeout(function(){network.redraw();network.fit();},50);document.getElementById("side").innerHTML='<p class="hint">Click a family node to load data, structure, trees and downloads.</p>';}
}
var searchInput=document.getElementById("searchinput"),searchTimer=null;
if(searchInput){searchInput.addEventListener("input",function(){var q=this.value;clearTimeout(searchTimer);searchTimer=setTimeout(function(){applyActiveSearch(q);},80);});searchInput.addEventListener("keydown",function(e){if(e.key==="Escape"){clearAtlasSearch();e.preventDefault();}else if(e.key==="Enter"&&atlasMode==="clusters"){var found=applyNetworkSearch(this.value);if(found.length===1){network.selectNodes(found);network.focus(found[0],{scale:1.15,animation:true});showFamily(found[0]);}else if(found.length>1){network.fit({nodes:found,animation:true});}e.preventDefault();}else if(e.key==="Enter"&&atlasMode==="domains"){var query=this.value,df=applyDomainSearch(query);if(df.length===1){var did=df[0].domain_family;domainNetwork.selectNodes([did]);domainNetwork.focus(did,{scale:1.15,animation:true});showDomain(did,query);}else if(df.length>1){domainNetwork.fit({nodes:df.map(function(x){return x.domain_family;}),animation:true});}e.preventDefault();}});}
var clearSearchButton=document.getElementById("clearsearch");if(clearSearchButton)clearSearchButton.addEventListener("click",clearAtlasSearch);document.getElementById("modeclusters").textContent="Families ("+NET.nodes.length+")";document.getElementById("modesingletons").textContent="Singletons ("+SINGLETONS.length+")";document.getElementById("modedomains").textContent="Families ("+DOMAIN_FAMILIES.length+")";document.getElementById("modeunassigned").textContent="Unassigned / filtered ("+DOMAIN_UNASSIGNED.length+")";if(ANALYSIS_SCOPE==="full"){document.getElementById("scopedomain").style.display="none";document.getElementById("domainModeTabs").style.display="none";}if(ANALYSIS_SCOPE==="domain"){document.getElementById("scopefull").style.display="none";document.getElementById("fullModeTabs").style.display="none";setAtlasMode("domains");}else setAtlasMode("clusters");
var curFam=null,glviewer=null,structMode="quality",repMode="cartoon",viewerBackground="white",selMembers={},curTree=null,pockMethod="fpocket";

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
 h+='<table>'+row("Members",d.n)+row("Foldseek TM · all pairs",fnum(d.tm))+row("Foldseek TM · retained edges",fnum(d.retained_tm))+row("BLAST best-HSP identity · all pairs",fnum(d.id_pct*100,1)+"%")+row("BLAST best-HSP identity · reported pairs",d.id_detected==null?"not detected":fnum(d.id_detected*100,1)+"%")+row("Max BLAST best-HSP identity",fnum(d.maxid*100,1)+"%")+row("Mean pLDDT",fnum(d.plddt,1))+'</table>';
 var domainLinks=((EXTRA[id]||{}).domain_families||[]).map(function(domain){return'<button class="domain-link" onclick="openLinkedDomain(\''+String(domain).replace(/'/g,"\\'")+'\')">'+esc(domain)+'</button>';}).join(" ");if(domainLinks)h+='<div class="domain-parent-links"><b>Domain families:</b> '+domainLinks+'</div>';
 var an=ANN[id];
 if(hasS){
  h+='<div class="tabs"><div class="tab on" onclick="tab(0)">Structure + Tree</div><div class="tab" onclick="tab(1)">FoldTree (figure)</div><div class="tab" onclick="tab(2)">Struct sim (TM)</div><div class="tab" onclick="tab(3)">Sequence + MSA</div><div class="tab" onclick="tab(4)">RNAseq</div><div class="tab" onclick="tab(5)">Annotation</div></div>';
  h+='<div id="p0" class="pane on"></div><div id="p1" class="pane"></div><div id="p2" class="pane"></div><div id="p3" class="pane"></div><div id="p4" class="pane"></div><div id="p5" class="pane"></div>';
 } else {
  h+='<div class="tabs"><div class="tab on" onclick="tab(5)">Annotation</div></div><div id="p5" class="pane on"></div>';
  h+='<p class="hint" style="margin-top:6px">Structure/tree/matrix layers computed for demo families F0 &amp; F12; annotation below covers all 39 families.</p>';
 }
 document.getElementById("side").innerHTML=h;
 if(hasS){ buildStructPane(id); setTimeout(function(){initViewer();var treebox=document.getElementById("treebox");if(PAY[id].newick){renderTree(parseNewick(PAY[id].newick),treebox);}else if(treebox){treebox.innerHTML='<p class="hint">FoldTree was not available for this family.</p>';}},50); } else { document.getElementById("p5").innerHTML=annHTML(id); }
}
function openLinkedDomain(id){setAtlasMode("domains");ensureDomainNetwork();if(domainNetwork){domainNetwork.selectNodes([id]);domainNetwork.focus(id,{scale:1.1,animation:true});}showDomain(id);}
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
   row("Pocket",s.pocket?esc(s.pocket_method||"detected")+" \u00b7 "+(s.pocket_metric==="probability"?"prob ":"score ")+fnum(s.pocket_value,3):"not detected")+
   row("RNA-seq peak",s.rna_condition?esc(s.rna_condition)+" \u00b7 "+fnum(s.rna_peak,2):"not available")+'</table>'+
   '<div class="tabs"><div class="tab on" onclick="singletonTab(0)">Structure</div><div class="tab" onclick="singletonTab(1)">RNA-seq</div><div class="tab" onclick="singletonTab(2)">Annotation</div><div class="tab" onclick="singletonTab(3)">Sequence</div></div>'+
   '<div id="sp0" class="pane on"></div><div id="sp1" class="pane"></div><div id="sp2" class="pane"></div><div id="sp3" class="pane"></div>';
 var linkedDomains=(ex.domain_families||[]).map(function(domain){return'<button class="domain-link" onclick="openLinkedDomain(\''+String(domain).replace(/'/g,"\\'")+'\')">'+esc(domain)+'</button>';}).join(" ");if(linkedDomains)h=h.replace('<div class="tabs">','<div class="domain-parent-links"><b>Domain families:</b> '+linkedDomains+'</div><div class="tabs">');
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
  if(i===1&&!document.getElementById("p1").innerHTML){var status=(EXTRA[curFam]||{}).sequence_analysis_status||{},treeHtml='<h3>FoldTree structural relationship tree</h3>'+(a.tree_svg?'<img src="'+a.tree_svg+'">':'<p class="hint">FoldTree was not available.</p>')+'<div class="hint">Built from structural distances; '+((EXTRA[curFam]&&EXTRA[curFam].foldtree_rooting_label)||"rooting status unavailable")+'. This is not a sequence phylogeny.</div>';treeHtml+='<h3 style="margin-top:14px">Sequence relationship tree</h3>'+(a.sequence_tree_svg?'<img src="'+a.sequence_tree_svg+'"><div class="hint">MAFFT L-INS-i homologous subgroup \u2192 FastTree WAG. Interpret as an exploratory sequence relationship tree.</div>':'<p class="hint">Not available: '+esc(status.reason||status.tree_status||"no sufficiently large homologous sequence subgroup")+'.</p>');document.getElementById("p1").innerHTML=treeHtml+dlbtn();}
  if(i===2&&!document.getElementById("p2").innerHTML){var ex2=EXTRA[curFam]||{};var h2='<h3>Structural similarity (Foldseek TM)</h3>'+(a.tm_svg?'<img src="'+a.tm_svg+'">':'<p class="hint">Foldseek TM matrix was not available for this family.</p>')+dlbtn();if(a.tmus_svg){var cr=(ex2.tm_cons_r!=null)?(' \u00b7 Foldseek\u2194US-align r='+fnum(ex2.tm_cons_r,3)+' over '+esc(ex2.tm_cons_n_pairs||0)+' mutually measured pair(s), max\u0394='+fnum(ex2.tm_cons_maxdiff,2)+', '+(ex2.tm_disagree||0)+' pair(s) disagree >0.1'):'';h2+='<h3 style="margin-top:14px">Structural similarity (US-align TM \u2014 independent algorithm)</h3><div class="hint" style="margin:2px 0 6px">Foldseek builds the families; US-align (TM-align successor) recomputes true TM within the family as an algorithm-independent cross-check. Missing Foldseek pairs are excluded from agreement statistics.'+cr+'</div><img src="'+a.tmus_svg+'">';}else{h2+='<h3 style="margin-top:14px">Structural similarity (US-align TM)</h3><p class="hint">Independent US-align matrix was not available.</p>';}document.getElementById("p2").innerHTML=h2;}
  if(i===3&&!document.getElementById("p3").innerHTML){document.getElementById("p3").innerHTML=(a.id_svg?'<h3>Sequence identity (BLASTp)</h3><img src="'+a.id_svg+'">':'<p class="hint">BLASTp identity matrix was not available for this cluster.</p>')+'<div id="familySequencePane"></div>'+dlbtn();buildSequencePane(curFam,true,"familySequencePane");}
  if(i===4&&!document.getElementById("p4").innerHTML)document.getElementById("p4").innerHTML=a.rna_svg?'<h3>RNA-seq expression</h3><img src="'+a.rna_svg+'">'+dlbtn():'<p class="hint">RNA-seq data were not available for this family.</p>'+dlbtn();
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
function backendJobId(){try{return new URLSearchParams(window.location.search).get("id")||"";}catch(e){return"";}}
function artifactUrl(kind,name){return"/artifact?id="+encodeURIComponent(backendJobId())+"&kind="+encodeURIComponent(kind)+"&name="+encodeURIComponent(name);}
var structureFetchFailed={};
function fetchStructure(acc,callback){var workspace=singletonById(acc)?acc:fullFamilyForAcc(acc),p=PAY[workspace]||{};p.struct=p.struct||{};if(p.struct[acc]){callback(p.struct[acc]);return;}if(!BACKEND.enabled||structureFetchFailed[acc]){callback(null);return;}fetch(artifactUrl("structure",acc)).then(function(r){if(!r.ok)throw new Error("structure "+r.status);return r.text();}).then(function(text){p.struct[acc]=text;callback(text);}).catch(function(){structureFetchFailed[acc]=true;callback(null);});}
var referenceFetchFailed={};
function hasReference(key){return!!REFPDB[key]||!!REFAVAIL[key];}
function fetchReference(key,callback){if(REFPDB[key]){callback(REFPDB[key]);return;}if(!BACKEND.enabled||referenceFetchFailed[key]||!REFAVAIL[key]){callback(null);return;}fetch(artifactUrl("reference",key)).then(function(r){if(!r.ok)throw new Error("reference "+r.status);return r.text();}).then(function(text){REFPDB[key]=text;callback(text);}).catch(function(){referenceFetchFailed[key]=true;callback(null);});}
function ensureFamilyStructures(fam,callback){var members=(PAY[fam]||{}).members||[],pending=members.length;if(!pending){callback();return;}members.forEach(function(m){fetchStructure(m,function(){pending--;if(!pending)callback();});});}
function basePdb(fam){var p=PAY[fam]||{},member=(p.members||[])[0];return REFPDB[fam+"_base"]||REFPDB[fam+"_cons"]||(p.struct||{})[fam]||(p.struct||{})[member]||"";}
function pdbWithValues(pdb,values,def){values=values||{};def=def==null?0:def;return pdb.split("\n").map(function(l){if((l.substring(0,4)==="ATOM"||l.substring(0,6)==="HETATM")&&l.length>=66){var ri=parseInt(l.substring(22,26)),v=Object.prototype.hasOwnProperty.call(values,ri)?Number(values[ri]):def;return l.substring(0,60)+v.toFixed(2).padStart(6," ")+l.substring(66);}return l;}).join("\n");}
function esmPdb(fam){var d=EXTRA[fam]||{};return REFPDB[fam+"_esm"]||(d.esm_values?pdbWithValues(basePdb(fam),d.esm_values,0):"");}
function alignedPdb(fam,m){var pdb=(PAY[fam].struct||{})[m],tr=(PAY[fam].transforms||{})[m];if(!pdb||!tr)return null;var r=tr.rotation,t=tr.translation;return pdb.split("\n").map(function(l){if((l.substring(0,4)==="ATOM"||l.substring(0,6)==="HETATM")&&l.length>=54){var x=parseFloat(l.substring(30,38)),y=parseFloat(l.substring(38,46)),z=parseFloat(l.substring(46,54));if(isFinite(x)&&isFinite(y)&&isFinite(z)){var nx=x*r[0][0]+y*r[1][0]+z*r[2][0]+t[0],ny=x*r[0][1]+y*r[1][1]+z*r[2][1]+t[1],nz=x*r[0][2]+y*r[1][2]+z*r[2][2]+t[2];return l.substring(0,30)+nx.toFixed(3).padStart(8," ")+ny.toFixed(3).padStart(8," ")+nz.toFixed(3).padStart(8," ")+l.substring(54);}}return l;}).join("\n");}
function dlStruct(kind){var fam=curFam,d=EXTRA[fam];
 if(kind==="quality"){var base=basePdb(fam),member=(PAY[fam].members||[])[0];if(!base&&BACKEND.enabled&&member&&!structureFetchFailed[member]){fetchStructure(member,function(){dlStruct(kind);});return;}if(!base){alert("No structure for "+fam);return;}dlText(base,fam+"_AlphaFold_structure.pdb");return;}
 if(kind==="structcons"){var sk=fam+"_struct";if(!hasReference(sk)){alert("Structural conservation was not available for "+fam);return;}if(!REFPDB[sk]&&BACKEND.enabled&&!referenceFetchFailed[sk]){fetchReference(sk,function(){dlStruct(kind);});return;}if(REFPDB[sk])dlText(REFPDB[sk],fam+"_structural_conservation.pdb");return;}
 if(kind==="evolcons"){var ck=fam+"_cons";if(!hasReference(ck)){alert("Evolutionary conservation was not applicable or unavailable for "+fam);return;}if(!REFPDB[ck]&&BACKEND.enabled&&!referenceFetchFailed[ck]){fetchReference(ck,function(){dlStruct(kind);});return;}if(REFPDB[ck])dlText(REFPDB[ck],fam+"_evolutionary_conservation.pdb");return;}
 if(kind==="esm"){var ek=fam+"_esm";if(!d.has_esm){alert("No ESM scan for "+fam);return;}if(hasReference(ek)&&!REFPDB[ek]&&BACKEND.enabled&&!referenceFetchFailed[ek]){fetchReference(ek,function(){dlStruct(kind);});return;}var epdb=esmPdb(fam);if(!epdb){alert("ESM structure unavailable for "+fam);return;}dlText(epdb,fam+"_ESM_tolerance.pdb");return;}
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
function dlXlsx(){var fam=curFam,b64=PAY[fam].assets.xlsx_b64;if(!b64&&BACKEND.enabled){window.location.href=artifactUrl("xlsx",fam);return;}var blob=b64toBlob(b64,"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet");var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fam+"_data.xlsx";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function fastaText(records,members){var out=[],n=0;(members||[]).forEach(function(m){var s=records[m];if(s){out.push(">"+m);for(var i=0;i<s.length;i+=60)out.push(s.substring(i,i+60));n++;}});return{text:out.join("\n")+(out.length?"\n":""),count:n};}
function dlSeqs(){var fam=curFam,p=PAY[fam],fa=fastaText(p.seq||{},p.members);if(!fa.count){alert("No sequences available for "+fam);return;}dlText(fa.text,fam+"_members_"+fa.count+"seqs.fasta","text/plain");}
function alignmentRecords(kind){var p=PAY[curFam]||{};if(kind==="sequence")return p.sequence_msa||{};if(kind==="threed")return p.three_di_msa||{};return p.structural_msa||p.msa||{};}
function dlAlignment(kind){var fam=curFam,p=PAY[fam],records=alignmentRecords(kind),members=Object.keys(records),fa=fastaText(records,members),names={sequence:"MAFFT_sequence_MSA",structural:"FoldMason_AA_MSA",threed:"FoldMason_3Di_MSA"};if(fa.count<2){alert("This alignment is not applicable or unavailable for "+fam);return;}dlText(fa.text,fam+"_"+names[kind]+"_"+fa.count+"seqs.fasta","text/plain");}
function dlMsa(){dlAlignment("structural");}
function dlAllStruct(){var fam=curFam,b64=PAY[fam].assets.structures_zip_b64;if(!b64&&BACKEND.enabled){window.location.href=artifactUrl("structures",fam);return;}if(!b64){alert("No structures available for "+fam);return;}var blob=b64toBlob(b64,"application/zip");var url=URL.createObjectURL(blob);var a=document.createElement("a");a.href=url;a.download=fam+"_member_structures.zip";document.body.appendChild(a);a.click();setTimeout(function(){document.body.removeChild(a);URL.revokeObjectURL(url);},1500);}
function dlMemberStruct(){var fam=curFam,sel=document.getElementById("memSel"),m=sel?sel.value:(PAY[fam].kind==="singleton"?fam:"");var st=(PAY[fam].struct||{})[m];if(!st&&BACKEND.enabled){window.location.href=artifactUrl("structure",m);return;}if(!st){alert("No structure for "+(m||fam));return;}dlText(st,m+".pdb");}
function dlMemberSeq(){var fam=curFam,sel=document.getElementById("memSel"),m=sel?sel.value:(PAY[fam].kind==="singleton"?fam:"");var s=(PAY[fam].seq||{})[m];if(!s){alert("No sequence for "+(m||fam));return;}var fa=fastaText((function(){var x={};x[m]=s;return x;})(),[m]);dlText(fa.text,m+".fasta","text/plain");}
function memberDlBar(fam){var p=PAY[fam],st=p.struct||{},sq=p.seq||{},sm=p.structural_msa||p.msa||{},qm=p.sequence_msa||{},di=p.three_di_msa||{},mem=p.members;var ns=BACKEND.enabled?mem.length:0,nq=0;mem.forEach(function(m){if(!BACKEND.enabled&&st[m])ns++;if(sq[m])nq++;});var opts=mem.map(function(m){return '<option value="'+m+'">'+m+'</option>';}).join("");return '<div class="dlbar">'+'<b>Download members</b> &middot; <span class="hint">'+nq+' sequences, '+ns+' structures</span><br>'+'<button class="dl" onclick="dlSeqs()">\u2b07 All sequences (FASTA)</button> '+(Object.keys(qm).length>1?'<button class="dl" onclick="dlAlignment(\'sequence\')">\u2b07 MAFFT sequence MSA</button> ':"")+(Object.keys(sm).length>1?'<button class="dl" onclick="dlAlignment(\'structural\')">\u2b07 FoldMason AA MSA</button> ':"")+(Object.keys(di).length>1?'<button class="dl" onclick="dlAlignment(\'threed\')">\u2b07 FoldMason 3Di MSA</button> ':"")+'<button class="dl" onclick="dlAllStruct()">\u2b07 All structures (ZIP)</button>'+'<br><span class="hint">single member:</span> <select id="memSel">'+opts+'</select> '+'<button class="dl" onclick="dlMemberSeq()">\u2b07 .fasta</button> '+'<button class="dl" onclick="dlMemberStruct()">\u2b07 .pdb</button>'+'</div>';}
function dlbtn(){return '<br><button class="dl" onclick="dlXlsx()">\u2b07 Download all '+curFam+' data (Excel: pockets / FoldTree / Foldseek / US-align / sequence / RNA-seq / per-site)</button>';}
function singletonDlbtn(){return '<div class="dlbar"><b>Download singleton</b><br><button class="dl" onclick="dlMemberSeq()">\u2b07 Sequence (FASTA)</button><button class="dl" onclick="dlMemberStruct()">\u2b07 Structure (PDB)</button><button class="dl" onclick="dlXlsx()">\u2b07 All evidence (Excel)</button></div>';}
var seqViewMode="member";
function sequenceMsaText(msa,members){var present=(members||[]).filter(function(m){return !!msa[m];});if(!present.length)return"No alignment available.";var nlen=Math.min(18,Math.max.apply(null,present.map(function(m){return m.length;}))),alen=Math.max.apply(null,present.map(function(m){return msa[m].length;})),width=60,out=[];for(var pos=0;pos<alen;pos+=width){out.push("Alignment columns "+(pos+1)+"-"+Math.min(pos+width,alen));present.forEach(function(m){var label=(m.length>nlen?m.slice(0,nlen):m).padEnd(nlen," ");out.push(label+"  "+msa[m].substring(pos,pos+width));});out.push("");}return out.join("\n");}
function setSequenceView(mode){seqViewMode=mode;["member","sequence","structural","threed"].forEach(function(x){var b=document.getElementById("seqmode_"+x);if(b)b.className=(x===mode?"on":"");});var sel=document.getElementById("seqMember");if(sel)sel.style.display=mode==="member"?"inline-block":"none";renderSequenceViewer();}
function renderSequenceViewer(){var p=PAY[curFam]||{},box=document.getElementById("sequenceViewer");if(!box)return;if(seqViewMode!=="member"){var records=alignmentRecords(seqViewMode),members=Object.keys(records);box.textContent=sequenceMsaText(records,members);return;}var sel=document.getElementById("seqMember"),m=sel?sel.value:(p.members||[])[0],s=(p.seq||{})[m]||"";box.textContent=s?(">"+m+" | "+s.length+" aa\n"+fastaText((function(){var x={};x[m]=s;return x;})(),[m]).text.split("\n").slice(1).join("\n")):"Sequence unavailable.";}
function buildSequencePane(fam,allowMsa,targetId){var p=PAY[fam]||{},mem=p.members||[],sequence=p.sequence_msa||{},structural=p.structural_msa||p.msa||{},threed=p.three_di_msa||{},target=document.getElementById(targetId);if(!target)return;seqViewMode="member";var opts=mem.map(function(m){return'<option value="'+esc(m)+'">'+esc(m)+'</option>';}).join(""),downloadFn=p.kind==="singleton"?"dlMemberSeq()":"dlSeqs()",status=(EXTRA[fam]||{}).sequence_analysis_status||{},buttons="";if(allowMsa&&Object.keys(sequence).length>1)buttons+='<button id="seqmode_sequence" onclick="setSequenceView(\'sequence\')">Sequence MSA</button>';if(allowMsa&&Object.keys(structural).length>1)buttons+='<button id="seqmode_structural" onclick="setSequenceView(\'structural\')">Structural MSA (AA)</button>';if(allowMsa&&Object.keys(threed).length>1)buttons+='<button id="seqmode_threed" onclick="setSequenceView(\'threed\')">Structural MSA (3Di)</button>';var note=allowMsa?'<div class="hint">Sequence MSA: MAFFT on the hub\u2019s BLAST-homologous subgroup ('+esc(status.n_sequences||0)+' sequences; '+esc(status.rate4site_status||"not run")+'). Structural MSA: FoldMason across the full structural family.</div>':'<div class="hint">Multiple alignments are not applicable to a singleton.</div>';target.innerHTML='<h3>Sequence and alignment viewer</h3><div class="sequence-toolbar"><button id="seqmode_member" class="on" onclick="setSequenceView(\'member\')">Member</button>'+buttons+'<select id="seqMember" onchange="renderSequenceViewer()">'+opts+'</select><button onclick="'+downloadFn+'">\u2b07 '+(p.kind==="singleton"?"Sequence":"All sequences")+' (FASTA)</button></div>'+note+'<pre id="sequenceViewer" class="sequence-view"></pre>';renderSequenceViewer();}
function buildStructPane(id){
 var ex=EXTRA[id],mem=PAY[id].members;pockMethod="fpocket";selMembers={};mem.forEach(function(m){selMembers[m]=true;});
 var hasStructural=hasReference(id+"_struct"),hasEvolutionary=hasReference(id+"_cons")&&ex.cons_min!=null&&ex.cons_max!=null&&(ex.sequence_scored_resi||[]).length>0;
 var h='<div><button id="bquality" class="on" onclick="setMode(\'quality\')">pLDDT</button>'+(hasStructural?'<button id="bstructcons" onclick="setMode(\'structcons\')">Structural conservation</button>':"")+
   (hasEvolutionary?'<button id="bevolcons" onclick="setMode(\'evolcons\')">Evolutionary conservation</button>':"")+
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
   '<button id="r_line" onclick="setRep(\'line\')">Line</button><span class="hint" style="margin-left:7px">Background:</span><button id="bg_white" class="on" onclick="setViewerBackground(\'white\')">White</button><button id="bg_black" onclick="setViewerBackground(\'black\')">Black</button></div>'+
   '<div id="v3d"></div><div id="leg" class="hint"></div>'+'<div style="margin-top:6px"><span class="hint">Download structure (ChimeraX/PyMOL):</span><br>'+(hasStructural?'<button onclick="dlStruct(\'structcons\')">Structural conservation PDB</button>':"")+(hasEvolutionary?'<button onclick="dlStruct(\'evolcons\')">Evolutionary conservation PDB</button>':"")+'<button onclick="dlStruct(\'esm\')">ESM PDB</button>'+'<button onclick="dlStruct(\'pocket\')" title="pocket residues in B-factor column (999=pocket)">Pocket-annotated PDB</button>'+'<button onclick="dlStruct(\'super\')">Superpose selected PDB</button></div>';
 h+=memberDlBar(id);
 h+='<h3>FoldTree \u2014 click tips or internal nodes to pick members</h3>'+'<div class="hint" style="margin:2px 0 6px">\u2605 <b>gold star</b> = family hub (best Foldseek pair coverage, then highest mean measured TM '+fnum(EXTRA[curFam]&&EXTRA[curFam].hub_meanTM,2)+'; '+esc((EXTRA[curFam]&&EXTRA[curFam].hub_pairs_measured)||0)+'/'+esc((EXTRA[curFam]&&EXTRA[curFam].hub_pairs_expected)||0)+' pairs measured). '+'The hub is the most representative measured fold; pocket/conservation are currently projected on '+(EXTRA[curFam]?EXTRA[curFam].ref_used:"?")+' (see note).</div>'+
    '<div style="margin:3px 0"><button onclick="allMem(1)">All</button><button onclick="allMem(0)">None</button>'+
    '<span class="hint" style="margin-left:8px">tip = one member &middot; orange node = whole clade</span></div>'+
    '<div id="treebox"></div>';
 function reslist(a){return (a&&a.length)?a.join(", "):"—";}
 h+='<h3>Selection &amp; sites</h3><table>'+row("Structural conservation (mean)",fnum(ex.structural_lddt_mean,3))+row("Evolutionary analysis",(ex.sequence_analysis_status||{}).rate4site_status||"not run")+row("Sequence subgroup",(ex.sequence_analysis_status||{}).sequence_subgroup||"not applicable")+row("fpocket",(ex.fpocket_score!=null?"score "+fnum(ex.fpocket_score,3)+", ":"")+ex.fpocket_resi.length+" res")+row("P2Rank",ex.p2rank_resi.length?("prob "+fnum(ex.p2rank_prob,3)+", "+ex.p2rank_resi.length+" res"):"no pocket")+row("Cys anchors",ex.n_cys)+row("Conserved-buried r",fnum(ex.cons_sasa_r,2))+row("ESM vs evolutionary conservation r",ex.esm_vs_cons_r!=null?fnum(ex.esm_vs_cons_r,2):"n/a")+row("ESM vs SASA r",ex.esm_vs_sasa_r!=null?fnum(ex.esm_vs_sasa_r,2):"n/a")+(ex.tm_us_mean!=null?row("US-align TM (mean)",fnum(ex.tm_us_mean,3))+row("Foldseek\u2194US-align r",ex.tm_cons_r!=null?fnum(ex.tm_cons_r,3):"n/a"):"")+'</table>';
 h+='<div class="hint" style="margin-top:4px"><b>fpocket lining residues</b> (ref '+(ex.ref_used||"?")+'): <span style="color:#233">'+reslist(ex.fpocket_resi)+'</span></div>';
 if(ex.p2rank_resi&&ex.p2rank_resi.length)h+='<div class="hint"><b>P2Rank lining residues</b>: <span style="color:#233">'+reslist(ex.p2rank_resi)+'</span></div>';
 h+='<div style="margin-top:5px"><button onclick="dlPockResidues()">\u2b07 Pocket residues (CSV)</button></div>';
 document.getElementById("p0").innerHTML=h;
}
function buildSingletonStructPane(id){
 var ex=EXTRA[id]||{},p=PAY[id],hasStruct=BACKEND.enabled||!!((p.struct||{})[id]||REFPDB[id+"_base"]);
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
   '<button id="r_line" onclick="setRep(\'line\')">Line</button><span class="hint" style="margin-left:7px">Background:</span><button id="bg_white" class="on" onclick="setViewerBackground(\'white\')">White</button><button id="bg_black" onclick="setViewerBackground(\'black\')">Black</button></div>'+
   (hasStruct?'<div id="v3d"></div><div id="leg" class="hint"></div>':'<p class="hint">No structure was embedded for this protein.</p>')+
   '<div style="margin-top:6px"><span class="hint">Download structure:</span><br>'+
   '<button onclick="dlStruct(\'quality\')">Original PDB</button>'+
   '<button onclick="dlStruct(\'esm\')">ESM PDB</button>'+
   '<button onclick="dlStruct(\'pocket\')">Pocket-annotated PDB</button></div>';
 h+=singletonDlbtn();
 function reslist(a){return(a&&a.length)?a.join(", "):"\u2013";}
 h+='<h3>Pocket and sequence evidence</h3><table>'+
   row("fpocket",(ex.fpocket_resi||[]).length?("score "+fnum(ex.fpocket_score,3)+" \u00b7 "+ex.fpocket_resi.length+" residues"):"no pocket")+
   row("P2Rank",(ex.p2rank_resi||[]).length?("prob "+fnum(ex.p2rank_prob,3)+" \u00b7 "+ex.p2rank_resi.length+" residues"):"no pocket")+
   row("Cysteines",ex.n_cys||0)+row("ESM scan",ex.has_esm?"available":"not available")+'</table>'+
   '<div class="hint"><b>fpocket lining residues:</b> '+reslist(ex.fpocket_resi)+'</div>'+
   '<div class="hint"><b>P2Rank lining residues:</b> '+reslist(ex.p2rank_resi)+'</div>'+
   '<div style="margin-top:5px"><button onclick="dlPockResidues()">\u2b07 Pocket residues (CSV)</button></div>';
 document.getElementById("sp0").innerHTML=h;
}
function allMem(v){for(var m in selMembers)selMembers[m]=!!v;paintTree();if(structMode==="super")drawStruct();}
function initViewer(){var el=document.getElementById("v3d");if(!el)return;var member=(PAY[curFam].members||[])[0],preferred="quality",key="";if(BACKEND.enabled&&!basePdb(curFam)&&member&&!structureFetchFailed[member]){el.innerHTML='<p class="hint" style="padding:12px">Loading structure…</p>';fetchStructure(member,function(){initViewer();});return;}el.innerHTML="";glviewer=$3Dmol.createViewer(el,{backgroundColor:viewerBackground});structMode=preferred;repMode="cartoon";setViewerBackground(viewerBackground);setMode(structMode);}
function setPock(x){pockMethod=x;["fpocket","p2rank"].forEach(function(y){var b=document.getElementById("pk_"+y);if(b)b.className=(y===x?"on":"");});if(structMode==="pocket")drawStruct();}
function setMode(m){structMode=m;var pr=document.getElementById("pockrow");if(pr)pr.style.display=(m==="pocket")?"block":"none";["quality","structcons","evolcons","pocket","esm","super"].forEach(function(x){var b=document.getElementById("b"+x);if(b)b.className=(x===m?"on":"");});drawStruct();}
function setRep(r){repMode=r;["cartoon","surface","stick","sphere","line"].forEach(function(x){var b=document.getElementById("r_"+x);if(b)b.className=(x===r?"on":"");});drawStruct();}
function setViewerBackground(color){viewerBackground=color;["white","black"].forEach(function(value){var button=document.getElementById("bg_"+value);if(button)button.className=value===color?"on":"";});if(glviewer){glviewer.setBackgroundColor(color);glviewer.render();}}
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
  var missing=BACKEND.enabled&&PAY[curFam].members.some(function(m){return selMembers[m]&&!(PAY[curFam].struct||{})[m]&&!structureFetchFailed[m];});
  if(missing){document.getElementById("leg").innerHTML="Loading selected family structures…";ensureFamilyStructures(curFam,function(){drawStruct();});return;}
  var cols=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];
  var ci=0,shown=0,mem=PAY[curFam].members;
  mem.forEach(function(m){var pdb=alignedPdb(curFam,m);if(selMembers[m]&&pdb){glviewer.addModel(pdb,"pdb");applyStyle({model:ci},{color:cols[ci%cols.length]});ci++;shown++;}});
  glviewer.zoomTo();glviewer.render();
  document.getElementById("leg").innerHTML="Hub-referenced FoldMason/Kabsch superposition of "+shown+" selected members ("+repMode+"). Tight core = conserved scaffold; splayed loops = variable surface.";
  return;
 }
 var refKey=structMode==="structcons"?curFam+"_struct":(structMode==="evolcons"?curFam+"_cons":(structMode==="esm"&&d.has_esm?curFam+"_esm":""));
 if(BACKEND.enabled&&refKey&&hasReference(refKey)&&!REFPDB[refKey]&&!referenceFetchFailed[refKey]){document.getElementById("leg").innerHTML="Loading structure data…";fetchReference(refKey,function(){drawStruct();});return;}
 var baseMember=(PAY[curFam].members||[])[0],needsBase=structMode==="quality"||structMode==="pocket"||(structMode==="esm"&&!hasReference(refKey));
 if(BACKEND.enabled&&needsBase&&!basePdb(curFam)&&baseMember&&!structureFetchFailed[baseMember]){document.getElementById("leg").innerHTML="Loading structure data…";fetchStructure(baseMember,function(){drawStruct();});return;}
 var pdbtext=(structMode==="esm"&&d.has_esm)?esmPdb(curFam):(structMode==="quality"?basePdb(curFam):(structMode==="structcons"?(REFPDB[curFam+"_struct"]||basePdb(curFam)):(REFPDB[curFam+"_cons"]||basePdb(curFam))));
 if(!pdbtext){document.getElementById("leg").innerHTML="Structure unavailable.";return;}
 glviewer.addModel(pdbtext,"pdb");
 if(structMode==="esm"){
  if(d.has_esm){applyStyle({},{colorscheme:{prop:"b",gradient:"rwb",min:d.esm_max,max:d.esm_min}});
   document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>constrained <span class="swatch" style="background:#b2182b"></span>tolerant &middot; ESM-1b '+fnum(d.esm_min,1)+"\u2026"+fnum(d.esm_max,1)+" &middot; red = mutation-tolerant (variable), blue = constrained";}
  else{applyStyle({},{color:"lightgrey"});document.getElementById("leg").innerHTML="ESM scan unavailable for this family.";}
 }else if(structMode==="quality"){
  var qmin=d.plddt_min==null?50:d.plddt_min,qmax=d.plddt_max==null?100:d.plddt_max;
  applyStyle({},{colorscheme:{prop:"b",gradient:"rwb",min:qmax,max:qmin}});
  document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>lower confidence <span class="swatch" style="background:#b2182b"></span>higher confidence &middot; AlphaFold pLDDT '+fnum(qmin,1)+"\u2026"+fnum(qmax,1);
 }else if(structMode==="structcons"){
  applyStyle({},{color:"lightgrey"});applyStyle({resi:d.structural_scored_resi||[]},{colorscheme:{prop:"b",gradient:"rwb",min:100,max:0}});
  document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>structurally variable <span class="swatch" style="background:#b2182b"></span>structurally conserved <span class="swatch" style="background:lightgrey"></span>insufficient pair coverage &middot; official FoldMason column LDDT (pair threshold '+fnum(d.structural_pair_threshold==null ? .5 : d.structural_pair_threshold,2)+')';
 }else if(structMode==="evolcons"){
  applyStyle({},{color:"lightgrey"});applyStyle({resi:d.sequence_scored_resi||[]},{colorscheme:{prop:"b",gradient:"rwb",min:d.cons_max,max:d.cons_min}});
  document.getElementById("leg").innerHTML='<span class="swatch" style="background:#2166ac"></span>evolutionarily variable <span class="swatch" style="background:#b2182b"></span>evolutionarily conserved <span class="swatch" style="background:lightgrey"></span>not scored &middot; Rate4Site on sequence-homologous subgroup '+fnum(d.cons_min,1)+"\u2026"+fnum(d.cons_max,1);
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
function openAtlasTargetFromUrl(){
 var params;try{params=new URLSearchParams(window.location.search);}catch(e){return;}
 var protein=params.get("protein"),target=params.get("open"),segment=params.get("segment");
 if(protein){openProtein(protein);return;}
 if(!target)return;
 if(/^D\d+$/i.test(target)&&DOMAIN_FAMILIES.some(function(d){return String(d.domain_family).toLowerCase()===target.toLowerCase();})){
  var domainId=DOMAIN_FAMILIES.find(function(d){return String(d.domain_family).toLowerCase()===target.toLowerCase();}).domain_family;
  setAtlasMode("domains");
  setTimeout(function(){
   if(domainNetwork){domainNetwork.selectNodes([domainId]);domainNetwork.focus(domainId,{scale:1.15,animation:true});}
   showDomain(domainId,segment?("protein:"+segment):"");
  },80);
  return;
 }
 var familyId=Object.keys(PAY).find(function(id){return id.toLowerCase()===target.toLowerCase()&&PAY[id].kind!=="singleton";});
 if(familyId){
  setAtlasMode("clusters");
  network.selectNodes([familyId]);network.focus(familyId,{scale:1.15,animation:true});showFamily(familyId);
 }
}
setTimeout(openAtlasTargetFromUrl,0);
