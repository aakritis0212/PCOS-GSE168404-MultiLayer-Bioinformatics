import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px

st.set_page_config(page_title="PCOS Bioinformatics Analysis", page_icon="🧬", layout="wide")
ROOT = Path(__file__).resolve().parent

def find(name):
    x = list(ROOT.rglob(name))
    return x[0] if x else None

def load(name):
    p = find(name)
    if not p: return None
    try: return pd.read_csv(p)
    except: return None

def col(df, names):
    if df is None: return None
    return next((x for x in names if x in df.columns), None)

limma=load("PCOS_limma_results.csv")
gsea=load("PCOS_GSEA_significant_pathways.csv") or load("PCOS_GSEA_GO_BP_results.csv")
core=load("PCOS_core_genes_shared_pathways.csv")
regs=load("PCOS_FINAL_REGULATOR_MASTER_TABLE.csv")
sigregs=load("PCOS_FINAL_SIGNIFICANT_REGULATORS.csv")
edges=load("PCOS_significant_regulator_gene_network_edges.csv")

st.sidebar.title("🧬 PCOS Dashboard")
page=st.sidebar.radio("Navigate",["🏠 Overview","📊 Differential Expression","🧬 GSEA","🔬 Core Genes","🧠 Regulators","🕸️ Network"])

if page=="🏠 Overview":
    st.title("🧬 PCOS Transcriptomic Analysis")
    st.write("Interactive dashboard for GSE168404: differential expression, GO-BP GSEA, leading-edge/core genes and regulator enrichment.")
    a,b,c,d=st.columns(4)
    a.metric("Samples","10")
    b.metric("Genes analysed","14,712")
    c.metric("Significant GSEA pathways","446")
    d.metric("Core genes","1,774")
    st.subheader("Significant regulators")
    if sigregs is not None: st.dataframe(sigregs,use_container_width=True,hide_index=True)
    else: st.info("Regulator table not found.")
    st.info("Interpretation: enrichment supports candidate regulatory hypotheses and does not establish causality.")

elif page=="📊 Differential Expression":
    st.title("📊 Differential Expression")
    if limma is None: st.error("PCOS_limma_results.csv not found.")
    else:
        g=col(limma,["GeneSymbol","Gene","Symbol"]); lf=col(limma,["logFC"]); f=col(limma,["adj.P.Val","FDR","padj"])
        q=st.text_input("Search gene")
        v=limma.copy()
        if q and g: v=v[v[g].astype(str).str.contains(q,case=False,na=False)]
        st.dataframe(v.head(300),use_container_width=True,hide_index=True)
        if lf and f and g:
            v[lf]=pd.to_numeric(v[lf],errors="coerce"); v[f]=pd.to_numeric(v[f],errors="coerce")
            v["-log10(FDR)"]=-np.log10(v[f].clip(lower=1e-300))
            v=v.replace([np.inf,-np.inf],np.nan).dropna(subset=[lf,"-log10(FDR)"])
            fig=px.scatter(v,x=lf,y="-log10(FDR)",hover_name=g,title="Differential-expression overview")
            st.plotly_chart(fig,use_container_width=True)

elif page=="🧬 GSEA":
    st.title("🧬 GO Biological Process GSEA")
    if gsea is None: st.error("GSEA results not found.")
    else:
        t=col(gsea,["Term","Description","Pathway"]); n=col(gsea,["NES"]); f=col(gsea,["padj","FDR","fdr"])
        if n: gsea[n]=pd.to_numeric(gsea[n],errors="coerce")
        st.metric("Pathways",f"{len(gsea):,}")
        if t and n:
            top=pd.concat([gsea.nlargest(20,n),gsea.nsmallest(20,n)]).drop_duplicates().sort_values(n)
            fig=px.bar(top,x=n,y=t,orientation="h",title="Top enriched pathways",height=800)
            st.plotly_chart(fig,use_container_width=True)
        st.dataframe(gsea.head(300),use_container_width=True,hide_index=True)

elif page=="🔬 Core Genes":
    st.title("🔬 Core Genes")
    if core is None: st.error("Core-gene table not found.")
    else:
        g=col(core,["Gene","GeneSymbol","Symbol"]); n=col(core,["Pathway_Count","PathwayCount","Count"])
        st.metric("Core genes",f"{len(core):,}")
        q=st.text_input("Search core gene")
        v=core.copy()
        if q and g: v=v[v[g].astype(str).str.contains(q,case=False,na=False)]
        st.dataframe(v.head(300),use_container_width=True,hide_index=True)
        if g and n:
            top=core.nlargest(30,n).sort_values(n)
            fig=px.bar(top,x=n,y=g,orientation="h",title="Top core genes by pathway participation",height=700)
            st.plotly_chart(fig,use_container_width=True)

elif page=="🧠 Regulators":
    st.title("🧠 Regulator Enrichment")
    if sigregs is not None:
        st.subheader("Significant regulators")
        st.dataframe(sigregs,use_container_width=True,hide_index=True)
    if regs is not None:
        r=col(regs,["Regulator","TF","GeneSymbol","Gene"]); f=col(regs,["FDR","padj","AdjustedPValue"]); o=col(regs,["Overlap","OverlapGenes","TargetOverlap"])
        if r and f:
            x=regs.copy(); x[f]=pd.to_numeric(x[f],errors="coerce"); x["-log10(FDR)"]=-np.log10(x[f].clip(lower=1e-300))
            x=x.replace([np.inf,-np.inf],np.nan).dropna(subset=["-log10(FDR)"]).nlargest(30,"-log10(FDR)")
            fig=px.bar(x.sort_values("-log10(FDR)"),x="-log10(FDR)",y=r,orientation="h",title="Top enriched regulators",height=700)
            st.plotly_chart(fig,use_container_width=True)
    st.caption("Regulator enrichment identifies over-representation of known regulator targets in the PCOS core-gene set; it does not prove causality.")

elif page=="🕸️ Network":
    st.title("🕸️ Significant Regulator–Target Network")
    if edges is None: st.warning("PCOS_significant_regulator_gene_network_edges.csv not found.")
    else:
        st.success(f"{len(edges):,} regulator-target associations loaded.")
        st.dataframe(edges.head(500),use_container_width=True,hide_index=True)
        s=col(edges,["Regulator","TF","Source"]); t=col(edges,["Gene","Target","TargetGene"])
        if s and t:
            counts=edges[s].value_counts().head(20)
            x=counts.reset_index(); x.columns=["Regulator","Targets"]
            fig=px.bar(x.sort_values("Targets"),x="Targets",y="Regulator",orientation="h",title="Regulators by number of connected targets",height=650)
            st.plotly_chart(fig,use_container_width=True)

st.markdown("---")
st.caption("PCOS Bioinformatics Analysis | GSE168404 | R-derived results presented through Streamlit.")
