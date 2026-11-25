#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correction d'effet batch de comptages RNA-seq avec pyComBat-seq (inmoose).

Entierement parametrable via le fichier de design et des options en ligne
de commande : aucun echantillon n'est code en dur, le script est reutilisable
pour toute autre experience.

Fichier de design (CSV ou TSV, une ligne par echantillon) :
  - colonne echantillons  : nom du fichier de comptages sans extension
  - colonne batch         : variable a corriger
  - colonnes covariables  : variables biologiques a PRESERVER

Exemples :

  docker run --rm -v "$(pwd)":/data pycombat-seq \
      python /data/batch_correction.py --design design_Huh75.csv

  docker run --rm -v "$(pwd)":/data pycombat-seq \
      python /data/batch_correction.py \
      --design design_Huh75.csv \
      --counts-dir /data/counts \
      --batch-col repetition \
      --covar-cols condition \
      --out-matrix Huh75_counts_corrected.tsv
"""

import argparse
import os
import sys

import pandas as pd
from inmoose.pycombat import pycombat_seq


def parse_args():
    p = argparse.ArgumentParser(
        description="Correction d'effet batch pyComBat-seq sur des comptages HTSeq.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--design", required=True,
                   help="Fichier de design (CSV ou TSV, une ligne par echantillon).")
    p.add_argument("--counts-dir", default=".",
                   help="Repertoire contenant les fichiers de comptages htseq-count.")
    p.add_argument("--sample-col", default="sample",
                   help="Colonne du design avec les noms d'echantillons "
                        "(noms de fichiers sans extension).")
    p.add_argument("--batch-col", default="repetition",
                   help="Colonne du design definissant les batches a corriger.")
    p.add_argument("--covar-cols", nargs="+", default=None,
                   help="Colonnes du design a conserver comme covariables "
                        "biologiques (ex. : condition). Defaut : toutes les "
                        "colonnes sauf echantillon et batch.")
    p.add_argument("--file-ext", default=".tsv",
                   help="Extension des fichiers de comptages.")
    p.add_argument("--counts-sep", default="\t",
                   help="Separateur des fichiers de comptages htseq-count.")
    p.add_argument("--design-sep", default=None,
                   help="Separateur du fichier de design "
                        "(defaut : detection automatique).")
    p.add_argument("--min-count", type=int, default=10,
                   help="Comptage minimal pour qu'un gene soit considere exprime.")
    p.add_argument("--min-samples", type=int, default=3,
                   help="Nombre minimal d'echantillons ou le gene doit etre exprime.")
    p.add_argument("--out-matrix", default="counts_corrected.tsv",
                   help="Fichier de sortie : comptages corriges (matrice genes x echantillons).")
    p.add_argument("--out-raw", default=None,
                   help="Fichier de sortie optionnel : comptages brutes apres filtre.")
    p.add_argument("--out-metadata", default=None,
                   help="Fichier de sortie optionnel : design avec la colonne "
                        "batch renommee (utile comme colData DESeq2). "
                        "Chaine vide pour desactiver.")
    p.add_argument("--htseq-dir", default="counts_corrected_htseq",
                   help="Repertoire de sortie : un fichier de comptages corrige "
                        "par echantillon, au format htseq-count (gene\\tcomptage, "
                        "sans en-tete). Chaine vide pour desactiver.")
    p.add_argument("--no-htseq-stats", action="store_true",
                   help="Ne pas recopier en fin de fichier les lignes de resume "
                        "htseq-count (__no_feature, __ambiguous, ...).")
    return p.parse_args()


def load_design(path, sep):
    """Charge le design (BOM gere, separateur detecte si non fourni)."""
    if not os.path.isfile(path):
        sys.exit(f"[ERREUR] Fichier de design introuvable : {path}")
    try:
        design = pd.read_csv(path, sep=sep if sep else None,
                             engine="python", encoding="utf-8-sig",
                             dtype=str).dropna(axis=0, how="all")
    except Exception as e:
        sys.exit(f"[ERREUR] Impossible de lire le design {path} : {e}")
    design.columns = [c.strip() for c in design.columns]
    return design


def main():
    args = parse_args()

    # ------------------------------------------------------------- design
    design = load_design(args.design, args.design_sep)

    for col in (args.sample_col, args.batch_col):
        if col not in design.columns:
            sys.exit(f"[ERREUR] Colonne '{col}' absente du design. "
                     f"Colonnes disponibles : {list(design.columns)}")

    samples = [str(s).strip() for s in design[args.sample_col]]
    if len(samples) != len(set(samples)):
        sys.exit("[ERREUR] Noms d'echantillons dupliques dans le design.")

    batches = design[args.batch_col].astype(str).str.strip()

    # Covariables biologiques : celles fournies, sinon toutes sauf batch
    if args.covar_cols:
        covar_cols = args.covar_cols
        missing = [c for c in covar_cols if c not in design.columns]
        if missing:
            sys.exit(f"[ERREUR] Colonnes de covariables absentes du design : {missing}")
    else:
        # Colonnes identifiantes (une valeur par echantillon, ex. 'prefix')
        # : a exclure des covariables car elles surajustent le modele.
        skip_ident = [c for c in design.columns
                      if c not in (args.sample_col, args.batch_col)
                      and design[c].nunique() == len(design)]
        if skip_ident:
            print(f"[INFO] Colonnes identifiantes ignorees (1 valeur par "
                  f"echantillon) : {skip_ident}")
        covar_cols = [c for c in design.columns
                      if c not in (args.sample_col, args.batch_col)
                      and c not in skip_ident]
        if not covar_cols:
            sys.exit("[ERREUR] Aucune covariable biologique disponible : "
                     "specifiez --covar-cols.")

    # Alerte si une covariable est confondue avec le batch
    for c in covar_cols:
        tab = pd.crosstab(design[c].astype(str), batches)
        if (tab > 0).sum(axis=1).max() == 1:
            print(f"[ATTENTION] La covariable '{c}' est confondue avec le "
                  f"batch : pyComBat-seq pourrait en supprimer l'effet.")

    # ------------------------------------------------------------- comptages
    counts = {}
    htseq_stats = {}     # lignes de resume htseq-count (__no_feature, ...), conservees
    for s in samples:
        fname = os.path.join(args.counts_dir, s + args.file_ext)
        if not os.path.isfile(fname):
            sys.exit(f"[ERREUR] Fichier de comptages introuvable : {fname}")
        df = pd.read_csv(fname, sep=args.counts_sep, header=None,
                         names=["gene", "count"], index_col=0,
                         comment="#", dtype=str)
        # Conserver les lignes de resume HTSeq puis les exclure de la matrice
        stats = df[df.index.astype(str).str.startswith("__")]
        htseq_stats[s] = list(zip(stats.index.astype(str), stats["count"]))
        df = df[~df.index.astype(str).str.startswith("__")]
        counts[s] = df["count"].astype(int)

    mat = pd.DataFrame(counts)[samples]           # ordre = ordre du design
    mat = mat.fillna(0).astype(int)              # genes absents d'un echantillon

    keep = (mat >= args.min_count).sum(axis=1) >= args.min_samples
    mat = mat[keep]

    print(f"Matrice : {mat.shape[0]} genes x {mat.shape[1]} echantillons")
    print(f"Batches ('{args.batch_col}') : "
          f"{dict(batches.value_counts().sort_index())}")
    print(f"Covariables conservees : {covar_cols}")
    if batches.nunique() < 2:
        sys.exit("[ERREUR] Il faut au moins 2 batches pour corriger un effet batch.")

    # ------------------------------------------------------------- correction
    # pyComBat-seq : on corrige le batch et on PRESERVE les covariables.
    covar = design.set_index(pd.Index(samples, name=args.sample_col))[[c for c in covar_cols]]
    for c in covar_cols:
        covar[c] = covar[c].astype(str).str.strip()

    # inmoose >= 0.9 : argument 'covar_mod' (matrice de design des covariables
    # categorielles, une ligne par echantillon) ; versions plus anciennes : 'covar'.
    try:
        adjusted = pycombat_seq(mat, list(batches), covar_mod=covar)
    except TypeError:
        adjusted = pycombat_seq(mat, list(batches), covar=covar)

    adj = pd.DataFrame(adjusted, index=mat.index, columns=mat.columns)
    adj = adj.clip(lower=0).round().astype(int)

    # ------------------------------------------------------------- sorties
    adj.to_csv(args.out_matrix, sep="\t", index=True, index_label="gene")
    print(f"Ecriture de {args.out_matrix}")
    if args.out_raw:
        mat.to_csv(args.out_raw, sep="\t", index=True, index_label="gene")
        print(f"Ecriture de {args.out_raw}")
    if args.out_metadata:
        meta = design.copy()
        meta.rename(columns={args.batch_col: "batch"}, inplace=True)
        meta.to_csv(args.out_metadata, sep="\t", index=True, index_label="sample")
        print(f"Ecriture de {args.out_metadata}")

    # ------------------------------------------------ sorties au format htseq
    # Un fichier par echantillon, identique en forme a la sortie de htseq-count :
    # deux colonnes sans en-tete (gene \t comptage), lignes __ en fin de fichier.
    if args.htseq_dir:
        os.makedirs(args.htseq_dir, exist_ok=True)
        for s in samples:
            outpath = os.path.join(args.htseq_dir, s + args.file_ext)
            with open(outpath, "w", encoding="utf-8") as fh:
                for gene, cnt in adj[s].items():
                    fh.write(f"{gene}\t{cnt}\n")
                if not args.no_htseq_stats:
                    for stat_name, stat_val in htseq_stats.get(s, []):
                        fh.write(f"{stat_name}\t{stat_val}\n")
        print(f"Ecriture de {len(samples)} fichiers htseq par echantillon dans "
              f"{args.htseq_dir}/")
    print("Termine.")


if __name__ == "__main__":
    main()