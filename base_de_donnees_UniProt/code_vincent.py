# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 07:51:50 2026

@author: vincent.durr
"""

import pandas as pd
import requests
import time
import json
from pathlib import Path

# =====================================================
# 1. Fichiers
# =====================================================
input_file = Path(r"C:\\Users\\vincent.durr\\Documents\\DTT red+ ox- no mix .xlsx")
output_file = Path(r"C:\\Users\\vincent.durr\\Documents\\uniprot_localisation_from_GeneID.xlsx")
json_dir = Path(r"C:\\Users\\vincent.durr\\Documents\\uniprot_json_localisation")

json_dir.mkdir(exist_ok=True)

# =====================================================
# 2. Lecture du fichier (DataFrame)
#    On utilise la colonne 'GeneID' comme identifiant de gène
# =====================================================
df = pd.read_excel(input_file)

if "GeneID" not in df.columns:
    raise ValueError(f"La colonne 'GeneID' n'existe pas dans le fichier. Colonnes disponibles : {list(df.columns)}")

gene_ids = df["GeneID"].dropna().astype(str).unique().tolist()
print(f"Nombre de GeneID uniques à interroger : {len(gene_ids)}")

# =====================================================
# 3. Requête UniProt pour récupérer la localisation
#    On suppose que GeneID correspond au symbole de gène humain
#    (sinon il faudra adapter le champ de requête : par ex. id: si c'est un UniProt ID)
# =====================================================
def query_uniprot_localisation(gene_ids, save_json=True):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    fields = [
        "accession",
        "id",
        "gene_primary",
        "cc_subcellular_location",
        "go_c",
        "cc_function"
    ]

    results = []

    for gene in gene_ids:
        # si GeneID = symbole de gène humain (ex: SLC2A1), on utilise gene:
        queries = [
            f"gene:{gene} AND organism_id:9606 AND reviewed:true",
            f"gene:{gene} AND organism_id:9606",
            f"{gene} AND organism_id:9606"
        ]

        entry = None
        data_to_save = None

        for q in queries:
            params = {
                "query": q,
                "fields": ",".join(fields),
                "format": "json",
                "size": 1
            }

            try:
                r = requests.get(base_url, params=params, timeout=15)
            except requests.exceptions.RequestException as e:
                print(f"Erreur réseau pour {gene} : {e}")
                break

            if r.status_code == 200:
                data = r.json()
                hits = data.get("results", [])
                if hits:
                    entry = hits[0]
                    data_to_save = data
                    break
            else:
                print(f"Erreur UniProt {r.status_code} pour {gene}")

            time.sleep(0.3)

        if save_json and data_to_save:
            with open(json_dir / f"{gene}.json", "w", encoding="utf-8") as f:
                json.dump(data_to_save, f, indent=2, ensure_ascii=False)

        # ------------- Extraction des informations -----------------
        uniprot_id = ""
        accession = ""
        confirmed_loc_list = []      # localisation « consensus » UniProt
        er_subtype_list = []         # sous-types ER (optionnel, filtré)
        go_locterms = []             # GO cellular component (putative)
        function_texts = []          # fonction

        if entry is not None:
            accession = entry.get("primaryAccession", "")
            uniprot_id = entry.get("uniProtkbId", "")

            # Commentaires UniProt (fonction, localisation)
            for comment in entry.get("comments", []):
                ctype = comment.get("commentType")

                if ctype == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []):
                        loc_txt = loc.get("location", {}).get("value", "")
                        if loc_txt:
                            confirmed_loc_list.append(loc_txt)
                            # exemple de détection plus fine pour ER
                            if "endoplasmic reticulum" in loc_txt.lower():
                                er_subtype_list.append(loc_txt)

                if ctype == "FUNCTION":
                    for txt in comment.get("texts", []):
                        val = txt.get("value", "")
                        if val:
                            function_texts.append(val)

            # GO cellular component pour localisation putative
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "GO":
                    # GO term de type C (cellular component)
                    prop_type = {p.get("value", "") for p in xref.get("properties", [])
                                 if p.get("key") == "aspect"}
                    if "C" in prop_type:
                        go_id = xref.get("id", "")
                        term_name = ""
                        for p in xref.get("properties", []):
                            if p.get("key") == "term":
                                term_name = p.get("value", "")
                        if term_name:
                            go_locterms.append(f"{go_id} ({term_name})")

        # Assemblage des champs texte
        confirmed_loc = "; ".join(sorted(set(confirmed_loc_list))) if confirmed_loc_list else ""
        er_subtype = "; ".join(sorted(set(er_subtype_list))) if er_subtype_list else ""
        putative_loc = "; ".join(sorted(set(go_locterms))) if go_locterms else ""
        function_str = " ".join(function_texts) if function_texts else ""

        results.append({
            "GeneID": gene,
            "Uniprot ID": uniprot_id,
            "Confirmed localization (UniProt consensus)": confirmed_loc,
            "ER location subtype": er_subtype,
            "Putative localization GO annotation": putative_loc,
            "Function": function_str,
            "UniProt accession": accession
        })

        time.sleep(0.3)

    return pd.DataFrame(results)

# =====================================================
# 4. Pipeline principal
# =====================================================
df_loc = query_uniprot_localisation(gene_ids, save_json=True)

# On peut soit :
#   - Sauvegarder uniquement le tableau localisation
#   - Ou le fusionner avec le fichier d'origine
# Ici, on ne garde que les colonnes demandées dans la question

cols_finales = [
    "GeneID",
    "Uniprot ID",
    "Confirmed localization (UniProt consensus)",
    "ER location subtype",
    "Putative localization GO annotation",
    "Function"
]

df_output = df_loc[cols_finales]

# =====================================================
# 5. Export Excel
# =====================================================
df_output.to_excel(output_file, index=False, engine="openpyxl")

print(f"✅ Fichier exporté : {output_file}")
print(f"📁 JSON UniProt sauvegardés dans : {json_dir}")