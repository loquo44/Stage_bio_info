# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 07:51:50 2026

@author: teo.merenda
"""

import pandas as pd
import requests
import time
import json
from pathlib import Path

# =====================================================
# 1. Fichiers
# creer 3 chemin distinct afin de récuperer des fichiers utiles au code
# =====================================================
input_file = Path(r"C:\\Users\\teo.merenda\\Documents\\DTT red+ ox- no mix .xlsx")
output_file = Path(r"C:\\Users\\teo.merenda\\Documents\\uniprot_localisation_from_GeneID.xlsx")
json_dir = Path(r"C:\\Users\\teo.merenda\\Documents\\uniprot_json_localisation")

json_dir.mkdir(exist_ok=True) # crée le dossier json_dir s'il n'existe pas déjà

# =====================================================
# 2. Lecture du fichier (DataFrame)
#    On utilise la colonne 'GeneID' comme identifiant de gène
# =====================================================
df = pd.read_excel(input_file) # permet de lire le fichier excel

if "GeneID" not in df.columns: # si la colonne "GeneID" n'est pas dans les colonnes du DataFrame alors renvoie une erreur
    raise ValueError(f"La colonne 'GeneID' n'existe pas dans le fichier. Colonnes disponibles : {list(df.columns)}")

gene_ids = df["GeneID"].dropna().astype(str).unique().tolist() # génère une liste des toutes les lignes différentes de la colonne demandée et la supprime dès que l'information du Gene ID est manquante
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
    ] # liste des champs que l'on souhaite récupérer depuis UniProt

    results = [] # liste vide qui va contenir les résultats de la requête

    for gene in gene_ids:
        # si GeneID = symbole de gène humain (ex: SLC2A1), on utilise gene:
        queries = [
            f"gene:{gene} AND organism_id:9606 AND reviewed:true",
            f"gene:{gene} AND organism_id:9606",
            f"{gene} AND organism_id:9606"
        ] # liste des requêtes à tester

        entry = None
        data_to_save = None

        for q in queries: # prend le gene donné par queries 
            params = {
                "query": q,
                "fields": ",".join(fields), # appelle field pour la liste des donneés Uniprot
                "format": "json",
                "size": 1
            }

            try :
                r = requests.get(base_url, params=params, timeout=15)
            except requests.exceptions.RequestException as e : 
                print(f"Erreur réseau pour {gene} : {e}")
                break # si il y a une erreur réseau alors affiche l'erreur et sort de la boucle

            if r.status_code == 200:
                data = r.json() # creer un dictionnaire data qui prend en compte ce que renvoie l'appelle HTTP de r en format "json"
                hits = data.get("results", []) # hits va contenir ce que va trouver get dans le dictionnaire data à "results" sinon renvoie une liste vide
                if hits:
                    entry = hits[0] # entry prend la première valeur de hits
                    data_to_save = data # data_to_save prend le dictionnaire data
                    
            else:
                print(f"Erreur UniProt {r.status_code} pour {gene}") # sinon affiche une erreur pour le gene demandé en inquant lequel

            time.sleep(0.3)

        if save_json and data_to_save: # si save_json est vrai et que data_to_save n'est pas vide alors sauvegarde le dictionnaire data_to_save dans un fichier json
            with open(json_dir / f"{gene}.json", "w", encoding="utf-8") as f:
                json.dump(data_to_save, f, indent=2, ensure_ascii=False)

        # ------------- Extraction des informations -----------------
        uniprot_id = ""
        accession = ""
        confirmed_loc_list = []      # localisation « consensus » UniProt
        er_subtype_list = []         # sous-types ER (optionnel, filtré)
        go_locterms = []             # GO cellular component (putative)
        function_texts = []          # fonction

        if entry is not None:   # on vérifie que results est bien dans le dictionnaire data
            accession = entry.get("primaryAccession", "") # cherche la requete dans entry et la met dans accession sinon met une entré vide
            uniprot_id = entry.get("uniProtkbId", "")   # cherche la requete dans entry et la met dans uniprot_id sinon met une entré vide

            # Commentaires UniProt (fonction, localisation)
            for comment in entry.get("comments", []): # prend chaque comment où entry à la valeur demandé 
                ctype = comment.get("commentType")

                if ctype == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []): # prend les localisations pour chaque comment qui a un subcellularLocations sinon prend une liste vide
                        loc_txt = loc.get("location", {}).get("value", "") # donne la valeur de location sinon met une entrée vide
                        if loc_txt: # si loc_txt n'est pas vide alors ajoute loc_txt à la liste confirmed_loc_list
                            confirmed_loc_list.append(loc_txt)
                            # exemple de détection plus fine pour ER
                            if "endoplasmic reticulum" in loc_txt.lower(): #si endoplasmic reticulum est dans loc_txt alors ajoute loc_txt à la liste er_subtype_list
                                er_subtype_list.append(loc_txt)

                if ctype == "FUNCTION": # si ctype est égal à FUNCTION alors prend les textes pour chaque comment qui a un texts sinon prend une liste vide
                    for txt in comment.get("texts", []): # prend les textes pour chaque comment qui a un texts sinon prend une liste vide
                        val = txt.get("value", "") # donne la valeur de txt sinon met une entrée vide
                        if val: # si val n'est pas vide alors ajoute val à la liste function_texts
                            function_texts.append(val)

            # GO cellular component pour localisation putative
            for xref in entry.get("uniProtKBCrossReferences", []): # prend les xref pour chaque fois que la demande a été confirmé dans entry 
                if xref.get("database") == "GO":
                    # GO term de type C (cellular component)
                    # inclu dans prop_type, chaque 
                    prop_type = {p.get("value", "") for p in xref.get("properties", []) if p.get("key") == "aspect"} # donne la valeur de chaque propriété qui a pour clé "aspect" sinon met une entrée vide
                    if "C" in prop_type: # si "C" est dans prop_type alors ajoute à la liste go_locterms le GO ID et le nom du terme
                        go_id = xref.get("id", "")
                        term_name = ""
                        for p in xref.get("properties", []): # prend chaque propriété de xref
                            if p.get("key") == "term": # si la clé de p est "term" alors met la valeur de p dans term_name sinon met une entrée vide
                                term_name = p.get("value", "")
                        if term_name: # si term_name n'est pas vide alors ajoute à la liste go_locterms le GO ID et le nom du terme
                            go_locterms.append(f"{go_id} ({term_name})")
        else:
            print(f"I can't do anything because entry is {entry}")

        # Assemblage des champs texte
        confirmed_loc = "; ".join(sorted(set(confirmed_loc_list))) if confirmed_loc_list else "" # si confirmed_loc_list n'est pas vide alors join les valeurs de la liste triée et sans doublons avec un "; " sinon met une entrée vide
        er_subtype = "; ".join(sorted(set(er_subtype_list))) if er_subtype_list else "" # si er_subtype_list n'est pas vide alors join les valeurs de la liste triée et sans doublons avec un "; " sinon met une entrée vide
        putative_loc = "; ".join(sorted(set(go_locterms))) if go_locterms else "" # si go_locterms n'est pas vide alors join les valeurs de la liste triée et sans doublons avec un "; " sinon met une entrée vide
        function_str = " ".join(function_texts) if function_texts else "" # si function_texts n'est pas vide alors join les valeurs de la liste avec un " " sinon met une entrée vide

        results.append({
            "GeneID": gene,
            "Uniprot ID": uniprot_id,
            "Confirmed localization (UniProt consensus)": confirmed_loc,
            "ER location subtype": er_subtype,
            "Putative localization GO annotation": putative_loc,
            "Function": function_str,
            "UniProt accession": accession
        }) # ajoute à la liste results un dictionnaire avec les valeurs demandées pour chaque gene

        time.sleep(0.3)

    return pd.DataFrame(results) # renvoie un DataFrame avec les résultats de la liste results

# =====================================================
# 4. Pipeline principal
# =====================================================
df_loc = query_uniprot_localisation(gene_ids, save_json=True) # appelle la fonction query_uniprot_localisation avec la liste des gene_ids et save_json à True pour sauvegarder les fichiers json

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
] # liste des colonnes que l'on souhaite garder dans le DataFrame final

df_output = df_loc[cols_finales] # renvoie un DataFrame avec les colonnes demandées dans cols_finales

# =====================================================
# 5. Export Excel
# =====================================================
df_output.to_excel(output_file, index=False, engine="openpyxl") # exporte le DataFrame df_output dans un fichier excel avec le nom de fichier output_file, sans les index et en utilisant le moteur openpyxl

print(f"✅ Fichier exporté : {output_file}")
print(f"📁 JSON UniProt sauvegardés dans : {json_dir}")