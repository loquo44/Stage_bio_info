# -*- coding: utf-8 -*-
"""
    Created on Tue Jul 2 2026
    @author: meren
"""

import os
import requests
import json
from pathlib import Path
import time
import pandas as pd  # Pour lire xlsx

# 1. CONFIGURATION DES CHEMINS DE SORTIE

# Dossier ou seront sauvegardes les rapports Markdown 
out_dir = Path(r"C:\Users\meren\Desktop\stage_info\string_reports")
out_dir.mkdir(exist_ok=True, parents=True)

# Fichier xlsx source UniProt
EXCEL_FILE = Path(r"base_de_donnees_STRING/test_ms_uniport.xlsx")
UNIPROT_COLUMN = "PG,Genes"  # Nom de la colonne cible

# DICTIONNAIRE GLOBAL 

# Ce dictionnaire stockera : { "Nom_Prot_Origine": [liste_preferredNames_valides] }
dictionnaire_enrichissement = {}

# ENTREES UTILISATEUR

tax_input = input("TaxID de l'organisme, par défaut : 9606 (Humain) : ").strip()
if tax_input == "":
    tax_id = 9606
else:
    try:
        tax_id = int(tax_input)
    except ValueError:
        print(" ID invalide. Utilisation du TaxID par défaut : 9606")
        tax_id = 9606

# 3. IDENTIFICATION DE L'ORGANISME VIA UNIPROT

print(f"\nRecherche du nom de l'organisme pour le TaxID {tax_id}")
try:
    res_org = requests.get(f"https://rest.uniprot.org/taxonomy/{tax_id}", timeout=10)
    if res_org.status_code == 200:
        org_name = res_org.json().get("scientificName", f"ID {tax_id}")
    else:
        org_name = f"ID {tax_id}"
except Exception:
    org_name = f"ID {tax_id}"
print(f"Organisme retenu : {org_name}\n")

# LECTURE DU FICHIER EXCEL (Saute la 1ere ligne et cible '(nom_colonne)')

print(f"Lecture du fichier Excel : {EXCEL_FILE}")

# skiprows=1 permet de sauter la premiere ligne parasite
df = pd.read_excel(EXCEL_FILE, skiprows=1)
df = df[UNIPROT_COLUMN].dropna().astype(str).unique().tolist()

# URL et champs pour l'API STRING
url_string = "https://string-db.org/api/json/functional_terms"
fields = [
        "category", # c'est la categorie à laquelle appartient la proteine
        "term", # c'est l'identifiant unique en format str
        "description", # recupere la derscription de la prot STRING en format str
        "stringIds", # recupere l'id STRING en format str
        "proteinCount", # recupere le nombre de proteines dans la categorie
        "preferredNames", # recupere le nom preferentiel de la proteine STRING
    ] # liste des champs que l'on souhaite recuperer depuis STRING

# Boucle ligne par ligne sur le fichier Excel
for row in df:
    prot_name = row.strip()
    
    # Configuration des paramètres STRING pour cette protéine
    api_params = {
    "preferredNames": prot_name,
    "species": tax_id,
    "fields": ",".join(fields),
    "format": "json",
    "size": "1"
    }

    raw_json = None
    print(f" Connexion à STRING-DB pour : {prot_name}")
    try:
        response = requests.get(url_string, params=api_params, timeout=10)
        if response.status_code == 200:
            raw_json = response.json()
        else:
            print(f" Erreur STRING {response.status_code} pour {prot_name}")
    except requests.exceptions.RequestException as error:
        print(f" Erreur réseau lors de l'appel à STRING : {error}")

    # On prepare une liste vide pour stocker les preferredNames VALIDES de cette proteine precise
    liste_preferred_valides = []

    results = []
    one_block = []
    all_block = []

    # TRAITEMENT ET FILTRAGE
    if raw_json:
        print("ok")
        # filtre demandes
        filtre_input = input(f"[{prot_name}] Entrez le filtre pour les categories, par defaut(None) : ").strip()
        filtre = None if filtre_input == "" else filtre_input
        
        filtre_hypothetic = input("Voulez-vous filtrer les annotations hypothétiques ? (oui/non, par défaut : non) : ").strip().lower()
        filtre_h = "non" if filtre_hypothetic == "" else filtre_hypothetic
        
        filtre_unknow = input("Voulez-vous filtrer les annotations inconnues ? (oui/non, par défaut : non) : ").strip().lower()
        filtre_u = "non" if filtre_unknow == "" else filtre_unknow  
        
        for item in raw_json:
            category = item.get("category", "")
            term_id = item.get("term", "")
            desc = item.get("description", "")
            string_ids = item.get("stringIds", [])
            protein_count = item.get("proteinCount", 0)
            preferredNames = item.get("preferredNames", [])
            
            # --- INTERCEPTION POUR LE DICTIONNAIRE ---
            # On verifie les filtres textuels avant de valider les preferredNames pour le dictionnaire
            is_valid = True
            if filtre_h == "oui" and "hypothetical" in desc.lower():
                is_valid = False
            if filtre_u == "oui" and "unknown" in desc.lower():
                is_valid = False
                
            if is_valid:
                # Si l'annotation passe les criteres, on garde ses preferredNames pour le dictionnaire
                liste_preferred_valides.extend(preferredNames)
            # -----------------------------------------

            term = ""
            loc_go = f"{desc}"
            mot = ""
            
            if filtre == None :
                term = term_id
                block = (
                    f"### Category: {category}\n"
                    f"### Gene ID: {prot_name}\n"
                    f"### Complex/Term ID: {term}\n"
                    f"### Putative GO annotation: {loc_go}\n"
                    f"### Protein Count: {protein_count}\n"
                    "\n"+"| STRING IDs           | Preferred Names \n"
                )
                for i in range(len(string_ids)):
                    block += (f"| {string_ids[i]} | {preferredNames[i]} \n")
                all_block.append(block)
                
            else:
                for id in range(len(term_id)):
                    mot += term_id[id]
                    
                if filtre == mot:
                    term = mot
                    block = (
                        f"### Category: {category}\n"
                        f"### Gene ID: {prot_name}\n"
                        f"### Complex/Term ID: {term}\n"
                        f"### Putative GO annotation: {loc_go}\n"
                        f"### Protein Count: {protein_count}\n"
                        "\n"+"| STRING IDs           | Preferred Names \n"
                    )
                    for i in range(len(string_ids)):
                        block += (f"| {string_ids[i]} | {preferredNames[i]} \n")
                    one_block.append(block)
        time.sleep(0.2) # pause pour l'API (evite DDOS)
        results = one_block if one_block else all_block

if not results:
    results.append(f"Aucune annotation trouvée dans STRING pour {prot_name} (TaxID: {tax_id}).")

# --- REMPLISSAGE DU DICTIONNAIRE GLOBAL ---
# On applique un set() pour eviter d'avoir plusieurs fois la même proteine si elle apparaît dans plusieurs categories
if liste_preferred_valides:
    dictionnaire_enrichissement[prot_name] = list(set(liste_preferred_valides))

# 5. ECRITURE DU FICHIER TEXTE FINAL
safe_name = "".join([c if c.isalnum() else "_" for c in prot_name])
out_txt = out_dir / f"{safe_name}"

with open(out_txt, "w", encoding="utf-8") as f_txt:
    f_txt.write(f"# RÉSULTATS DE LA RECHERCHE STRING-DB\n")
    f_txt.write("\n" + f"## Protéine cible : {prot_name}\n")
    f_txt.write(f"## Organisme : {org_name} (TaxID: {tax_id})\n")
    
    for result_block in results:
        f_txt.write("\n" + result_block + "\n")
        
print(f" Fichier généré : {out_txt}\n")

# 6. SAUVEGARDE DU DICTIONNAIRE 

# On sauvegarde le dictionnaire dans un fichier JSON pour pouvoir le recuperer
dict_output_path = out_dir.parent / "dictionnaire_enrichissement.json"
with open(dict_output_path, "w", encoding="utf-8") as f_dict:
    json.dump(dictionnaire_enrichissement, f_dict, indent=4, ensure_ascii=False)

print("\n--- TRAITEMENT TERMINÉ ---")
print(f" le rapport de toutes les proteines sont dans : {out_dir}")
print(f" Le dictionnaire global a été sauvegardé ici : {dict_output_path}")