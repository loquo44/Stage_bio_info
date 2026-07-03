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

# 1. CONFIGURATION DES CHEMINS DE SORTIE

# Dossier et fichiers où seront sauvegardés les résultats
out_txt = Path(r"C:\Users\meren\Desktop\stage_info\string_results.txt")
json_dir = Path(r"C:\Users\meren\Desktop\stage_info\string_json_output")

# Crée le dossier pour les JSON s'il n'existe pas
json_dir.mkdir(exist_ok=True, parents=True)

# Crée le dossier parent du fichier texte s'il n'existe pas
out_txt.parent.mkdir(exist_ok=True, parents=True)

# ENTRÉES UTILISATEUR
# Demande le nom de la protéine à chercher
prot_name = input("Nom de la protéine : ").strip()

# Bloque le script si l'utilisateur ne tape rien
if not prot_name:
    raise ValueError("Le nom de la protéine ne peut pas être vide.")

# Demande le code de l'espèce, avec l'humain (9606) par défaut
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
    # Interroge l'API UniProt pour avoir le nom scientifique de l'espèce
    res_org = requests.get(f"https://rest.uniprot.org/taxonomy/{tax_id}", timeout=10)
    if res_org.status_code == 200:
        # Extrait le nom scientifique du JSON reçu
        org_name = res_org.json().get("scientificName", f"ID {tax_id}")
    else:
        org_name = f"ID {tax_id}"
except Exception:
    # Solution de secours si la connexion internet échoue
    org_name = f"ID {tax_id}"
print(f"Organisme retenu : {org_name}\n")

# 4. REQUÊTE ET EXTRACTION DEPUIS L'API STRING

# URL de l'API STRING pour récupérer les annotations fonctionnelles
url_string = "https://string-db.org/api/json/functional_terms"
fields = [
        "category", # c'est la categorie à laquelle appartient la protéine
        "term", # c'est l'identifiant unique en format str
        "description", # recupère la derscription de la prot STRING en format str
        "StringIds", # récupère lid STRING en format str
    ] # liste des champs que l'on souhaite récupérer depuis STRING

# Paramètres à envoyer à l'API STRING
api_params = {
    "identifiers": prot_name,
    "species": tax_id,
    "fields": ",".join(fields),
    "format": "json"
}

# Liste pour stocker les lignes de résultats formatées
results = []
raw_json = None

print(f" Connexion à STRING-DB pour : {prot_name}")
try:
    # Envoi de la requête GET à STRING
    response = requests.get(url_string, params=api_params, timeout=15)
    if response.status_code == 200:
        # Récupère le contenu JSON si la requête a réussi
        raw_json = response.json()
    else:
        print(f" Erreur STRING {response.status_code} pour {prot_name}")
except requests.exceptions.RequestException as error:
    print(f" Erreur réseau lors de l'appel à STRING : {error}")

# Sauvegarde du fichier JSON brut s'il contient des données
if raw_json:
    # Nettoie le nom de la protéine pour éviter les caractères interdits dans les fichiers
    safe_name = "".join([c if c.isalnum() else "_" for c in prot_name])
    with open(json_dir / f"{safe_name}_string.json", "w", encoding="utf-8") as f_json:
        # Écrit le JSON de manière lisible (indentation de 2 espaces)
        json.dump(raw_json, f_json, indent=2, ensure_ascii=False)

# Traitement des données extraites du JSON
if raw_json:
    # Parcourt chaque annotation renvoyée par STRING
    for item in raw_json:
        category = item.get("category", "")
        term_id = item.get("term", "")
        desc = item.get("description", "")
        
        # Initialise les variables de texte pour cette annotation
        loc_conf = ""
        loc_er = ""
        loc_go = ""
        func_txt = ""
        
        # Filtre les annotations liées aux composants cellulaires (Localisation)
        if category == "Component" or "GO:" in term_id:
            loc_go = f"{term_id} ({desc})"
            
            # Recherche des mots-clés pour valider une localisation générale
            if any(key in desc.lower() for key in ["membrane", "nucleus", "cytoplasm", "mitochondrion"]):
                loc_conf = desc
            
            # Repère si l'annotation concerne le Réticulum Endoplasmique (ER)
            if "endoplasmic reticulum" in desc.lower():
                loc_er = desc
                
        # Filtre les annotations liées aux processus biologiques (Fonctions)
        if category == "Process":
            func_txt = desc

        # Crée un bloc de texte propre pour cette annotation spécifique
        block = (
            f"Gene ID: {prot_name}\n"
            f"Complex/Term ID: {term_id}\n"
            f"Confirmed localization: {loc_conf}\n"
            f"ER location subtype: {loc_er}\n"
            f"Putative GO annotation: {loc_go}\n"
            f"Function: {func_txt}\n"
        )
        # Ajoute ce bloc à notre liste de résultats
        results.append(block)

# Si l'API n'a retourné aucun résultat pour cette protéine
if not results:
    results.append(f"Aucune annotation trouvée dans STRING pour {prot_name} (TaxID: {tax_id}).")

# 5. ÉCRITURE DU FICHIER TEXTE FINAL

# Ouvre le fichier texte en mode écriture 
with open(out_txt, "w", encoding="utf-8") as f_txt:
    # Écrit en-tête général au début du fichier
    f_txt.write(f"RÉSULTATS DE LA RECHERCHE STRING-DB\n")
    f_txt.write(f"Protéine cible : {prot_name}\n")
    f_txt.write(f"Organisme : {org_name} (TaxID: {tax_id})\n")
    
    # Écrit chaque bloc d'annotation l'un après l'autre
    for result_block in results:
        f_txt.write(result_block + "\n")

# Messages de confirmation de fin de script
print("\n--- TRAITEMENT TERMINÉ ---")
print(f" Fichier texte généré : {out_txt}")
print(f" Fichier JSON brut sauvegardé dans : {json_dir}")
print(f"Statistiques globales enregistrées pour l'organisme : {org_name}")