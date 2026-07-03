# -*- coding: utf-8 -*-

"""
    Created on Thu Jul 7 2026
    
    @author: meren
"""

import pandas as pd
import os
import requests
import json
from pathlib import Path

output_file = Path(r"C:\Users\meren\Desktop\stage_info\complexportal_localisation_from_complexID.xlsx")
json_dir=Path(r"C:\Users\meren\Desktop\stage_info\complexportal_json_localisation")

json_dir.mkdir(exist_ok=True) # crée le dossier json_dir s'il n'existe pas déjà

Complexe_id = "CPX-5977" # valeur à donner par l'utilisateur, peut être modifiée pour tester d'autres complexes
organisme = 9606 # valeur à donner par l'utilisateur, peut être modifiée pour tester d'autres organismes

print(f"recherche du nom de l'organisme : {organisme}")
try: 
    org_res = requests.get(f"https://rest.uniprot.org/taxonomy/{organisme}")
    data_org = org_res.json()
    nom_organisme = data_org.get("scientificName", f"ID {organisme}")
except Exception:
    nom_organisme = f"ID {organisme}" #solution de secour si le réseau est problématique
print(f"Nom de l'organisme : {nom_organisme}")

print(f"Recherche du complexe : {Complexe_id}")
try:
    comp_res = requests.get(f"https://www.ebi.ac.uk/ebisearch/search?db=allebi&sortignorenull=true&query={Complexe_id}&size=15&requestFrom=masthead-black-bar")
    data_comp = comp_res.json()
except Exception:
    print(f"Erreur lors de la récupération du complexe : {Complexe_id}")

