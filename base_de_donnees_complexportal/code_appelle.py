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

json_dir=Path(r"C:\Users\meren\Desktop\stage_info\complexportal_json_localisation")
json_dir.mkdir(exist_ok=True) # crée le dossier json_dir s'il n'existe pas déjà

Complexe = "CPX-5977"
organisme = 9606

print(f"recherche du nom de l'organisme : {organisme}")
try: 
    org_res = requests.get(f"https://rest.uniprot.org/taxonomy/{organisme}")
    data_org = org_res.json()
    nom_organisme = data_org.get("scientificName", f"ID {organisme}")
except Exception:
    nom_organisme = f"ID {organisme}" #solution de secour si le réseau est problématique
print(f"Nom de l'organisme : {nom_organisme}")

print(f"Recherche du complexe : {Complexe}")
try:
    comp_res = requests.get(f"https://www.ebi.ac.uk/complexportal/complexes/{Complexe}")
    data_comp = comp_res.json()
except Exception:
    print(f"Erreur lors de la récupération du complexe : {Complexe}")
    
        
