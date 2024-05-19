import json

def format_mol_for_json(mol_file_path):
    with open(mol_file_path, 'r') as file:
        mol_data = file.read()

    # Escape new lines for JSON
    mol_data_json_friendly = mol_data.replace('\n', '\\n')
    return mol_data_json_friendly

# Assume you have the path to your MOL files
mol1 = format_mol_for_json('benzene-3D-structure-CT1001419667.sdf')
mol2 = format_mol_for_json('ISOAMYL-ACETATE-3D-structure-CT1002615369.sdf')


# Prepare JSON payload
json_payload = json.dumps({
    "mol1": mol1,
    "mol2": mol2
})

print(json_payload)
