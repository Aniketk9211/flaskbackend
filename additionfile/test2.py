# from rdkit import Chem
# from rdkit.Chem import Draw
# import matplotlib.pyplot as plt
# from io import BytesIO

# # stringWithMolData=open('benzene-3D-structure-CT1001419667.sdf','r').read()
# stringWithMolData= '''CT1001419667


#  12 12  0  1  0               999 V2000
#    -0.0167    1.3781    0.0096 C   0 00  0  0  0  0  0  0  0  0  0  0
#     0.0021   -0.0041    0.0020 C   0 00  0  0  0  0  0  0  0  0  0  0
#     1.1709    2.0855    0.0021 C   0 00  0  0  0  0  0  0  0  0  0  0
#     1.2084   -0.6789   -0.0132 C   0 00  0  0  0  0  0  0  0  0  0  0
#     2.3960    0.0285   -0.0212 C   0 00  0  0  0  0  0  0  0  0  0  0
#     2.3773    1.4107   -0.0131 C   0 00  0  0  0  0  0  0  0  0  0  0
#    -0.9592    1.9054    0.0170 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.9258   -0.5567    0.0083 H   0  0  0  0  0  0  0  0  0  0  0  0
#     1.1563    3.1654    0.0077 H   0  0  0  0  0  0  0  0  0  0  0  0
#     1.2231   -1.7588   -0.0184 H   0  0  0  0  0  0  0  0  0  0  0  0
#     3.3385   -0.4987   -0.0324 H   0  0  0  0  0  0  0  0  0  0  0  0
#     3.3051    1.9634   -0.0197 H   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2 02  0  1  0  0
#   1  3 01  0  1  0  0
#   1  7  1  0  0  0  0
#   2  4 01  0  1  0  0
#   2  8  1  0  0  0  0
#   3  6 02  0  1  0  0
#   3  9  1  0  0  0  0
#   4  5 02  0  1  0  0
#   4 10  1  0  0  0  0
#   5  6 01  0  1  0  0
#   5 11  1  0  0  0  0
#   6 12  1  0  0  0  0
# M  END
# $$$$'''
# print(type(stringWithMolData))
# m = Chem.MolFromMolBlock(stringWithMolData)
# img = Draw.MolToImage(m)
# img_byte_arr = BytesIO()
# img.save(img_byte_arr, format='PNG')
# img_byte_arr.seek(0)

# plt.imshow(img_byte_arr)
# plt.axis('off')  # Hide axes
# plt.show()



from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

# Define two SMILES strings (representing chemical structures)
smiles1 = 'CCOCN(C)(C)'
smiles2 = 'CCO'

# Convert SMILES strings to RDKit Mol objects
mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)

# Generate fingerprints for each molecule
fp1 = FingerprintMols.FingerprintMol(mol1)
fp2 = FingerprintMols.FingerprintMol(mol2)

# Calculate Tanimoto similarity
similarity = DataStructs.FingerprintSimilarity(fp1, fp2)

# Convert similarity to percentage
similarity_percentage = similarity * 100

print(f"Similarity percentage: {similarity_percentage:.2f}%")


from urllib.parse import unquote

@app.route('/convert', methods=['GET'])
def convert():
    mol_data = request.args.get('moldata', '')
    decoded_mol_data = unquote(mol_data)  # Decode URL-encoded data
    try:
        mol = Chem.MolFromMolBlock(decoded_mol_data)
        if mol:
            img = Draw.MolToImage(mol)
            img_byte_arr = BytesIO()
            img.save(img_byte_arr, format='PNG')
            img_byte_arr.seek(0)
            return send_file(img_byte_arr, mimetype='image/png', as_attachment=False)
        else:
            return jsonify({'error': 'Could not parse MOL data'}), 400
    except Exception as e:
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500
