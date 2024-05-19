from flask import Flask, request, send_file, jsonify
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw
import traceback
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import DataStructs

app = Flask(__name__)

@app.route('/', methods = ['get'])
def landingpage():
    return ("hello world")

@app.route('/convert', methods=['POST'])
def convert():
    data = request.get_json()
    mol_data = data['mol'].replace('\\n', '\n')
    try:
        mol = Chem.MolFromMolBlock(mol_data)
        if mol:
            img = Draw.MolToImage(mol)
            img_byte_arr = BytesIO()
            img.save(img_byte_arr, format='PNG')
            img_byte_arr.seek(0)
            return send_file(img_byte_arr, mimetype='image/png')
        else:
            return jsonify({'error': 'Could not parse MOL data'}), 400
    except Exception as e:
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500
    

@app.route('/compare', methods=['POST'])
def compare_molecules():
    data = request.get_json(force = True)
    mol1_data = data['mol1'].replace('\\n', '\n')  # Decode JSON newline escapes
    mol2_data = data['mol2'].replace('\\n', '\n')
    print(mol1_data, mol2_data)

    try:
        mol1 = Chem.MolFromMolBlock(mol1_data)
        mol2 = Chem.MolFromMolBlock(mol2_data)


        if not mol1 or not mol2:
            return jsonify({'error': 'Invalid MOL data provided for one or both molecules.'}), 400

        # Calculate similarity
        fp1 = FingerprintMols.FingerprintMol(mol1)
        fp2 = FingerprintMols.FingerprintMol(mol2)
        similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
        
        return jsonify({"similarity": f"{similarity * 100:.2f}%"}), 200

    except Exception as e:
        return jsonify({'error': str(e), 'trace': traceback.format_exc()}), 500


if __name__ == '__main__':
    app.run(debug=True, port=5000)