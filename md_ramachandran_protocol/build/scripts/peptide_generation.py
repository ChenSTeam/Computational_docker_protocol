from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw, AllChem

#useful smiles
ACETYL_SMILES = "CC(=O)O"
AMINE_SMILES = "NC"
Alanine_SMILES = "N[C@@H](C)C(=O)O"
Arginine_SMILES = "N[C@@H](CCCNC(=N)N)C(=O)O"
Asparagine_SMILES = "N[C@@H](CC(=O)N)C(=O)O"
Aspartate_SMILES = "N[C@@H](CC(=O)O)C(=O)O"
Cysteine_SMILES = "N[C@@H](CS)C(=O)O"
Glutamine_SMILES = "N[C@@H](CCC(=O)N)C(=O)O"
Glutamate_SMILES = "N[C@@H](CCC(=O)O)C(=O)O"
Glycine_SMILES = "NCC(=O)O"
Histidine_SMILES = "N[C@@H](Cc1c[nH]cn1)C(=O)O"
Isoleucine_SMILES = "N[C@@H]([C@@H](C)CC)C(=O)O"
Leucine_SMILES = "N[C@@H](CC(C)C)C(=O)O"
Lysine_SMILES = "N[C@@H](CCCCN)C(=O)O"
Methionine_SMILES = "N[C@@H](CCSC)C(=O)O"
Phenylalanine_SMILES = "N[C@@H](Cc1ccccc1)C(=O)O"
Proline_SMILES = "O=C(O)[C@@H]1CCCN1"
Serine_SMILES = "N[C@@H](CO)C(=O)O"
Threonine_SMILES = "N[C@@H]([C@H](O)C)C(=O)O"
Tryptophan_SMILES = "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"
Tyrosine_SMILES = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"
Valine_SMILES = "N[C@@H](C(C)C)C(=O)O"

aa1tosmiles = {
    "A":Alanine_SMILES,
    "R":Arginine_SMILES,
    "N":Asparagine_SMILES,
    "D":Aspartate_SMILES,
    "C":Cysteine_SMILES,
    "Q":Glutamine_SMILES,
    "E":Glutamate_SMILES,
    "G":Glycine_SMILES,
    "H":Histidine_SMILES,
    "I":Isoleucine_SMILES,
    "L":Leucine_SMILES,
    "K":Lysine_SMILES,
    "M":Methionine_SMILES,
    "F":Phenylalanine_SMILES,
    "P":Proline_SMILES,
    "S":Serine_SMILES,
    "T":Threonine_SMILES,
    "W":Tryptophan_SMILES,
    "Y":Tyrosine_SMILES,
    "V":Valine_SMILES
}

#pattern of amide bond reaction on alpha COO and N
#amide_rxn = rdChemReactions.ReactionFromSmarts("[N:1][C:2][C:3](=[O:4])[O:5].[N;H2:6][C:7][C:8](=[O:9])>>[N:1][C:2][C:3](=[O:4])[N;H:6][C:7][C:8](=[O:9])")
amide_rxn = rdChemReactions.ReactionFromSmarts("[N:1][C:2][C:3](=[O:4])[O:5].[N;H1,H2:6][C:7][C:8](=[O:9])>>[N:1][C:2][C:3](=[O:4])[N:6][C:7][C:8](=[O:9])")
acylation_rxn = rdChemReactions.ReactionFromSmarts("[C:3](=[O:4])[O:5].[N;H2:6][C:7][C:8](=[O:9])>>[C:3](=[O:4])[N;H:6][C:7][C:8](=[O:9])")
amidation_rxn = rdChemReactions.ReactionFromSmarts("[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7]>>[N:1][C:2][C:3](=[O:4])[N:6][C:7]")

def is_alpha_amino_acid(mol):
    if mol is None:
        return False
    # Any N that can act as amine (non-quaternary)
    patt_N = Chem.MolFromSmarts("[N;X3;!$([N+])]")
    has_N = mol.HasSubstructMatch(patt_N)
    # Carboxylate (either protonated or deprotonated)
    patt_COO = Chem.MolFromSmarts("C(=O)[O,O-]")
    has_COO = mol.HasSubstructMatch(patt_COO)
    if not (has_N and has_COO):
        return False
    # Check alpha carbon connectivity (C bonded to N and C(=O))
    patt_alpha = Chem.MolFromSmarts("[C;X4]([N])C(=O)[O,O-]")
    has_alpha = mol.HasSubstructMatch(patt_alpha)
    return has_alpha

def connect_residues(mol1, mol2):
    prods = amide_rxn.RunReactants((mol1, mol2))
    if not prods:
        raise ValueError("Failed to connect residues")
    prod = prods[0][0]
    Chem.SanitizeMol(prod)
    return prod

def peptide_formation(mol_list):
    prod = mol_list[0]
    for mol in mol_list[1:]:
        prod = connect_residues(prod, mol)
    Chem.SanitizeMol(prod)
    return prod

def peptide_protect(mol):
    acetyl = Chem.MolFromSmiles(ACETYL_SMILES)
    amine = Chem.MolFromSmiles(AMINE_SMILES)
    prod = acylation_rxn.RunReactants((acetyl, mol))[0][0]
    prod = amidation_rxn.RunReactants((prod, amine))[0][0]
    Chem.SanitizeMol(prod)
    return prod

def read_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)
    return mol

def natural_peptide_generate(seq, protect=True):
    seq = seq.upper()
    mol_list = []
    for res in seq:
        mol_list.append(read_mol(aa1tosmiles[res]))
    pept = peptide_formation(mol_list)
    if protect:
        pept = peptide_protect(pept)
    return pept

def uaa_tripeptide_generate(uaa_smiles, flanking_res='A', protect=True):
    uaa = read_mol(uaa_smiles)
    if not is_alpha_amino_acid(uaa):
        raise ValueError('Not alpha amino acid!')

    mol_list = [read_mol(aa1tosmiles[flanking_res]), uaa, read_mol(aa1tosmiles[flanking_res])]
    pept = peptide_formation(mol_list)
    if protect:
        pept = peptide_protect(pept)
    return pept

def generate_3d_struct(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    return mol
