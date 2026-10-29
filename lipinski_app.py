import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

st.title("Lipinski Rule of Five Checker")

smiles = st.text_input("Type a SMILES:", "CC(=O)OC1=CC=CC=C1C(=O)O")
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    st.error("That SMILES is not valid. Check it and try again.")
else:
    weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)

    checks = {
        "Weight under 500": weight < 500,
        "LogP under 5": logp < 5,
        "H-bond donors under 5": donors < 5,
        "H-bond acceptors under 10": acceptors < 10,
    }

    st.write("Weight:", round(weight, 2))
    st.write("LogP:", round(logp, 2))
    st.write("Donors:", donors)
    st.write("Acceptors:", acceptors)

    passed = 0
    for name, ok in checks.items():
        if ok:
            st.success(name + ": PASS")
            passed = passed + 1
        else:
            st.error(name + ": FAIL")

    if passed >= 3:
        st.balloons()
        st.write("Drug-like: this molecule follows the Rule of Five.")
    else:
        st.write("Not drug-like: it breaks too many rules.")
