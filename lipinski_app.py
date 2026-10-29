import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw

st.title("Lipinski Rule of Five Checker")

smiles = st.text_input("Type a SMILES:", "CC(=O)OC1=CC=CC=C1C(=O)O")
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    st.error("That SMILES is not valid. Check it and try again.")
else:
    # Render the molecular structure image
    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, caption="Molecular Structure")

    weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    donors = Lipinski.NumHDonors(mol)
    acceptors = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)

    checks = {
        "Weight under 500": weight < 500,
        "LogP under 5": logp < 5,
        "H-bond donors under 5": donors < 5,
        "H-bond acceptors under 10": acceptors < 10,
    }

    st.write("**Computed Properties:**")
    st.write(f"- **Weight:** {round(weight, 2)} g/mol")
    st.write(f"- **LogP:** {round(logp, 2)}")
    st.write(f"- **H-Bond Donors:** {donors}")
    st.write(f"- **H-Bond Acceptors:** {acceptors}")
    st.write(f"- **Rotatable Bonds:** {rot_bonds}")

    passed = 0
    st.write("---")
    st.write("**Rule Compliance:**")
    for name, ok in checks.items():
        if ok:
            st.success(name + ": PASS")
            passed += 1
        else:
            st.error(name + ": FAIL")

    if passed >= 3:
        st.balloons()
        st.success("Drug-like: this molecule follows the Rule of Five.")
    else:
        st.warning("Not drug-like: it breaks too many rules.")
