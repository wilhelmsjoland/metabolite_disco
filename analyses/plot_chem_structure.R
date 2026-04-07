library(reticulate)
library(readxl)

draw_out_path <- paste0(
  "/Users/wilhelm/Library/CloudStorage",
  "/ProtonDrive-wilhelm.sjoland@proton.me-folder",
  "/01_juniper/01_arbete/01_projekt/03_psm",
  "/molekyler/vektorbilder"
)

glycone_list <- readxl::read_xlsx(
  path = paste0(
    "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/",
    "molekyler/information/glycone_list.xlsx"
  ),
  sheet = "molecules"
)

reticulate::use_condaenv("chem_env", required = TRUE)
Chem <- reticulate::import("rdkit.Chem")
Draw <- reticulate::import("rdkit.Chem.Draw")
mol_draw <- reticulate::import("rdkit.Chem.Draw.rdMolDraw2D")
inchi_mod <- reticulate::import("rdkit.Chem.inchi")

reticulate::py_run_string("
import re

def draw_from_inchi(inchi_str, name, output_path):
    from rdkit.Chem.inchi import MolFromInchi
    from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG

    mol = MolFromInchi(inchi_str)
    drawer = MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 5.5
    opts.fixedFontSize = 35

    green = (0.001, 0.557, 0.000)

    # Set ALL atoms to green, no highlights
    opts.updateAtomPalette({i: green for i in range(119)})

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_text = drawer.GetDrawingText()

    # Find hydroxyl oxygens and their bonds
    oh_atom_ids = set()
    oh_bond_ids = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0:
            oh_atom_ids.add(atom.GetIdx())
            for bond in atom.GetBonds():
                oh_bond_ids.add(bond.GetIdx())

    red = '#f00'

    # Recolor OH bond paths in SVG
    for bid in oh_bond_ids:
        svg_text = re.sub(
            r\"(<path\\s+class='bond-\"
            + str(bid)
            + r\"\\s[^']*'[^>]*?)stroke:#[0-9A-Fa-f]{6}\",
            r'\\1stroke:' + red,
            svg_text,
        )

    # Recolor OH atom label paths in SVG
    for aid in oh_atom_ids:
        svg_text = re.sub(
            r\"(<path\\s+class='atom-\"
            + str(aid)
            + r\"'.*?)fill='#[0-9A-Fa-f]{6}'\",
            r\"\\1fill='\" + red + r\"'\",
            svg_text,
            flags=re.DOTALL,
        )

    with open(output_path, 'w') as f:
        f.write(svg_text)
")

draw_from_inchi <- function(inchi, name, output_path) {
  reticulate::py$draw_from_inchi(inchi, name, output_path)
}

for (i in seq_len(nrow(glycone_list))) {
  draw_from_inchi(
    inchi = glycone_list$InChI[i],
    name = glycone_list$molecule[i],
    output_path = file.path(
      draw_out_path,
      paste0(glycone_list$molecule[i], ".svg")
    )
  )
}
