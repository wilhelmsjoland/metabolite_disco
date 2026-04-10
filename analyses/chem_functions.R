library(reticulate)
library(dplyr)
library(readxl)

reticulate::use_condaenv("chem_env", required = TRUE)
Chem <- reticulate::import("rdkit.Chem")
Draw <- reticulate::import("rdkit.Chem.Draw")
mol_draw <- reticulate::import("rdkit.Chem.Draw.rdMolDraw2D")
inchi_mod <- reticulate::import("rdkit.Chem.inchi")

reticulate::py_run_string("
import re
import math

from rdkit.Chem import rdDepictor
from rdkit.Chem.inchi import MolFromInchi
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, MeanBondLength

FLAVONOID_CLOCKWISE_ROTATION = -22.85348480697037

def rotate_2d(mol, degrees):
    radians = math.radians(degrees)
    cos_theta = math.cos(radians)
    sin_theta = math.sin(radians)

    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    center_x = sum(xs) / len(xs)
    center_y = sum(ys) / len(ys)

    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x = pos.x - center_x
        y = pos.y - center_y
        pos.x = center_x + x * cos_theta - y * sin_theta
        pos.y = center_y + x * sin_theta + y * cos_theta
        conf.SetAtomPosition(i, pos)

def straighten_template(template, fallback_rotation=None):
    conf = template.GetConformer()

    for bond in template.GetBonds():
        begin_atom = template.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = template.GetAtomWithIdx(bond.GetEndAtomIdx())

        if bond.GetBondTypeAsDouble() != 2:
            continue

        if begin_atom.GetAtomicNum() == 8:
            oxygen_idx = begin_atom.GetIdx()
            carbon_idx = end_atom.GetIdx()
        elif end_atom.GetAtomicNum() == 8:
            oxygen_idx = end_atom.GetIdx()
            carbon_idx = begin_atom.GetIdx()
        else:
            continue

        oxygen_pos = conf.GetAtomPosition(oxygen_idx)
        carbon_pos = conf.GetAtomPosition(carbon_idx)
        angle = math.degrees(
            math.atan2(oxygen_pos.y - carbon_pos.y, oxygen_pos.x - carbon_pos.x)
        )
        rotate_2d(template, 90 - angle)
        return template

    if fallback_rotation is not None:
        rotate_2d(template, fallback_rotation)

    return template

def flipped_scaffold_template(inchi_str, fallback_rotation=None):
    mol = MolFromInchi(inchi_str)
    template = MurckoScaffold.GetScaffoldForMol(mol)
    rdDepictor.Compute2DCoords(template)

    # Flip vertically so the flavonoid cores match the conventional
    # left-fused-ring/right-phenyl orientation in the reference image.
    conf = template.GetConformer()
    for i in range(template.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        pos.y = -pos.y
        conf.SetAtomPosition(i, pos)

    return straighten_template(template, fallback_rotation)

alignment_templates = [
    flipped_scaffold_template(
        'InChI=1S/C15H10O5/c16-9-3-1-8(2-4-9)13-7-12(19)15-11(18)5-10(17)6-14(15)20-13/h1-7,16-18H'
    ),
    flipped_scaffold_template(
        'InChI=1S/C15H12O5/c16-9-3-1-8(2-4-9)13-7-12(19)15-11(18)5-10(17)6-14(15)20-13/h1-6,13,16-18H,7H2'
    ),
    flipped_scaffold_template(
        'InChI=1S/C15H10O5/c16-9-3-1-8(2-4-9)11-7-20-13-6-10(17)5-12(18)14(13)15(11)19/h1-7,16-18H'
    ),
    flipped_scaffold_template(
        'InChI=1S/C15H14O/c1-2-6-12(7-3-1)15-11-10-13-8-4-5-9-14(13)16-15/h1-9,15H,10-11H2',
        FLAVONOID_CLOCKWISE_ROTATION
    ),
]

def orient_to_template(mol):
    for template in alignment_templates:
        if mol.HasSubstructMatch(template):
            rdDepictor.GenerateDepictionMatching2DStructure(mol, template)
            return

def draw_from_inchi(inchi_str, name, output_path):
    target_bond_length = 22
    canvas_padding = 3 * target_bond_length
    min_canvas_size = 220

    mol = MolFromInchi(inchi_str)
    rdDepictor.Compute2DCoords(mol)
    orient_to_template(mol)

    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]

    mean_bond_length = MeanBondLength(mol)
    if mean_bond_length <= 0:
        mean_bond_length = 1.5

    scale = target_bond_length / mean_bond_length
    width = max(
        min_canvas_size,
        int((max(xs) - min(xs)) * scale + 2 * canvas_padding),
    )
    height = max(
        min_canvas_size,
        int((max(ys) - min(ys)) * scale + 2 * canvas_padding),
    )

    drawer = MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 4.0
    opts.fixedBondLength = target_bond_length
    opts.fixedFontSize = 28

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