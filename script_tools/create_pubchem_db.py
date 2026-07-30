import urllib.request
from pathlib import Path
import gzip
import shutil
import sqlite3
import hashlib

pubchem_dir = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/"
database_dir = "/Volumes/bluecub/databases/pubchem"
database_name = "260730_pubchem.db"

# built on 260730 from data fetched on 260730

files = [
    "CID-InChI-Key.gz",
    "CID-Mass.gz",
    "CID-Title.gz",
    "CID-SMILES.gz"
]

for i in files:
    filenm = f"{database_dir}/{i}"
    md5nm = f"{filenm}.md5"
    if not Path(filenm).exists():
        urllib.request.urlretrieve(
            url=f"{pubchem_dir}/{i}",
            filename=filenm
        )
    if not Path(md5nm).exists():
        urllib.request.urlretrieve(
            url=f"{pubchem_dir}/{i}.md5",
            filename=md5nm
        )
    with open(f"{database_dir}/{i}.md5", "r") as f:
        md5_pubchem = f.read().split()[0]

    with open(f"{database_dir}/{i}", "rb") as k:
        md5_dl = hashlib.md5(k.read()).hexdigest()

    if md5_dl != md5_pubchem:
        raise ValueError(f"md5's do not match for: {database_dir}/{i}.md5")
        

con = sqlite3.connect(f"{database_dir}/{database_name}")
cur = con.cursor()

cur.execute("""
    CREATE TABLE IF NOT EXISTS compounds (
        cid INTEGER PRIMARY KEY,
        inchikey TEXT,
        formula TEXT,
        monoisotopic_mass REAL,
        exact_mass REAL,
        title TEXT,
        smiles TEXT
    )
""")
con.commit()

BATCH_SIZE = 200_000

# CID-InChI-Key.gz
batch = []
with gzip.open(f"{database_dir}/CID-InChI-Key.gz", "rt") as fh:
    for line in fh:
        cid, _inchi, inchikey = line.rstrip("\n").split("\t")
        batch.append((int(cid), inchikey))
        if len(batch) >= BATCH_SIZE:
            cur.executemany(
                """
                INSERT INTO compounds (cid, inchikey) VALUES (?, ?)
                ON CONFLICT(cid) DO UPDATE SET inchikey = excluded.inchikey
                """,
                batch
            )
            con.commit()
            batch = []
if batch:
    cur.executemany(
        """
        INSERT INTO compounds (cid, inchikey) VALUES (?, ?)
        ON CONFLICT(cid) DO UPDATE SET inchikey = excluded.inchikey
        """,
        batch
    )
    con.commit()

# CID-Mass.gz
batch = []
with gzip.open(f"{database_dir}/CID-Mass.gz", "rt") as fh:
    for line in fh:
        cid, formula, mono_mass, exact_mass = line.rstrip("\n").split("\t")
        batch.append((int(cid), formula, float(mono_mass), float(exact_mass)))
        if len(batch) >= BATCH_SIZE:
            cur.executemany(
                """
                INSERT INTO compounds (cid, formula, monoisotopic_mass, exact_mass)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(cid) DO UPDATE SET
                    formula = excluded.formula,
                    monoisotopic_mass = excluded.monoisotopic_mass,
                    exact_mass = excluded.exact_mass
                """,
                batch
            )
            con.commit()
            batch = []
if batch:
    cur.executemany(
        """
        INSERT INTO compounds (cid, formula, monoisotopic_mass, exact_mass)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(cid) DO UPDATE SET
            formula = excluded.formula,
            monoisotopic_mass = excluded.monoisotopic_mass,
            exact_mass = excluded.exact_mass
        """,
        batch
    )
    con.commit()

# CID-Title.gz
batch = []
with gzip.open(f"{database_dir}/CID-Title.gz", "rt") as fh:
    for line in fh:
        cid, title = line.rstrip("\n").split("\t", maxsplit=1)
        batch.append((int(cid), title))
        if len(batch) >= BATCH_SIZE:
            cur.executemany(
                """
                INSERT INTO compounds (cid, title) VALUES (?, ?)
                ON CONFLICT(cid) DO UPDATE SET title = excluded.title
                """,
                batch
            )
            con.commit()
            batch = []
if batch:
    cur.executemany(
        """
        INSERT INTO compounds (cid, title) VALUES (?, ?)
        ON CONFLICT(cid) DO UPDATE SET title = excluded.title
        """,
        batch
    )
    con.commit()

# CID-SMILES.gz
batch = []
with gzip.open(f"{database_dir}/CID-SMILES.gz", "rt") as fh:
    for line in fh:
        cid, smiles = line.rstrip("\n").split("\t", maxsplit=1)
        batch.append((int(cid), smiles))
        if len(batch) >= BATCH_SIZE:
            cur.executemany(
                """
                INSERT INTO compounds (cid, smiles) VALUES (?, ?)
                ON CONFLICT(cid) DO UPDATE SET smiles = excluded.smiles
                """,
                batch
            )
            con.commit()
            batch = []
if batch:
    cur.executemany(
        """
        INSERT INTO compounds (cid, smiles) VALUES (?, ?)
        ON CONFLICT(cid) DO UPDATE SET smiles = excluded.smiles
        """,
        batch
    )
    con.commit()

cur.execute("CREATE INDEX IF NOT EXISTS idx_inchikey ON compounds (inchikey)")
con.commit()