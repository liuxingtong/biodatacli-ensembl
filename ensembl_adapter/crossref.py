"""Tiny readers for the PDB side of the cross-database check.

Only what is needed for the 4DNK / P31946 example: C-alpha observations from
an mmCIF file and the SIFTS UniProt<->PDB residue mapping.
"""

from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from typing import Dict, Iterable, List, Optional, Set

NS = {"s": "http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"}


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _atom_columns(path: str):
    with open_text(path) as handle:
        cols = []
        rows = []
        in_atom = False
        for line in handle:
            if line.startswith("_atom_site."):
                cols.append(line.split(".", 1)[1].strip())
                in_atom = True
                continue
            if not in_atom:
                continue
            if line.startswith("#") or line.startswith("loop_"):
                if rows:
                    break
                continue
            if line.startswith("_") or line.startswith("data_"):
                break
            if not line.strip():
                continue
            rows.append(line.rstrip("\n").split())
    return cols, rows


def load_ca_positions(path: str, model: str = "1") -> Dict[str, Dict[int, Dict[str, object]]]:
    """Return {chain: {auth_seq_id: {label_seq_id, comp, altloc}}} for C-alpha."""
    cols, rows = _atom_columns(path)
    idx = {name: i for i, name in enumerate(cols)}
    out: Dict[str, Dict[int, Dict[str, object]]] = {}
    for row in rows:
        if len(row) < len(cols):
            continue
        if row[idx["group_PDB"]] != "ATOM":
            continue
        if row[idx["pdbx_PDB_model_num"]] != model:
            continue
        if row[idx["label_atom_id"]].strip('"') != "CA":
            continue
        chain = row[idx["auth_asym_id"]]
        try:
            auth = int(row[idx["auth_seq_id"]])
        except ValueError:
            continue
        altloc = row[idx["label_alt_id"]]
        entry = out.setdefault(chain, {})
        # Prefer the blank/primary altloc; otherwise keep the first.
        if auth in entry and altloc not in (".", "?", "A"):
            continue
        entry[auth] = {
            "label_seq_id": row[idx["label_seq_id"]],
            "comp": row[idx["label_comp_id"]],
            "altloc": altloc,
            "ins_code": row[idx["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in idx else "?",
        }
    return out


def load_sifts(path: str, uniprot_acc: str) -> Dict[str, Dict[int, Dict[str, object]]]:
    """Return {chain: {uniprot_pos: {pdb_auth, pdb_label, resname, evidence}}}."""
    with open_text(path) as handle:
        try:
            root = ET.fromstring(handle.read())
        except ET.ParseError:
            handle.seek(0)
            root = ET.fromstring(handle.read())
    out: Dict[str, Dict[int, Dict[str, object]]] = {}
    for entity in root.findall("s:entity", NS):
        chain = entity.get("entityId")
        for residue in entity.iter("{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residue"):
            up = None
            pdb = None
            for cross in residue.findall("s:crossRefDb", NS):
                if cross.get("dbSource") == "UniProt" and cross.get("dbAccessionId") == uniprot_acc:
                    up = cross
                elif cross.get("dbSource") == "PDB":
                    pdb = cross
            if up is None:
                continue
            try:
                up_pos = int(up.get("dbResNum"))
            except (TypeError, ValueError):
                continue
            pdb_auth = pdb.get("dbResNum") if pdb is not None else "null"
            out.setdefault(chain, {})[up_pos] = {
                "pdb_auth": None if pdb_auth in (None, "null", "?") else int(pdb_auth),
                "pdb_label": int(residue.get("dbResNum")),
                "resname": residue.get("dbResName"),
                "evidence": up.get("dbEvidence", ""),
            }
    return out


def sifts_metadata(path: str) -> Dict[str, str]:
    """Pull the db version block from the SIFTS XML header."""
    with open_text(path) as handle:
        head = handle.read(6000)
    meta = {}
    for token in ("dbVersion", "date"):
        start = head.find(f'{token}="')
        if start >= 0:
            start += len(token) + 2
            end = head.find('"', start)
            meta[token] = head[start:end]
    return meta
