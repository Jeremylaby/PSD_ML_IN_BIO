import torch
import pytest

from predictor.fasta import read_fasta
from predictor.model import DummyStructureModel
from predictor.pdb_writer import format_ca_atom
from predictor.predict import predict
from predictor.vocab import PAD_INDEX, UNK_INDEX, VOCAB_SIZE, encode_batch


def test_vocab_size():
    assert VOCAB_SIZE == 22


def test_encode_pads_truncates_and_masks():
    x, mask, lengths = encode_batch(["ACD", "A" * 10, "BZ"], max_len=5)
    assert x.shape == (3, 5, VOCAB_SIZE)
    assert lengths == [3, 5, 2]
    assert mask.sum(dim=1).tolist() == [3, 5, 2]
    assert x[0, 3:, PAD_INDEX].eq(1).all()
    assert x[2, :2, UNK_INDEX].eq(1).all()
    assert x.sum(dim=-1).eq(1).all()  # exactly one-hot everywhere


def test_model_shapes_and_masking():
    x, mask, lengths = encode_batch(["ACDE", "A" * 20], max_len=8)
    model = DummyStructureModel()
    coords = model(x, mask)
    assert coords.shape == (2, 8, 3)
    assert coords[0, 4:].eq(0).all()
    per_seq = model.predict(x, mask)
    assert [c.shape for c in per_seq] == [(4, 3), (8, 3)]


def test_model_rejects_wrong_shape():
    with pytest.raises(ValueError):
        DummyStructureModel()(torch.zeros(1, 4, 5), torch.ones(1, 4, dtype=torch.bool))


def test_atom_record_columns():
    line = format_ca_atom(1, "MET", "A", 1, 1.5, -22.25, 333.125)
    assert len(line) == 80
    assert line[0:6] == "ATOM  "
    assert line[6:11] == "    1"
    assert line[12:16] == " CA "
    assert line[17:20] == "MET"
    assert line[21] == "A"
    assert line[22:26] == "   1"
    assert float(line[30:38]) == 1.5
    assert float(line[38:46]) == -22.25
    assert float(line[46:54]) == 333.125
    assert line[54:60] == "  1.00"
    assert line[76:78] == " C"


def test_end_to_end(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">a\nACDEFG\n>b\nMKXBU\n>c\n" + "A" * 50 + "\n")
    out = tmp_path / "out.pdb"
    assert predict(str(fasta), str(out), max_len=20, batch_size=2, seed=0) == 3

    lines = out.read_text().splitlines()
    atoms = [l for l in lines if l.startswith("ATOM")]
    assert len(atoms) == 6 + 5 + 20  # min(len, D) per protein
    assert sum(l.startswith("MODEL") for l in lines) == 3
    assert lines[-1] == "END"
    assert all(len(l) == 80 and l[12:16] == " CA " for l in atoms)
    assert sum(l[17:20] == "UNK" for l in atoms) == 3  # X, B, U


def test_read_fasta(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">p1 desc\nACD\nEF\n>p2\nGH\n")
    assert read_fasta(fasta) == [("p1", "ACDEF"), ("p2", "GH")]


def test_make_run_dir_creates_unique_folders(tmp_path):
    from predictor.predict import make_run_dir

    first = make_run_dir(tmp_path)
    second = make_run_dir(tmp_path)  # same second -> must not reuse the folder
    assert first.is_dir() and second.is_dir()
    assert first != second
    assert first.name.startswith("run_")


def test_info_json(tmp_path):
    import json

    fasta = tmp_path / "in.fasta"
    # 7ZX8_2 is a real dataset ID (test split), "mine" is a custom sequence.
    fasta.write_text(">mine\nACDEFG\n>7ZX8_2\n" + "A" * 30 + "\n")
    out, info_path = tmp_path / "out.pdb", tmp_path / "info.json"
    predict(str(fasta), str(out), max_len=20, batch_size=1, info_path=str(info_path))

    info = json.loads(info_path.read_text())
    assert isinstance(info["parameters"]["seed"], int)  # random seed is recorded
    assert info["n_proteins"] == 2 and info["n_ca_atoms"] == 6 + 20
    mine, real = info["proteins"]
    assert mine["source"] == {"in_dataset": False}
    assert (mine["fasta_index"], mine["pdb_model"]) == (0, 1)
    assert (real["fasta_index"], real["pdb_model"], real["batch_index"]) == (1, 2, 1)
    assert real["truncated"] and real["sequence_length"] == 30 and real["processed_length"] == 20
    assert real["source"]["split"] == "test"
    assert real["source"]["pdb_id"] == "7ZX8" and real["source"]["entity_id"] == "2"


def test_same_seed_same_coordinates(tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">a\nACDEFGHIK\n")
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    predict(str(fasta), str(a), seed=42)
    predict(str(fasta), str(b), seed=42)
    assert a.read_text() == b.read_text()
