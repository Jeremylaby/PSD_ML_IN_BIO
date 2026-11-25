from predictor.fasta import read_fasta
from utils.cif_to_fasta import entity_sequences, write_fasta

# Minimal mmCIF: entity 1 is a protein with a multi-line sequence and a
# selenomethionine (MSE -> M in the canonical code), entity 2 is DNA.
CIF = """data_TEST
#
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code
_entity_poly.pdbx_seq_one_letter_code_can
1 'polypeptide(L)'
;MKT(MSE)AYIA
KQR
;
;MKTMAYIA
KQR
;
2 polydeoxyribonucleotide '(DA)(DC)' AC
#
"""


def test_entity_sequences(tmp_path):
    cif = tmp_path / "TEST.cif"
    cif.write_text(CIF)
    assert entity_sequences(cif) == {"1": "MKTMAYIAKQR"}  # DNA skipped, lines joined


def test_write_fasta_round_trip(tmp_path):
    out = tmp_path / "out" / "x.fasta"
    records = [("10ZN_1", "A" * 130), ("1ABC_2", "MKV")]
    write_fasta(out, records)
    assert read_fasta(out) == records
    assert max(len(l) for l in out.read_text().splitlines()) == 60
