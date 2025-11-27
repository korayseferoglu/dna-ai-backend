from flask import Flask, request, jsonify
from textwrap import wrap

app = Flask(__name__)

@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "POST, OPTIONS"
    return response


# ---- Genetic code (mRNA codons) ----
CODON_TABLE = {
    "AUG": "M",
    "UAA": "*", "UAG": "*", "UGA": "*",

    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",

    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    "UAU": "Y", "UAC": "Y",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",

    "UGU": "C", "UGC": "C",
    "UGG": "W",

    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DNA_STOP = {"TAA", "TAG", "TGA"}
MRNA_STOP = {"UAA", "UAG", "UGA"}

HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQCH")
CHARGED = set("KRDE")


def aa_class(a: str) -> str:
    if a in HYDROPHOBIC:
        return "hydrophobic"
    if a in POLAR:
        return "polar"
    if a in CHARGED:
        return "charged"
    return "other"


def clean_dna(seq: str) -> str:
    return "".join([b for b in (seq or "").upper() if b in "ATGC"])


def dna_to_mrna(dna: str) -> str:
    return dna.replace("T", "U")


def dna_codon_to_mrna_codon(dna_codon: str) -> str:
    return dna_codon.replace("T", "U")


def find_first_atg(dna: str) -> int:
    return dna.find("ATG")


def find_next_inframe_atg(dna: str, start: int) -> int:
    """start'tan itibaren aynı frame (start, start+3, start+6...) içinde ATG ara."""
    if start < 0:
        return -1
    for i in range(start, len(dna) - 2, 3):
        if dna[i:i+3] == "ATG":
            return i
    return -1


def translate_from_start_dna(dna: str, start: int) -> dict:
    """
    Coding-strand DNA varsayımıyla (ATG start / TAA-TAG-TGA stop),
    verilen start indeksinden (0-based) çeviri yapar.
    STOP görüldüğünde durur.
    """
    if start < 0 or start > len(dna) - 3:
        return {
            "start_dna_index": start,
            "stop_dna_index": None,
            "stop_codon": None,
            "ended_by_stop": False,
            "codons_in_frame_dna": [],
            "codons_in_frame_mrna": [],
            "protein": "",
        }

    codons_dna = []
    codons_mrna = []
    protein = []
    stop_idx = None
    stop_codon = None

    for i in range(start, len(dna) - 2, 3):
        codon_dna = dna[i:i+3]
        codon_mrna = dna_codon_to_mrna_codon(codon_dna)
        codons_dna.append(codon_dna)
        codons_mrna.append(codon_mrna)

        aa = CODON_TABLE.get(codon_mrna, "?")
        if aa == "*":
            stop_idx = i
            stop_codon = codon_dna
            break
        protein.append(aa)

    return {
        "start_dna_index": start,
        "stop_dna_index": stop_idx,              # stop codon start index (0-based, DNA)
        "stop_codon": stop_codon,                # DNA stop codon (TAA/TAG/TGA)
        "ended_by_stop": stop_idx is not None,
        "codons_in_frame_dna": codons_dna,       # includes stop codon if encountered
        "codons_in_frame_mrna": codons_mrna,
        "protein": "".join(protein),
    }


def protein_domains(protein_len: int):
    if protein_len <= 0:
        return []
    if protein_len >= 6:
        return [
            {"start": 2, "end": 4, "name": "Catalytic core"},
            {"start": 5, "end": protein_len, "name": "Regulatory tail"},
        ]
    return [{"start": 1, "end": protein_len, "name": "Peptide region"}]


def translate_dna(dna: str, forced_start_dna_index: int | None = None) -> dict:
    dna_clean = clean_dna(dna)
    mrna = dna_to_mrna(dna_clean)

    if forced_start_dna_index is None:
        start = find_first_atg(dna_clean)
    else:
        start = forced_start_dna_index

    tr = translate_from_start_dna(dna_clean, start)

    protein = tr["protein"]
    aa_props = [aa_class(a) for a in protein]
    domains = protein_domains(len(protein))

    return {
        "clean_dna": dna_clean,
        "mrna": mrna,
        "codons_all_mrna": wrap(mrna, 3),
        "start_dna_index": tr["start_dna_index"],
        "stop_dna_index": tr["stop_dna_index"],
        "stop_codon": tr["stop_codon"],
        "ended_by_stop": tr["ended_by_stop"],
        "codons_in_frame_dna": tr["codons_in_frame_dna"],
        "codons_in_frame_mrna": tr["codons_in_frame_mrna"],
        "protein": protein,
        "aa_properties": aa_props,
        "domains": domains,
    }


def apply_point_mutation(dna: str, index_1based: int, new_base: str) -> dict:
    dna_clean = clean_dna(dna)
    if not dna_clean:
        raise ValueError("DNA boş veya geçersiz.")

    new_base = (new_base or "").upper()
    if new_base not in "ATGC":
        raise ValueError("new_base sadece A/T/G/C olabilir.")

    i = index_1based - 1
    if i < 0 or i >= len(dna_clean):
        raise ValueError("Pozisyon DNA uzunluğu dışında.")

    old_base = dna_clean[i]
    if old_base == new_base:
        return {"mutated_dna": dna_clean, "old_base": old_base, "changed": False}

    dna_list = list(dna_clean)
    dna_list[i] = new_base
    return {"mutated_dna": "".join(dna_list), "old_base": old_base, "changed": True}


def classify_point_mutation(wt: dict, mut: dict, dna_pos_1based: int) -> dict:
    """
    Biyolojik sınıflandırma: tek baz değişimi -> WT ORF frame referans alınır.
    """
    dna_index0 = dna_pos_1based - 1
    wt_start = wt.get("start_dna_index", -1)

    if wt_start is None or wt_start < 0:
        return {"type": "no_start_codon", "position": None, "aa_from": None, "aa_to": None}

    # WT stop varsa coding sınırı: stop codon dahil
    wt_stop = wt.get("stop_dna_index", None)

    # Bölge tayini (WT ORF bazlı)
    if dna_index0 < wt_start:
        return {"type": "noncoding", "region": "5_UTR_or_upstream", "position": None, "aa_from": None, "aa_to": None}

    if wt_stop is not None and dna_index0 >= wt_stop + 3:
        return {"type": "noncoding", "region": "3_UTR_or_downstream", "position": None, "aa_from": None, "aa_to": None}

    # Frame içinde hangi codon?
    rel = dna_index0 - wt_start
    codon_index0 = rel // 3
    pos_in_codon = (rel % 3) + 1  # 1..3
    codon_start = wt_start + codon_index0 * 3

    wt_codon_dna = wt["clean_dna"][codon_start:codon_start+3]
    mut_codon_dna = mut["clean_dna"][codon_start:codon_start+3]
    wt_codon_mrna = dna_codon_to_mrna_codon(wt_codon_dna)
    mut_codon_mrna = dna_codon_to_mrna_codon(mut_codon_dna)

    wt_aa = CODON_TABLE.get(wt_codon_mrna, "?")
    mut_aa = CODON_TABLE.get(mut_codon_mrna, "?")
    aa_pos_1based = codon_index0 + 1

    # Start codon bozuldu mu?
    if aa_pos_1based == 1 and wt_codon_dna == "ATG" and mut_codon_dna != "ATG":
        return {
            "type": "start_lost",
            "position": 1,
            "aa_from": "M",
            "aa_to": None,
            "details": {"wt_codon": wt_codon_dna, "mut_codon": mut_codon_dna, "pos_in_codon": pos_in_codon},
        }

    # Silent
    if wt_aa == mut_aa:
        return {
            "type": "silent",
            "position": aa_pos_1based,
            "aa_from": wt_aa,
            "aa_to": mut_aa,
            "details": {"wt_codon": wt_codon_dna, "mut_codon": mut_codon_dna, "pos_in_codon": pos_in_codon},
        }

    # Stop gain = nonsense
    if mut_aa == "*" and wt_aa != "*":
        return {
            "type": "nonsense",
            "position": aa_pos_1based,
            "aa_from": wt_aa,
            "aa_to": "*",
            "details": {"wt_codon": wt_codon_dna, "mut_codon": mut_codon_dna, "pos_in_codon": pos_in_codon},
        }

    # Stop loss
    if wt_aa == "*" and mut_aa != "*":
        return {
            "type": "stop_loss",
            "position": aa_pos_1based,
            "aa_from": "*",
            "aa_to": mut_aa,
            "details": {"wt_codon": wt_codon_dna, "mut_codon": mut_codon_dna, "pos_in_codon": pos_in_codon},
        }

    # Missense
    return {
        "type": "missense",
        "position": aa_pos_1based,
        "aa_from": wt_aa,
        "aa_to": mut_aa,
        "details": {"wt_codon": wt_codon_dna, "mut_codon": mut_codon_dna, "pos_in_codon": pos_in_codon},
    }


def predict_impact(mut_info: dict) -> str:
    t = mut_info.get("type")
    if t in ("silent",):
        return "benign"
    if t in ("nonsense", "stop_loss", "start_lost"):
        return "probably_damaging"
    if t == "missense":
        a, b = mut_info.get("aa_from"), mut_info.get("aa_to")
        if not a or not b or a == "?" or b == "?":
            return "unknown"
        return "possibly_benign" if aa_class(a) == aa_class(b) else "possibly_damaging"
    return "unknown"


def domain_effect(mut_info: dict, domains: list) -> dict:
    pos = mut_info.get("position")
    if not pos or not domains:
        return {"in_domain": False, "domain_name": None}
    for d in domains:
        if d["start"] <= pos <= d["end"]:
            return {"in_domain": True, "domain_name": d["name"]}
    return {"in_domain": False, "domain_name": None}


@app.route("/api/translate", methods=["POST"])
def api_translate():
    data = request.get_json() or {}
    dna = data.get("dna", "")
    return jsonify(translate_dna(dna))


@app.route("/api/mutate", methods=["POST"])
def api_mutate():
    data = request.get_json() or {}
    dna = data.get("dna", "")
    pos = int(data.get("position", 1))     # DNA baz pozisyonu (1-based)
    new_base = data.get("new_base", "A")

    wt = translate_dna(dna)

    try:
        mut_res = apply_point_mutation(dna, pos, new_base)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    dna_mut = mut_res["mutated_dna"]

    # Mutant: WT start’ı sabitle (frame tutarlı olsun)
    mut = translate_dna(dna_mut, forced_start_dna_index=wt["start_dna_index"])

    # Eğer WT start bozulduysa, alternatif start (aynı frame) arayıp ayrı protein de verelim
    alt_start_used = False
    alt = None
    wt_start = wt.get("start_dna_index", -1)
    if wt_start is not None and wt_start >= 0:
        if mut["clean_dna"][wt_start:wt_start+3] != "ATG":
            next_start = find_next_inframe_atg(mut["clean_dna"], wt_start + 3)
            if next_start != -1:
                alt = translate_dna(dna_mut, forced_start_dna_index=next_start)
                alt_start_used = True

    mut_info = classify_point_mutation(wt, mut, pos)
    impact = predict_impact(mut_info)
    dom_eff = domain_effect(mut_info, wt["domains"])

    aa_change = None
    if mut_info.get("position") and mut_info.get("aa_from") is not None and mut_info.get("aa_to") is not None:
        aa_change = f"{mut_info['aa_from']}{mut_info['position']}{mut_info['aa_to']}"
    elif mut_info.get("type") == "start_lost":
        aa_change = f"M1?"

    response = {
        "wt": wt,
        "mut": mut,
        "mutation_info": {**mut_info, "aa_change": aa_change},
        "impact": impact,
        "mutated_dna": dna_mut,
        "domain_effect": dom_eff,
        "debug": {
            "changed": mut_res["changed"],
            "old_base": mut_res["old_base"],
            "dna_position_1based": pos,
            "dna_index_0based": pos - 1,
            "wt_start_dna_index": wt.get("start_dna_index"),
            "wt_stop_dna_index": wt.get("stop_dna_index"),
            "alt_start_used": alt_start_used,
        },
        "alt_orf_same_frame": alt,  # start bozulduysa burada alternatif ORF proteinini de görürsün
    }
    return jsonify(response)


if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)


