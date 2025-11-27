from flask import Flask, request, jsonify
from textwrap import wrap

app = Flask(__name__)

@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "POST, OPTIONS"
    return response


# ---------------- GENETİK KOD TABLOSU (mRNA codons) ----------------
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

STOP_CODONS = {"UAA", "UAG", "UGA"}

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
    if a == "*":
        return "stop"
    return "other"


def clean_dna(seq: str) -> str:
    return "".join([b for b in seq.upper() if b in "ATGC"])


def dna_to_mrna(dna: str) -> str:
    return dna.replace("T", "U")


def apply_point_mutation(dna: str, index_1based: int, new_base: str) -> str:
    new_base = new_base.upper()
    if new_base not in "ATGC":
        raise ValueError("new_base sadece A/T/G/C olabilir.")

    dna_list = list(clean_dna(dna))
    i = index_1based - 1
    if i < 0 or i >= len(dna_list):
        raise ValueError("Pozisyon DNA uzunluğu dışında.")
    dna_list[i] = new_base
    return "".join(dna_list)


def find_first_aug(mrna: str) -> int:
    return mrna.find("AUG")


def find_next_inframe_aug(mrna: str, start: int) -> int:
    """start'tan itibaren aynı frame içinde AUG ara (start, start+3, start+6...)."""
    if start < 0:
        return -1
    for i in range(start, len(mrna) - 2, 3):
        if mrna[i:i+3] == "AUG":
            return i
    return -1


def translate_from_start(mrna: str, start: int) -> dict:
    """
    Belirli bir start indeksinden (0-based, mrna üzerinde) translate et.
    Stop codonu görünce DURUR.
    """
    if start < 0 or start > len(mrna) - 3:
        return {
            "start": start,
            "stop": None,
            "stop_codon": None,
            "ended_by_stop": False,
            "codons_in_frame": [],
            "protein": "",
        }

    codons = []
    prot = []
    stop_idx = None
    stop_codon = None

    for i in range(start, len(mrna) - 2, 3):
        codon = mrna[i:i+3]
        codons.append(codon)
        aa = CODON_TABLE.get(codon, "?")
        if aa == "*":
            stop_idx = i
            stop_codon = codon
            break
        prot.append(aa)

    return {
        "start": start,
        "stop": stop_idx,                   # mrna index (0-based) of STOP codon start
        "stop_codon": stop_codon,
        "ended_by_stop": stop_idx is not None,
        "codons_in_frame": codons,           # includes stop codon if encountered
        "protein": "".join(prot),
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


def translate_dna(dna: str, forced_start: int | None = None) -> dict:
    dna_clean = clean_dna(dna)
    mrna = dna_to_mrna(dna_clean)

    if forced_start is None:
        start = find_first_aug(mrna)
    else:
        start = forced_start

    tr = translate_from_start(mrna, start)

    # UI için: tüm mrna'yı 3'lü parçalamak yerine,
    # frame-codons ayrıca verelim
    protein = tr["protein"]
    aa_props = [aa_class(a) for a in protein]
    domains = protein_domains(len(protein))

    return {
        "clean_dna": dna_clean,
        "mrna": mrna,
        "codons_all": wrap(mrna, 3),
        "codons_in_frame": tr["codons_in_frame"],
        "protein": protein,
        "aa_properties": aa_props,
        "domains": domains,
        "start_mrna_index": tr["start"],
        "stop_mrna_index": tr["stop"],
        "stop_codon": tr["stop_codon"],
        "ended_by_stop": tr["ended_by_stop"],
    }


def classify_point_mutation_biologically(wt: dict, mut: dict, dna_pos_1based: int) -> dict:
    """
    Biyolojik sınıflandırma: point mutation -> codon/AA bazında karar.
    WT'nin start+frame’i referans alınır.
    """
    mrna_index = dna_pos_1based - 1  # sense DNA varsayımı (T->U), indeks birebir gider

    if wt["start_mrna_index"] is None or wt["start_mrna_index"] < 0:
        return {
            "type": "no_start_codon",
            "region": "unknown",
            "position": None,
            "aa_from": None,
            "aa_to": None,
        }

    start = wt["start_mrna_index"]

    # coding bölgenin WT sınırları (stop kodonu dahil):
    # coding codonları start'tan başlar, stop varsa stop_codon(start..start+2)
    stop = wt["stop_mrna_index"]

    if mrna_index < start:
        return {
            "type": "noncoding",
            "region": "5_UTR_or_upstream",
            "position": None,
            "aa_from": None,
            "aa_to": None,
        }

    # stop yoksa coding’i mrna sonuna kadar kabul ediyoruz (basit model)
    if stop is not None and mrna_index >= stop + 3:
        return {
            "type": "noncoding",
            "region": "3_UTR_or_downstream",
            "position": None,
            "aa_from": None,
            "aa_to": None,
        }

    # frame içinde hangi codon?
    rel = mrna_index - start
    codon_index_0 = rel // 3            # 0-based codon number from start
    codon_pos_in_codon = (rel % 3) + 1  # 1..3

    codon_start = start + codon_index_0 * 3
    wt_codon = wt["mrna"][codon_start:codon_start+3]
    mut_codon = mut["mrna"][codon_start:codon_start+3]

    wt_aa = CODON_TABLE.get(wt_codon, "?")
    mut_aa = CODON_TABLE.get(mut_codon, "?")

    aa_pos_1based = codon_index_0 + 1

    # Start codonu etkisi (AA1)
    if aa_pos_1based == 1 and wt_codon == "AUG" and mut_codon != "AUG":
        return {
            "type": "start_lost",
            "region": "start_codon",
            "position": 1,
            "aa_from": "M",
            "aa_to": None,
            "details": {
                "wt_codon": wt_codon,
                "mut_codon": mut_codon,
                "codon_pos_in_codon": codon_pos_in_codon,
            }
        }

    # Stop gain/loss & silent/missense
    if wt_aa == mut_aa:
        return {
            "type": "silent" if wt_aa != "*" else "stop_preserved",
            "region": "coding",
            "position": aa_pos_1based,
            "aa_from": wt_aa,
            "aa_to": mut_aa,
            "details": {
                "wt_codon": wt_codon,
                "mut_codon": mut_codon,
                "codon_pos_in_codon": codon_pos_in_codon,
            }
        }

    if mut_aa == "*" and wt_aa != "*":
        # gerçek nonsense = STOP GAIN
        return {
            "type": "nonsense",
            "region": "coding",
            "position": aa_pos_1based,
            "aa_from": wt_aa,
            "aa_to": "*",
            "details": {
                "wt_codon": wt_codon,
                "mut_codon": mut_codon,
                "codon_pos_in_codon": codon_pos_in_codon,
            }
        }

    if wt_aa == "*" and mut_aa != "*":
        # stop loss -> readthrough (uzama olabilir)
        return {
            "type": "stop_loss",
            "region": "coding",
            "position": aa_pos_1based,
            "aa_from": "*",
            "aa_to": mut_aa,
            "details": {
                "wt_codon": wt_codon,
                "mut_codon": mut_codon,
                "codon_pos_in_codon": codon_pos_in_codon,
            }
        }

    return {
        "type": "missense",
        "region": "coding",
        "position": aa_pos_1based,
        "aa_from": wt_aa,
        "aa_to": mut_aa,
        "details": {
            "wt_codon": wt_codon,
            "mut_codon": mut_codon,
            "codon_pos_in_codon": codon_pos_in_codon,
        }
    }


def predict_impact(mut_info: dict) -> str:
    t = mut_info["type"]
    if t in ("silent", "stop_preserved"):
        return "benign"
    if t in ("nonsense", "stop_loss", "start_lost"):
        return "probably_damaging"
    if t == "missense":
        a = mut_info.get("aa_from")
        b = mut_info.get("aa_to")
        if not a or not b or a == "?" or b == "?":
            return "unknown"
        return "possibly_benign" if aa_class(a) == aa_class(b) else "possibly_damaging"
    if t == "noncoding":
        return "unknown"
    return "unknown"


def domain_effect(mut_info: dict, domains: list) -> dict:
    pos = mut_info.get("position")
    if not pos or not domains:
        return {"in_domain": False, "domain_name": None}
    for d in domains:
        if d["start"] <= pos <= d["end"]:
            return {"in_domain": True, "domain_name": d["name"]}
    return {"in_domain": False, "domain_name": None}


# ---------------- API ENDPOINTLERİ ----------------

@app.route("/api/translate", methods=["POST"])
def api_translate():
    data = request.get_json() or {}
    dna = data.get("dna", "")
    result = translate_dna(dna)
    return jsonify(result)


@app.route("/api/mutate", methods=["POST"])
def api_mutate():
    data = request.get_json() or {}
    dna = data.get("dna", "")
    pos = int(data.get("position", 1))         # DNA baz pozisyonu (1-based)
    new_base = data.get("new_base", "A")

    try:
        dna_mut = apply_point_mutation(dna, pos, new_base)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    # WT'yi bul; WT start+frame referans olsun
    wt = translate_dna(dna)

    # Mutant çeviriyi aynı start ile yap (start bozulursa "start_lost" raporlarız)
    mut = translate_dna(dna_mut, forced_start=wt["start_mrna_index"])

    mut_info = classify_point_mutation_biologically(wt, mut, pos)
    impact = predict_impact(mut_info)
    dom_eff = domain_effect(mut_info, wt["domains"])

    # ekstra: mutant stop AA pozisyonu (visual/debug)
    mut_stop_aa_pos = (len(mut["protein"]) + 1) if mut["ended_by_stop"] else None

    response = {
        "wt": wt,
        "mut": mut,
        "mutation_info": mut_info,
        "impact": impact,
        "mutated_dna": dna_mut,
        "domain_effect": dom_eff,
        "debug": {
            "dna_position_1based": pos,
            "mrna_index_0based": pos - 1,
            "wt_start_mrna_index": wt["start_mrna_index"],
            "wt_stop_mrna_index": wt["stop_mrna_index"],
            "mut_stop_aa_position": mut_stop_aa_pos,
        }
    }
    return jsonify(response)


if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)


