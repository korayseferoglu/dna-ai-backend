from flask import Flask, request, jsonify
from textwrap import wrap

app = Flask(__name__)

# --- CORS: index.html dosyasından erişim için ---
@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "POST, OPTIONS"
    return response


# ---------------- GENETİK KOD TABLOSU ----------------
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

# Amino asit sınıfları (hidrofobik / polar / yüklü)
HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQCH")
CHARGED = set("KRDE")


def aa_class(a: str) -> str:
    """AA'yı sınıfa çevir (UI renklendirme ve AI etkisi için)."""
    if a in HYDROPHOBIC:
        return "hydrophobic"
    if a in POLAR:
        return "polar"
    if a in CHARGED:
        return "charged"
    return "other"


# ---------------- TEMEL FONKSİYONLAR ----------------

def clean_dna(seq: str) -> str:
    return "".join([b for b in seq.upper() if b in "ATGC"])


def dna_to_mrna(dna: str) -> str:
    return dna.replace("T", "U")


def translate_mrna(mrna: str) -> str:
    """mRNA → protein. İlk AUG'dan başla, stop'ta dur."""
    start = mrna.find("AUG")
    if start == -1:
        return ""
    protein = []
    for i in range(start, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]
        aa = CODON_TABLE.get(codon, "?")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def protein_domains(protein_len: int):
    """
    Örnek domain haritası:
      - Eğer uzunluk ≥ 6 ise: 2–4 arası 'Catalytic core', geri kalanı 'Regulatory tail'
      - Aksi durumda tek 'Peptide region'
    Bu tamamen görselleştirme için, gerçek Pfam yerine educational amaçlı.
    """
    if protein_len <= 0:
        return []
    if protein_len >= 6:
        return [
            {"start": 2, "end": 4, "name": "Catalytic core"},
            {"start": 5, "end": protein_len, "name": "Regulatory tail"},
        ]
    else:
        return [{"start": 1, "end": protein_len, "name": "Peptide region"}]


def translate_dna(dna: str) -> dict:
    """DNA → mRNA → protein + AA sınıfları + domain haritası."""
    dna_clean = clean_dna(dna)
    mrna = dna_to_mrna(dna_clean)
    protein = translate_mrna(mrna)
    codons = wrap(mrna, 3)
    aa_props = [aa_class(a) for a in protein]
    domains = protein_domains(len(protein))

    return {
        "clean_dna": dna_clean,
        "mrna": mrna,
        "codons": codons,
        "protein": protein,
        "aa_properties": aa_props,
        "domains": domains,
    }


def apply_point_mutation(dna: str, index_1based: int, new_base: str) -> str:
    dna_list = list(clean_dna(dna))
    i = index_1based - 1
    if i < 0 or i >= len(dna_list):
        raise ValueError("Pozisyon DNA uzunluğu dışında.")
    dna_list[i] = new_base.upper()
    return "".join(dna_list)


def classify_mutation(wt_protein: str, mut_protein: str) -> dict:
    """
    WT vs mutant protein karşılaştırması:
      - Aynı ise → silent
      - Mutant daha kısa ise → truncation (nonsense veya frameshift+nonsense)
            → AA değişimi: Xn* şeklinde raporlanır
      - Mutant daha uzun ise → insertion/frameshift benzeri
      - Uzunluk aynı ve fark varsa → missense
    """
    # Tamamen aynı → silent
    if wt_protein == mut_protein:
        return {
            "type": "silent",
            "position": None,
            "aa_from": None,
            "aa_to": None,
        }

    len_wt = len(wt_protein)
    len_mut = len(mut_protein)
    min_len = min(len_wt, len_mut)

    # İlk farklı AA pozisyonunu bul
    first_diff = None
    for i in range(min_len):
        if wt_protein[i] != mut_protein[i]:
            first_diff = i  # 0-based
            break

    # --- 1) Mutant daha kısa → erken STOP (truncation) ---
    if len_mut < len_wt:
        # Eğer tüm prefix aynıysa: WT’nin ilk len_mut AA’sı mutantla aynı
        # demek ki STOP, bir sonraki pozisyonda oluşmuş
        if first_diff is None:
            pos = len_mut + 1           # 1-based
        else:
            # Hem AA değişimi hem erken STOP olabilir,
            # yine de etkiyi STOP üzerinden raporlayalım
            pos = first_diff + 1

        aa_from = wt_protein[pos - 1]   # WT’de bu AA olacaktı
        aa_to = "*"                     # Mutantta burada STOP var (nonsense)
        return {
            "type": "nonsense",
            "position": pos,
            "aa_from": aa_from,
            "aa_to": aa_to,
        }

    # --- 2) Mutant daha uzun → insertion / frameshift benzeri ---
    if len_mut > len_wt:
        if first_diff is None:
            # WT tamamen prefix, sonra uzama var
            pos = len_wt + 1
            aa_from = None
            aa_to = mut_protein[pos - 1]
        else:
            pos = first_diff + 1
            aa_from = wt_protein[pos - 1]
            aa_to = mut_protein[pos - 1]

        return {
            "type": "frameshift_or_insertion",
            "position": pos,
            "aa_from": aa_from,
            "aa_to": aa_to,
        }

    # --- 3) Uzunluklar aynı → klasik missense ---
    if first_diff is not None:
        pos = first_diff + 1
        aa_from = wt_protein[pos - 1]
        aa_to = mut_protein[pos - 1]
        return {
            "type": "missense",
            "position": pos,
            "aa_from": aa_from,
            "aa_to": aa_to,
        }

    # Teorik olarak buraya düşmez ama güvence olsun
    return {
        "type": "unknown",
        "position": None,
        "aa_from": None,
        "aa_to": None,
    }


def predict_impact(mut_info: dict) -> str:
    """
    Mini "AI" uzman sistemi:
      - silent → benign
      - nonsense, frameshift_or_insertion → probably_damaging
      - missense:
            aynı sınıf → possibly_benign
            farklı sınıf → possibly_damaging
    """
    mtype = mut_info["type"]

    if mtype == "silent":
        return "benign"

    if mtype in ("nonsense", "frameshift_or_insertion"):
        return "probably_damaging"

    aa_from = mut_info["aa_from"]
    aa_to = mut_info["aa_to"]
    if not aa_from or not aa_to:
        return "unknown"

    cls_from = aa_class(aa_from)
    cls_to = aa_class(aa_to)

    if cls_from == cls_to:
        return "possibly_benign"
    else:
        return "possibly_damaging"



def domain_effect(mut_info: dict, domains: list) -> dict:
    """Mutasyon domain içinde mi, hangi domain?"""
    pos = mut_info.get("position")
    if not pos:
        return {"in_domain": False, "domain_name": None}

    for d in domains:
        if d["start"] <= pos <= d["end"]:
            return {"in_domain": True, "domain_name": d["name"]}
    return {"in_domain": False, "domain_name": None}


# ---------------- API ENDPOINTLERİ ----------------

@app.route("/api/translate", methods=["POST"])
def api_translate():
    data = request.get_json()
    dna = data.get("dna", "")
    result = translate_dna(dna)
    return jsonify(result)


@app.route("/api/mutate", methods=["POST"])
def api_mutate():
    data = request.get_json()
    dna = data.get("dna", "")
    pos = int(data.get("position", 1))
    new_base = data.get("new_base", "A")

    try:
        dna_mut = apply_point_mutation(dna, pos, new_base)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    wt = translate_dna(dna)
    mut = translate_dna(dna_mut)

    mut_info = classify_mutation(wt["protein"], mut["protein"])
    impact = predict_impact(mut_info)
    dom_eff = domain_effect(mut_info, wt["domains"])

    response = {
        "wt": wt,
        "mut": mut,
        "mutation_info": mut_info,
        "impact": impact,
        "mutated_dna": dna_mut,
        "domain_effect": dom_eff,
    }
    return jsonify(response)


if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)

