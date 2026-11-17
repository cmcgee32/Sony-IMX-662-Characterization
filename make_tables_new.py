#!/usr/bin/env python3
"""Booktabs tables for the IMX226 and IMX290 analyses and the three-sensor comparison.

usage: make_tables_new.py <IMX226_results> <IMX290_results> <IMX662 booktabs folder> <out>
Reads the CSVs written by imx226_tables.m / imx290_ecc_analysis.m and writes
.tex, .pdf, .png and .csv for every table (same style as make_tables.py).
"""
import csv, math, os, re, shutil, subprocess, sys
R226, R290, T662, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
ME, LU, BS = "Carissma (2026)", "Lupo (2025)", "Schneider et al. (2024)"
LINE = {"C K-alpha": "C K$\\alpha$", "F K-alpha": "F K$\\alpha$", "Mn K-alpha": "Mn K$\\alpha$", "Mn K-beta": "Mn K$\\beta$"}
T = []


def rd(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return list(csv.reader(f))


def num(x, d=1):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return esc(str(x))
    if math.isnan(v):
        return "---"
    return f"{v:,.{d}f}".replace(",", "{,}")


def sci(x, d=2):
    v = float(x)
    if v == 0 or math.isnan(v):
        return "---"
    e = int(math.floor(math.log10(abs(v))))
    return f"${v/10**e:.{d}f} \\times 10^{{{e}}}$"


def pm(s):
    s = str(s).strip()
    if "+-" in s:
        a, b = [x.strip() for x in s.split("+-")]
        return f"${a} \\pm {b}$"
    return esc(s)


def esc(s):
    s = (s.replace("Mn K-alpha", "Mn K$\\alpha$").replace("Mn K-beta", "Mn K$\\beta$").replace("F K-alpha", "F K$\\alpha$")
          .replace("C K-alpha", "C K$\\alpha$").replace("5-95%", "5--95%").replace("1000-4000", "1000--4000"))
    return (s.replace("%", "\\%").replace("&", "\\&").replace("_", "\\_").replace("~", "$\\sim$")
             .replace(">=", "$\\geq$").replace("> ", "$>$\\,").replace("< ", "$<$\\,"))


def add(key, cap, colspec, head, body, notes=()):
    T.append((key, cap, colspec, head, body, list(notes)))


def unit_split(q):
    m = re.match(r"(.*?)\s*\[([^\]]+)\](.*)$", q)
    return ((m.group(1) + m.group(3)).strip(), m.group(2)) if m else (q, "")


def unit_tex(u):
    return esc(u).replace("e-/pixel/s", "e$^-$/pixel/sec").replace("e-", "e$^-$")


# ------------------------------------------------------------------ IMX226
r = rd(os.path.join(R226, "Table_IMX226_hot_pixels.csv"))
if r:
    body = [f"{esc(a)} & {esc(d)} & {num(b,0)} & {num(c,2)} & {esc(e)} \\\\" for a, b, c, d, e in r[1:]]
    add("Table_IMX226_01_Hot_pixels", "Hot pixels identified in the IMX226 X-ray exposures.", "llrrl",
        [f" & & \\multicolumn{{2}}{{c}}{{{ME}}} & \\\\", "\\cmidrule(lr){3-4}",
         f"Source & Criterion & Pixels & \\% of sensor & {LU} \\\\"], body)
for k, f, lab in (("02", "single_pixel", "single-pixel (n1)"), ("03", "multi_pixel", "multi-pixel (n1234)")):
    r = rd(os.path.join(R226, f"Table_IMX226_energy_resolution_{f}.csv"))
    if r:
        body = [f"{LINE.get(a,a)} & {num(b,0)} & {num(c,1)} & {pm(d)} & {pm(g)} & {num(e,1)} \\\\" for a, b, c, d, e, g in r[1:]]
        add(f"Table_IMX226_{k}_Energy_resolution_{f}", f"Energy resolution of the IMX226 for {lab} events.", "lrrccr",
            [" & \\multicolumn{2}{c}{Line energy [eV]} & \\multicolumn{3}{c}{FWHM [eV]} \\\\", "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
             f"Emission line & Expected & Measured & {ME} & {LU} & Fano limit \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_total_counts.csv"))
if r:
    body = [f"{LINE.get(a,a)} & {sci(b)} & {sci(c)} & {sci(d)} & {sci(e)} \\\\" for a, b, c, d, e in r[1:]]
    add("Table_IMX226_04_Total_counts", "IMX226 total counts in each emission line [Counts].", "lrrrr",
        [" & \\multicolumn{2}{c}{Single-pixel (n1)} & \\multicolumn{2}{c}{Multi-pixel (n1234)} \\\\", "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
         f"Emission line & {ME} & {LU} & {ME} & {LU} \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_count_rates.csv"))
if r:
    fm = lambda x: f"{float(x):.4f}" if float(x) < 10 else f"{float(x):.2f}"
    body = [f"{LINE.get(a,a)} & {fm(b)} & {fm(c)} & {fm(d)} & {fm(e)} \\\\" for a, b, c, d, e in r[1:]]
    add("Table_IMX226_05_Count_rates", "IMX226 count rates in each emission line [Counts/Sec/mm$^2$].", "lrrrr",
        [" & \\multicolumn{2}{c}{Single-pixel (n1)} & \\multicolumn{2}{c}{Multi-pixel (n1234)} \\\\", "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
         f"Emission line & {ME} & {LU} & {ME} & {LU} \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_readout_noise.csv"))
if r:
    body = [f"{a} & {float(b):.2f} & {num(c,2)} \\\\" for a, b, c in r[1:]]
    add("Table_IMX226_06_Readout_noise", "IMX226 readout noise [electron RMS].", "lrr", [f"Statistic & {ME} & {LU} \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_summary_gain_linearity_dark.csv"))
if r:
    body = []
    for q, me, lu in r[1:]:
        n, u = unit_split(q)
        try:
            v = float(me); mv = f"{v:,.0f}".replace(",", "{,}") if v > 100 else (f"{v:.4f}" if v < 1 else f"{v:.2f}")
        except ValueError:
            mv = esc(me)
        body.append(f"{esc(n)} & {unit_tex(u)} & {mv} & {esc(lu) if lu != '-' else '---'} \\\\")
    add("Table_IMX226_07_Summary_gain_linearity_dark", "Summary of IMX226 gain, linearity and dark-current results.", "llrr",
        [f"Quantity & Unit & {ME} & {LU} \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_energy_resolution_budget.csv"))
if r:
    h = r[0]; ix = {k: h.index(k) for k in h}
    body, prev = [], None
    for row in r[1:]:
        l = row[ix["Line"]]
        if l in ("Si K-alpha", "Si escape (Mn K-alpha)"):
            continue
        if prev and l != prev:
            body.append("\\addlinespace[2pt]")
        lab = LINE.get(l, esc(l)) if l != prev else ""; prev = l
        body.append(f"{lab} & {row[ix['Events']]} & {float(row[ix['Measured_FWHM_eV']]):.1f} & {float(row[ix['Fano_FWHM_eV']]):.1f} & "
                    f"{float(row[ix['Calc_FWHM_Fano_Read_eV']]):.1f} & {float(row[ix['Calc_FWHM_Fano_Read_Gain_eV']]):.1f} & "
                    f"{float(row[ix['Unaccounted_sigma_eV']]):.1f} & {num(row[ix['Lupo2025_unaccounted_eV']],1)} \\\\")
    add("Table_IMX226_08_Energy_resolution_budget", "IMX226 energy-resolution budget.", "llrrrrrr",
        [" & & & \\multicolumn{3}{c}{Calculated FWHM [eV]} & \\multicolumn{2}{c}{Unaccounted $\\sigma$ [eV]} \\\\",
         "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
         f"Line & Events & Measured FWHM [eV] & Fano & + read noise & + gain var. & {ME} & {LU} \\\\"], body)
r = rd(os.path.join(R226, "Table_IMX226_beamline_background.csv"))
if r:
    body = [f"{int(float(a))} & {int(float(b))} & {int(float(c))} & {float(d):.1f} & {float(e):.2f} & {float(g):.3f} \\\\" for a, b, c, d, e, g in r[1:]]
    add("Table_IMX226_09_Beamline_background", "IMX226 event rates on the MIT polarimetry beamline (13 Nov 2024, 23\\,$^\\circ$C, 0.2 s frames).", "rrrrrr",
        [" & & & \\multicolumn{2}{c}{Event rate [Counts/Sec/cm$^2$]} & \\\\", "\\cmidrule(lr){4-5}",
         "Cube & Frames & Warm pixels & All ($>$\\,40 ADU) & $>$\\,300 eV & Residual offset [ADU] \\\\"], body)

# ------------------------------------------------------------------ IMX290
r = rd(os.path.join(R290, "Table_IMX290_ECC_summary.csv"))
if r:
    body = []
    for q, me, bs in r[1:]:
        n, u = unit_split(q)
        def f(x):
            try:
                v = float(x)
            except ValueError:
                return esc(x)
            if math.isnan(v): return "---"
            return f"{v:,.0f}".replace(",", "{,}") if v >= 1000 else (f"{v:.5f}" if v < 1 else f"{v:.1f}")
        body.append(f"{esc(n)} & {unit_tex(u)} & {f(me)} & {f(bs)} \\\\")
    add("Table_IMX290_01_ECC_summary", "IMX290 per-pixel gain calibration (Fe-55, ECC) compared with Schneider et al. (2024).", "llrr",
        [f"Quantity & Unit & {ME} & {BS} \\\\"], body)
r = rd(os.path.join(R290, "Table_IMX290_accuracy_vs_events.csv"))
if r:
    body = [f"{int(float(a)):,} & {float(b):.3f} & {float(c):.1f} & ${float(d):.1f} \\pm {float(e):.1f}$ \\\\".replace(",", "{,}") for a, b, c, d, e in r[1:]]
    add("Table_IMX290_02_Accuracy_vs_events", "IMX290 calibration accuracy versus the number of events per pixel.", "rrrc",
        ["Events per pixel & RMS gain error [\\%] & RMS gain error at 5.9 keV [eV] & Corrected Mn K$\\alpha$ FWHM [eV] \\\\"], body)

# ------------------------------------------------------------------ comparison
def getsum(path, key):
    r = rd(path)
    if not r: return None
    for row in r[1:]:
        if row[0].startswith(key): return row[1]
    return None
s662 = os.path.join(T662, "Table_07_Summary_gain_linearity_dark.csv"); s226 = os.path.join(R226, "Table_IMX226_summary_gain_linearity_dark.csv")
e662 = rd(os.path.join(T662, "Table_02_Energy_resolution_single_pixel.csv")); e226 = rd(os.path.join(R226, "Table_IMX226_energy_resolution_single_pixel.csv"))
rn662 = rd(os.path.join(T662, "Table_06_Readout_noise.csv")); rn226 = rd(os.path.join(R226, "Table_IMX226_readout_noise.csv"))
s290 = rd(os.path.join(R290, "Table_IMX290_ECC_summary.csv"))
def fw(t, line):
    if not t: return "---"
    for row in t[1:]:
        if row[0] == line: return pm(row[3])
    return "---"
def f4(x, d=4):
    try: return f"{float(x):.{d}f}"
    except (TypeError, ValueError): return "---"
g290 = s290[1][1] if s290 else None; fw290 = (s290[4][1], s290[5][1]) if s290 else ("", "")
rows = [("Sensor family", "", "Sony STARVIS", "Sony STARVIS", "Sony STARVIS"),
        ("Pixel size", "$\\mu$m", "2.9", "1.85", "2.9"),
        ("Conversion gain", "ADU/eV", f4(getsum(s662, "Conversion gain from Mn")), f4(getsum(s226, "Conversion gain from F")), f4(g290, 4)),
        ("Readout noise, mode", "e$^-$ rms", f4(rn662[3][1], 2) if rn662 else "---", f4(rn226[3][1], 2) if rn226 else "---", "1.9$^{a}$"),
        ("C K$\\alpha$ FWHM, n1", "eV", fw(e662, "C K-alpha"), fw(e226, "C K-alpha"), "---"),
        ("F K$\\alpha$ FWHM, n1", "eV", fw(e662, "F K-alpha"), fw(e226, "F K-alpha"), "---"),
        ("Mn K$\\alpha$ FWHM, n1 (single gain)", "eV", fw(e662, "Mn K-alpha"), fw(e226, "Mn K-alpha"), f4(fw290[0], 1)),
        ("Mn K$\\alpha$ FWHM, per-pixel gain", "eV", "---", "---", f4(fw290[1], 1)),
        ("Pixel-to-pixel gain variation", "\\%", f4(getsum(s662, "Pixel-to-pixel"), 2), f4(getsum(s226, "Pixel-to-pixel"), 2), f4(s290[2][1], 2) if s290 else "---"),
        ("Dark current, lab method", "e$^-$/pixel/sec", f4(getsum(s662, "Dark current, lab method, full"), 3), f4(getsum(s226, "Dark current, lab method, full"), 3), "---")]
body = [" & ".join(x) + " \\\\" for x in rows]
add("Table_Comparison_IMX662_IMX226_IMX290", "Comparison of the three Sony STARVIS sensors characterized with the lab's X-ray setup (all values Carissma 2026).",
    "llccc", ["Quantity & Unit & IMX662 & IMX226 & IMX290 \\\\"], body, ["$^{a}$\\,Schneider et al. (2024)."])

# ------------------------------------------------------------------ writing
PRE = r"""\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{booktabs,threeparttable,array}
\usepackage[font=small,labelfont=bf,skip=4pt,justification=raggedright,singlelinecheck=false]{caption}
\renewcommand{\arraystretch}{1.15}
\setlength{\tabcolsep}{7pt}
"""
def table_tex(num, cap, colspec, head, body, notes, standalone):
    s = [f"\\setcounter{{table}}{{{num-1}}}" if standalone else "\\begin{table}[htbp]\\centering", "\\begin{threeparttable}",
         ("\\captionof{table}{" if standalone else "\\caption{") + cap + "}", f"\\begin{{tabular}}{{@{{}}{colspec}@{{}}}}", "\\toprule"]
    s += head + ["\\midrule"] + body + ["\\bottomrule", "\\end{tabular}"]
    if notes:
        s += ["\\begin{tablenotes}[flushleft]\\footnotesize"] + ["\\item " + n for n in notes] + ["\\end{tablenotes}"]
    s += ["\\end{threeparttable}"] + ([] if standalone else ["\\end{table}"])
    return "\n".join(s)

build = os.path.join(OUT, "_build"); os.makedirs(build, exist_ok=True)
def run(tex):
    p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex], cwd=build, capture_output=True, text=True)
    if p.returncode:
        print(p.stdout[-2000:]); raise SystemExit("LaTeX failed: " + tex)
for i, (key, cap, cs, head, body, notes) in enumerate(T, 1):
    doc = "\\documentclass[varwidth=26cm,border=10pt]{standalone}\n" + PRE + "\\begin{document}\n" + table_tex(i, cap, cs, head, body, notes, True) + "\n\\end{document}\n"
    open(os.path.join(build, key + ".tex"), "w").write(doc); run(key + ".tex")
    for ext in ("pdf", "tex"): shutil.copy(os.path.join(build, f"{key}.{ext}"), OUT)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile", os.path.join(build, key + ".pdf"), os.path.join(OUT, key)], check=True)
    print("built", key)
shutil.rmtree(build)
