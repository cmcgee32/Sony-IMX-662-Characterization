#!/usr/bin/env python3
r"""Build academic (booktabs) summary tables for the IMX662 characterization.

Reads the CSVs written by the MATLAB pipeline (imx662_tables.m) and writes, for
each table: a standalone LaTeX source (.tex), a vector PDF, a 300-dpi PNG and a
clean CSV. Also writes All_summary_tables.tex/.pdf (all tables in one document,
ready to paste into Overleaf).

Style: three horizontal rules only (\toprule, \midrule, \bottomrule), no
vertical lines, \cmidrule for grouped headers, caption above, notes below.
"""
import csv, math, os, re, shutil, subprocess, sys

SRC = sys.argv[1]          # folder with the pipeline CSVs
OUT = sys.argv[2]          # output folder
os.makedirs(OUT, exist_ok=True)
ME, LU = "Carissma (2026)", "Lupo (2025)"


def rd(name):
    with open(os.path.join(SRC, name), newline="") as f:
        return list(csv.reader(f))[1:]


def pm(s):
    """'43.6 +- 0.9' -> '$43.6 \\pm 0.9$'"""
    s = s.strip()
    if "+-" in s:
        a, b = [x.strip() for x in s.split("+-")]
        return f"${a} \\pm {b}$"
    return s


def f(x, d=1):
    try:
        v = float(x)
    except ValueError:
        return x
    if math.isnan(v):
        return "---"
    return f"{v:,.{d}f}".replace(",", "{,}")


def sci(x, d=2):
    v = float(x)
    e = int(math.floor(math.log10(abs(v))))
    return f"${v/10**e:.{d}f} \\times 10^{{{e}}}$"


def esc(s):
    s = (s.replace("Mn K-alpha", "Mn K$\\alpha$").replace("5-95%", "5--95%")
          .replace("2000-4000", "2000--4000"))
    if s.isdigit() and len(s) > 3:
        s = f"{int(s):,}".replace(",", "{,}")
    return (s.replace("%", "\\%").replace("&", "\\&").replace("_", "\\_")
             .replace("~", "$\\sim$").replace(">=", "$\\geq$")
             .replace("> ", "$>$\\,").replace("< ", "$<$\\,"))


LINE = {"C K-alpha": "C K$\\alpha$", "F K-alpha": "F K$\\alpha$",
        "Mn K-alpha": "Mn K$\\alpha$", "Mn K-beta": "Mn K$\\beta$"}

TABLES = []   # (key, number, caption, colspec, header_lines, body_lines, notes)


def add(key, cap, colspec, head, body, notes):
    TABLES.append((key, len(TABLES) + 1, cap, colspec, head, body, notes))


# 1 ---------------------------------------------------------------- hot pixels
r = rd("Table_01_hot_pixels.csv")
body = []
for src, n, pct, thr, lup in r:
    m = re.match(r"(\d+) \(([\d.]+)%\)", lup)
    body.append(f"{src.replace('-', '-')} & {esc(thr)} & {f(n,0)} & {float(pct):.2f} & "
                f"{f(m.group(1),0)} & {m.group(2)} \\\\")
add("Table_01_Hot_pixels",
    "Hot pixels identified in the IMX662 X-ray exposures.",
    "llrrrr",
    [" & & \\multicolumn{2}{c}{" + ME + "} & \\multicolumn{2}{c}{" + LU + "} \\\\",
     "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
     "Source & Hot-pixel criterion & Pixels & \\% of sensor & Pixels & \\% of sensor \\\\"],
    body,
    ["A hot pixel is a pixel that records far more events than an X-ray source "
     "can produce, usually because of a defect or high leakage current; these "
     "pixels are removed before building spectra. Sensor: $1920\\times1080$ "
     "= 2{,}073{,}600 pixels. Fe-55: Gaussian centre + $3\\sigma$ of the "
     "events-per-pixel histogram. Po-210: $\\geq 10$ events per pixel."])


# 2,3 --------------------------------------------------------- energy resolution
def eres(fname, key, label, evdesc):
    r = rd(fname)
    body = [f"{LINE[l]} & {f(e,0)} & {f(c,1)} & {pm(m)} & {pm(lu)} & {f(fa,1)} \\\\"
            for l, e, c, m, fa, lu in r]
    add(key,
        f"Energy resolution of the IMX662 for {label} events.",
        "lrrccr",
        [" & \\multicolumn{2}{c}{Line energy [eV]} & \\multicolumn{3}{c}{FWHM [eV]} \\\\",
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
         f"Emission line & Expected & Measured & {ME} & {LU} & Fano limit \\\\"],
        body,
        [f"{evdesc} Energy resolution is the full width at half maximum (FWHM) "
         "of a Gaussian fit to each emission line: smaller values mean the sensor "
         "separates nearby X-ray energies better. The Fano limit, "
         "$2.3548\\sqrt{F\\,w\\,E}$ with $F=0.115$ and $w=3.65$~eV, is the best "
         "resolution physically possible in silicon. Uncertainties are $1\\sigma$ "
         "fit errors. C and F lines from Po-210 fluorescence; Mn lines from Fe-55."])


eres("Table_03_energy_resolution_single_pixel.csv", "Table_02_Energy_resolution_single_pixel",
     "single-pixel (n1)", "Single-pixel (n1) events deposit all charge in one pixel.")
eres("Table_04_energy_resolution_multi_pixel.csv", "Table_03_Energy_resolution_multi_pixel",
     "multi-pixel (n1234)", "Multi-pixel (n1234) events spread charge over up to four pixels.")


# 4,5 ------------------------------------------------------ counts, count rates
def grouped(fname, key, cap, fmt, unit, notes):
    r = rd(fname)
    body = [f"{LINE[l]} & {fmt(a)} & {fmt(b)} & {fmt(c)} & {fmt(d)} \\\\" for l, a, b, c, d in r]
    add(key, cap, "lrrrr",
        [" & \\multicolumn{2}{c}{Single-pixel (n1)} & \\multicolumn{2}{c}{Multi-pixel (n1234)} \\\\",
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
         f"Emission line & {ME} & {LU} & {ME} & {LU} \\\\"],
        body, notes)


grouped("Table_05_06_total_counts.csv", "Table_04_Total_counts",
        "Total counts in each emission line (Gaussian area $a\\sigma\\sqrt{2\\pi}$).",
        sci, "",
        ["Counts are the area under the fitted Gaussian of each line. Totals scale "
         "with the length of the exposure, so they are not expected to match "
         "between runs of different duration; Table~5 normalises by time and area."])
grouped("Table_07_08_count_rates.csv", "Table_05_Count_rates",
        "Count rates in each emission line [Counts/Sec/mm$^2$].",
        lambda x: (f"{float(x):.4f}" if float(x) < 1 else f"{float(x):.2f}"), "",
        ["Rate = counts in the line divided by exposure time and by the sensor "
         "area (17.44~mm$^2$), so runs of different length can be compared "
         "directly. Agreement with Lupo (2025) to within a few per cent "
         "shows that the source and detector response are stable."])

# 6 -------------------------------------------------------------- readout noise
r = rd("Table_09_readout_noise.csv")
body = [f"{p} & {float(a):.2f} & {f(b,2)} \\\\" for p, a, b in r]
add("Table_06_Readout_noise", "Readout noise of the IMX662 [electron RMS].", "lrr",
    [f"Statistic & {ME} & {LU} \\\\"], body,
    ["Readout noise is the random electronic noise added each time a pixel is "
     "read, expressed in electrons. Per-pixel noise from the bias maps of all "
     "92 Fe-55 and Po-210 runs, converted from ADU with the measured gain "
     "(0.2652~ADU\\,eV$^{-1}$) and $w=3.65$~eV per electron. "
     "Average, RMS and mode are taken over all pixels."])

# 7 ------------------------------------------------------------------ summary
r = rd("Table_summary_gain_linearity_dark.csv")
groups = [("Gain", 0, 6), ("Linearity", 6, 11), ("Dark current", 11, 16)]


def split_unit(q):
    m = re.match(r"(.*?)\s*\[([^\]]+)\](.*)$", q)
    if m:
        return (m.group(1) + m.group(3)).strip(), m.group(2)
    return q, ""


def unit_tex(u):
    return (u.replace("e-/pixel/s", "e$^-$/pixel/sec")
             .replace("%", "\\%").replace("eV", "eV"))


body = []
for gname, a, b in groups:
    body.append(f"\\addlinespace[3pt]\\multicolumn{{4}}{{l}}{{\\textit{{{gname}}}}} \\\\")
    for q, me, lu in r[a:b]:
        name, u = split_unit(q)
        extra = "(extra)" in name
        name = name.replace("(extra)", "").strip()
        try:
            v = float(me)
            mv = f"{v:,.0f}".replace(",", "{,}") if v > 100 else (f"{v:.4f}" if v < 1 else f"{v:.2f}")
        except ValueError:
            mv = me
        lu_t = esc(lu) if lu != "-" else "---"
        body.append(f"\\quad {esc(name)}{'$^{a}$' if extra else ''} & {unit_tex(u)} & {mv} & {lu_t} \\\\")
add("Table_07_Summary_gain_linearity_dark",
    "Summary of IMX662 gain, linearity and dark-current results.", "llrr",
    [f"Quantity & Unit & {ME} & {LU} \\\\"], body,
    ["$^{a}$\\,Additional measurement not reported by Lupo (2025). Conversion gain "
     "is the number of ADU (digital units) produced per eV of X-ray energy. "
     "Linearity error is the largest deviation of the data from a straight-line fit "
     "as a percentage of full range. The lab-method dark current is below what "
     "these 4.5~s exposures can measure; the shot-noise method gives an upper estimate."])

# 8 ------------------------------------------------------- resolution budget
r = rd("Table_10_15_energy_resolution_budget.csv")
body, prev = [], None
for l, ev, meas, fano, fr, frg, un, lu in r:
    lab = LINE[l] if l != prev else ""
    if prev and l != prev:
        body.append("\\addlinespace[2pt]")
    prev = l
    body.append(f"{lab} & {ev} & {float(meas):.1f} & {float(fano):.1f} & {float(fr):.1f} & "
                f"{float(frg):.1f} & {float(un):.1f} & {float(lu):.1f} \\\\")
add("Table_08_Energy_resolution_budget",
    "Energy-resolution budget: measured FWHM against the expected noise contributions.",
    "llrrrrrr",
    [" & & & \\multicolumn{3}{c}{Calculated FWHM [eV]} & \\multicolumn{2}{c}{Unaccounted $\\sigma$ [eV]} \\\\",
     "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
     "Line & Events & Measured FWHM [eV] & Fano & + read noise & + gain var. & " + ME + " & " + LU + " \\\\"],
    body,
    ["Each calculated column adds one noise source in quadrature: Fano statistics, "
     "then readout noise (3.66~e$^-$), then pixel-to-pixel gain variation "
     "(applied at the Mn lines only, where it was measured). The unaccounted "
     "$\\sigma$ is the extra broadening needed to reach the measured width "
     "(e.g.\\ charge sharing, residual calibration). n1 = single-pixel, "
     "n1234 = up to four-pixel events."])

# 9 -------------------------------------------------------------- dark signal
r = rd("Table_dark_signal_vs_exposure.csv")
body = [f"{float(t):g} & {float(a):.2f} & {float(s):.2f} & {float(tm):.2f} & {float(c):.2f} & {float(e):.2f} \\\\"
        for t, a, s, tm, c, e in r]
body[0] = body[0].replace("1.5e-05", "$1.5\\times10^{-5}$")
add("Table_09_Dark_signal_vs_exposure",
    "Dark signal of the IMX662 versus exposure time [e$^-$/pixel].",
    "rrrrrr",
    [" & \\multicolumn{3}{c}{Full detector} & \\multicolumn{2}{c}{Region} \\\\",
     "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
     "Exposure [s] & Signal & Spatial error & Temporal error & Centre & Edges \\\\"],
    body,
    ["Dark signal is the charge a pixel collects with the shutter closed, from "
     "thermal leakage current. Each value is the clipped mean of the 25-frame "
     "median minus the 15~$\\mu$s bias frame. Spatial error = RMS across pixels; "
     "temporal error = mean frame-to-frame scatter$/\\sqrt{25}$. "
     "Edges = outer 375 pixel columns on each side; centre = the rest. "
     "All values are consistent with zero within the errors."])

# ------------------------------------------------------------------ writing
PRE = r"""\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{booktabs,threeparttable,array}
\usepackage[font=small,labelfont=bf,skip=4pt,justification=raggedright,singlelinecheck=false]{caption}
\renewcommand{\arraystretch}{1.15}
\setlength{\tabcolsep}{7pt}
"""


def table_tex(num, cap, colspec, head, body, notes, standalone):
    notes = [n for n in notes if n.startswith("$^{a}$")]
    notes = [n.split(" Conversion gain")[0] for n in notes]
    s = []
    if standalone:
        s.append(f"\\setcounter{{table}}{{{num-1}}}")
    else:
        s.append("\\begin{table}[htbp]\\centering")
    s.append("\\begin{threeparttable}")
    s.append(("\\captionof{table}{" if standalone else "\\caption{") + cap + "}")
    s.append(f"\\begin{{tabular}}{{@{{}}{colspec}@{{}}}}")
    s.append("\\toprule")
    s += head
    s.append("\\midrule")
    s += body
    s.append("\\bottomrule")
    s.append("\\end{tabular}")
    if notes:
        s.append("\\begin{tablenotes}[flushleft]\\footnotesize")
        for n in notes:
            s.append("\\item " + n)
        s.append("\\end{tablenotes}")
    s.append("\\end{threeparttable}")
    if not standalone:
        s.append("\\end{table}")
    return "\n".join(s)


def run(tex, d):
    p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex],
                       cwd=d, capture_output=True, text=True)
    if p.returncode:
        print(p.stdout[-2500:])
        raise SystemExit(f"LaTeX failed: {tex}")


build = os.path.join(OUT, "_build")
os.makedirs(build, exist_ok=True)
for key, num, cap, colspec, head, body, notes in TABLES:
    doc = ("\\documentclass[varwidth=26cm,border=10pt]{standalone}\n" + PRE +
           "\\begin{document}\n" + table_tex(num, cap, colspec, head, body, notes, True) +
           "\n\\end{document}\n")
    with open(os.path.join(build, key + ".tex"), "w") as fh:
        fh.write(doc)
    run(key + ".tex", build)
    shutil.copy(os.path.join(build, key + ".pdf"), OUT)
    shutil.copy(os.path.join(build, key + ".tex"), OUT)
    subprocess.run(["pdftoppm", "-png", "-r", "300", "-singlefile",
                    os.path.join(build, key + ".pdf"), os.path.join(OUT, key)], check=True)
    print("built", key)

# all tables in one article (Overleaf-ready)
alldoc = ("\\documentclass[11pt]{article}\n\\usepackage[margin=2cm,landscape]{geometry}\n" + PRE +
          "\\begin{document}\n\\section*{IMX662 characterization: summary tables}\n"
          "\\noindent Carissma (2026), compared with Lupo (2025, MIT S.M.\\ thesis).\n\n" +
          "\n\n\\clearpage\n".join(table_tex(n, c, cs, h, b, no, False)
                                    for _, n, c, cs, h, b, no in TABLES) +
          "\n\\end{document}\n")
with open(os.path.join(build, "All_summary_tables.tex"), "w") as fh:
    fh.write(alldoc)
run("All_summary_tables.tex", build)
for ext in ("tex", "pdf"):
    shutil.copy(os.path.join(build, "All_summary_tables." + ext), OUT)

# clean CSVs with the renamed column
mapping = {"Table_01_hot_pixels.csv": "Table_01_Hot_pixels.csv",
           "Table_03_energy_resolution_single_pixel.csv": "Table_02_Energy_resolution_single_pixel.csv",
           "Table_04_energy_resolution_multi_pixel.csv": "Table_03_Energy_resolution_multi_pixel.csv",
           "Table_05_06_total_counts.csv": "Table_04_Total_counts.csv",
           "Table_07_08_count_rates.csv": "Table_05_Count_rates.csv",
           "Table_09_readout_noise.csv": "Table_06_Readout_noise.csv",
           "Table_summary_gain_linearity_dark.csv": "Table_07_Summary_gain_linearity_dark.csv",
           "Table_10_15_energy_resolution_budget.csv": "Table_08_Energy_resolution_budget.csv",
           "Table_dark_signal_vs_exposure.csv": "Table_09_Dark_signal_vs_exposure.csv"}
for a, b in mapping.items():
    with open(os.path.join(SRC, a), newline="") as fi:
        rows = list(csv.reader(fi))
    H = {"This analysis": ME,
         "Hot pixels": f"Hot pixels, {ME}",
         "Percent of sensor [%]": f"Percent of sensor, {ME} [%]",
         "Energy resolution FWHM [eV]": f"FWHM, {ME} [eV]",
         "Total counts, single-pixel": f"Total counts single-pixel, {ME}",
         "Total counts, multi-pixel": f"Total counts multi-pixel, {ME}",
         "Count rate, single [counts/s/mm^2]": f"Count rate single-pixel, {ME} [counts/s/mm^2]",
         "Count rate, multi [counts/s/mm^2]": f"Count rate multi-pixel, {ME} [counts/s/mm^2]",
         "Readout noise [e- RMS]": f"{ME} [e- RMS]",
         "Unaccounted sigma [eV]": f"Unaccounted sigma, {ME} [eV]"}
    rows[0] = [H.get(c, c) for c in rows[0]]
    with open(os.path.join(OUT, b), "w", newline="") as fo:
        csv.writer(fo).writerows(rows)
shutil.rmtree(build)
print("done")
