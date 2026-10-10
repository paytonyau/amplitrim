"""
Outputs: REPORT.md (the one file to read), summary.png, summary.json,
plus the short summary printed to the screen.
"""

import json
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import __version__
from .config import ee_pair
from .pipeline import RULE_CAP, RULE_Q, Result


def _fmt(x) -> str:
    return f"{x:g}"


def _ee_f(max_ee) -> str:
    return _fmt(ee_pair(max_ee)[0])


def _ee_r(max_ee) -> str:
    return _fmt(ee_pair(max_ee)[1])


def _ee_text(max_ee) -> str:
    """maxEE as written in the report: one number when both reads share it, otherwise 'F,R'."""
    f, r = ee_pair(max_ee)
    return _fmt(f) if f == r else f"{_fmt(f)},{_fmt(r)}"


LOW_PASS = 0.85      # below this share of pairs passing the filter, the report adds a note that a low share is not a problem in itself

_AUTHOR_NOTE = ("DADA2's author calls maxEE = 2 a good general purpose value and says that relaxing it somewhat is also OK, because DADA2's error "
                "model uses the quality scores and so tolerates some noisier reads; but pushing through more lower-quality reads may not improve the "
                "relative abundances and could make them worse (DADA2 issue #2136, about PacBio reads). The table shows what a looser limit costs: "
                "the expected errors in the kept reads.")


def _maxee_note(res: Result, short: bool = False) -> str:
    """Why this maxEE: DADA2's default of 2, or the limit that covers three quarters of the reads at the recommended cut."""
    r, rule = res.recommendation, res.rule
    ee = _ee_text(res.cfg.max_ee)
    if rule is None:                                                          # ITS, a maxEE given in the settings, or no cut to measure at
        if res.marker == "its":
            return "maxEE is 2, DADA2's default. AmpliTrim does not adjust it for ITS: reads are not cut at a common position, so there is no cut to measure them at (each read is only trimmed where its quality first falls to 2 or below), and the rule has not been tested on reads of different lengths."
        if ee_pair(res.cfg.max_ee) == (2.0, 2.0):
            return f"maxEE is 2, DADA2's default: {r.ee_pass_rate:.0%} of the read pairs pass it at these settings."
        return f"maxEE is {ee}: {r.ee_pass_rate:.0%} of the read pairs pass it at these settings."
    pf, pr = rule["percentile"]
    wanted = f"{rule['max_ee'][0]},{rule['max_ee'][1]}"
    if res.relaxed:
        had = "none of the read pairs pass" if not res.pass_rate_default else f"only {res.pass_rate_default:.0%} of the read pairs pass"
        if short:
            return (f"maxEE is {ee}, not 2: {had} at maxEE 2; each limit covers three quarters of the reads in its direction (never above {RULE_CAP}).")
        cap = f", and the limit was capped at {RULE_CAP}" if rule["capped"] else ""
        return (f"maxEE is {ee}, not DADA2's default of 2: at the recommended cut the {RULE_Q * 100:.0f}th percentile of the reads' expected errors is {pf:.1f} forward and {pr:.1f} "
                f"reverse, so each read's limit is set to cover three quarters of the reads in its direction (never below 2, never above {RULE_CAP}{cap}). {had} at maxEE 2 and "
                f"{r.ee_pass_rate:.0%} pass now; the reads this admits are noisier (see the table below). DADA2's author says that relaxing maxEE somewhat is acceptable because "
                f"DADA2's error model uses the quality scores; the {RULE_Q * 100:.0f}th percentile and the cap of {RULE_CAP} are AmpliTrim's own and were chosen after seeing the validation data.")
    if tuple(rule["max_ee"]) == (2, 2):
        return (f"maxEE is 2, DADA2's default: at the recommended cut the {RULE_Q * 100:.0f}th percentile of the reads' expected errors is {pf:.1f} forward and {pr:.1f} reverse, both under 2, "
                f"so the default already covers three quarters of the reads in each direction; {r.ee_pass_rate:.0%} of the read pairs pass.")
    return (f"maxEE stays 2: the rule would set {wanted} (the {RULE_Q * 100:.0f}th percentile of the reads' expected errors is {pf:.1f} forward and {pr:.1f} reverse), "
            f"but nothing could merge your longest amplicons at that setting either.")


def _load_text(x) -> str:
    """Expected errors in the average kept read, forward / reverse."""
    return "n/a" if x.get("ee_f") is None or x["ee_f"] != x["ee_f"] else f"{x['ee_f']:.2f} / {x['ee_r']:.2f}"


def _added_text(x) -> str:
    """Expected errors in the reads this maxEE adds beyond maxEE 2, forward / reverse (a dash when it adds none)."""
    a = x.get("add_f")
    return "-" if a is None or a != a else f"{x['add_f']:.2f} / {x['add_r']:.2f}"


def ee_advice(rows: list, current, relaxed: bool = False) -> List[str]:
    """What to say about maxEE next to the table. It states DADA2's guidance and never asks you to change anything."""
    cur = ee_pair(current)
    now = next((x for x in rows if ee_pair(x["ee"]) == cur), None)
    if now is None:
        return []
    lines = [f"At maxEE {_ee_text(cur)}, **{now['pass_rate']:.1%}** of the read pairs pass DADA2's quality filter at these truncation lengths.",
             _AUTHOR_NOTE]
    if now["pass_rate"] < LOW_PASS and not relaxed:
        lines.append("A low share is not a problem in itself: DADA2's author says that the number of reads is not the aim, and that keeping about 40% "
                     "after filtering is fine (DADA2 issue #2136, about PacBio reads).")
    return lines


def _suggestion(res: Result):
    from .ranking import quality_first_suggestion
    return quality_first_suggestion(res.ee_search, res.min_retention) if res.ee_search else (None, False)


def _ee_search_section(res: Result) -> List[str]:
    from .ranking import dominators
    rows, floor = res.ee_search, res.min_retention
    out = ["## maxEE: retention against read quality", "",
           "Each maxEE below was given its own best truncation (the same search as for your setting), so the comparison is not tied to "
           "one pair of truncation lengths:", "",
           "| maxEE | truncLenF | truncLenR | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | "
           "expected errors in the reads added beyond maxEE 2 (F / R) | beaten by |", "|---|---|---|---|---|---|---|---|"]
    beaten = dominators(rows)
    for x, better in zip(rows, beaten):
        label = _ee_text(x["ee"]) + (" (yours)" if x["current"] else "")
        label = f"**{label}**" if x["current"] else label
        by = _ee_text(max(better, key=lambda o: o["pass_rate"])["ee"]) if better else "-"
        out.append(f"| {label} | no setting | | | | | | |" if x["trim_f"] is None else
                   f"| {label} | {x['trim_f']} | {x['trim_r']} | {x['pass_rate']:.1%} | {x['yield_est']:.1%} | {_load_text(x)} | {_added_text(x)} | {by} |")
    now = next((x for x in rows if x["current"]), None)
    out.append("")
    if now is None or now["yield_est"] is None:
        return out + ["No maxEE in the table had a setting that could merge your longest amplicons.", ""]
    mine = beaten[rows.index(now)]
    if mine:
        top = max(mine, key=lambda o: o["pass_rate"])
        out += [f"Your setting is **beaten** by maxEE {_ee_text(top['ee'])} with truncLenF {top['trim_f']} and truncLenR {top['trim_r']}: it keeps more "
                f"pairs ({top['pass_rate']:.1%} against {now['pass_rate']:.1%}) with fewer expected errors per kept pair "
                f"({top['ee_f'] + top['ee_r']:.2f} against {now['ee_f'] + now['ee_r']:.2f}). \"Beaten by\" means another row is better on both counts; "
                f"it needs no weighting, but it looks only at retention and read noise, not at accuracy.", ""]
    else:
        out += ["No other maxEE in the table keeps more pairs with fewer expected errors per kept pair, so yours is on the trade-off "
                "frontier between retention and read noise.", ""]
    pick, reached = _suggestion(res)
    base = next((x for x in rows if ee_pair(x["ee"]) == (2.0, 2.0)), None)
    if pick is None:
        out += ["Quality first: no setting at or above DADA2's default maxEE of 2 could merge your longest amplicons.", ""]
    elif reached and ee_pair(pick["ee"]) == (2.0, 2.0):
        out += [f"**Quality first: keep maxEE 2** (DADA2's default). With its own best truncation (truncLenF {pick['trim_f']}, truncLenR {pick['trim_r']}) "
                f"it keeps an estimated {pick['yield_est']:.1%} of the pairs, above the {floor:.0%} floor (`--min-retention`), with the fewest "
                f"expected errors per kept pair of the settings that are not beaten.", ""]
    elif reached:
        had = f"At maxEE 2 only {base['yield_est']:.1%} of the pairs would be usable, below the {floor:.0%} floor (`--min-retention`). " if base and base["yield_est"] is not None else ""
        out += [f"**Quality first: {had}The least noisy setting that reaches the floor is maxEE {_ee_text(pick['ee'])}** with truncLenF {pick['trim_f']} and "
                f"truncLenR {pick['trim_r']}: an estimated {pick['yield_est']:.1%} of the pairs, with {pick['ee_f'] + pick['ee_r']:.2f} expected errors per "
                f"kept pair. Run again with `--max-ee {_fmt(ee_pair(pick['ee'])[0])},{_fmt(ee_pair(pick['ee'])[1])}` to use it.", ""]
    else:
        out += [f"Quality first: no setting at or above maxEE 2 keeps {floor:.0%} of the pairs (`--min-retention`); the most any keeps is "
                f"{pick['yield_est']:.1%}, at maxEE {_ee_text(pick['ee'])} with truncLenF {pick['trim_f']} and truncLenR {pick['trim_r']}.", ""]
    out += [_AUTHOR_NOTE, "",
            "The suggestion weighs how many pairs are kept only against the noise of the reads; it never changes your maxEE, and AmpliTrim does not "
            "measure what a looser maxEE does to the accuracy of your results. If you have a mock community, run DADA2 at the settings in this table and compare the variants it recovers.", ""]
    return out


def _coverage_section(res: Result) -> List[str]:
    rows = res.coverage_search
    if not rows or res.recommendation is None:
        return []
    cur = next((x for x in rows if x["current"]), None)
    out = ["## Coverage against quality", "",
           f"Your settings can merge amplicons up to {res.recommendation.max_mergeable_insert} nt. The longest amplicons come from taxa with the "
           f"longest 16S regions, and merging them needs long reads, but quality falls toward the read ends, so covering fewer of them "
           f"lets the reads be shorter and cleaner. Each row below has its own best truncation:", "",
           "| share of the longest amplicons left out | longest amplicon covered | total length | truncLenF | truncLenR | "
           "pairs passing the filter | estimated yield | expected errors per kept read (F / R) |", "|---|---|---|---|---|---|---|---|"]
    for x in rows:
        label = f"{x['level']:.1%}" + (" (yours)" if x["current"] else "")
        label = f"**{label}**" if x["current"] else label
        out.append(f"| {label} | no setting | | | | | | |" if x["total"] is None else
                   f"| {label} | {x['longest']} nt | {x['total']} | {x['trim_f']} | {x['trim_r']} | {x['pass_rate']:.1%} | "
                   f"{x['yield_est']:.1%} | {_load_text(x)} |")
    out.append("")
    later = [x for x in rows if x["total"] is not None and cur and cur["total"] is not None and x["level"] > cur["level"]]
    best = max(later, key=lambda x: x["yield_est"], default=None)
    if best and cur["yield_est"] is not None and best["yield_est"] - cur["yield_est"] >= 0.03:
        out.append(f"Leaving out the longest {best['level']:.0%} of amplicons (those over {best['longest']} nt) would allow a total of {best['total']} nt "
                   f"instead of {cur['total']}, and raise the estimated yield from {cur['yield_est']:.1%} to {best['yield_est']:.1%}. Those amplicons "
                   f"belong to the taxa with the longest 16S regions, which the shorter setting cannot merge at all. AmpliTrim cannot tell which taxa "
                   f"they are, so the trade is yours. To base the recommendation on another row, run again with `--allow-unmergeable "
                   f"{best['level']:g}`.")
    elif cur and cur["yield_est"] is not None:
        out.append("Leaving out more of the longest amplicons would barely improve the yield here, so the coverage rule is not what limits it.")
    return out + [""]


def _ee_section(res: Result) -> List[str]:
    if res.ee_search:
        return _ee_search_section(res)
    if not res.ee_table:
        return []
    r = res.recommendation
    out = ["## What maxEE keeps", "",
           f"At the recommended truncation ({r.trim_f} + {r.trim_r}), the share of read pairs that pass DADA2's quality filter, and that "
           f"then pass filtering and merge, for other maxEE values (forward, reverse):", "",
           "| maxEE | pairs passing the filter | estimated yield | expected errors per kept read (F / R) | "
           "expected errors in the reads added beyond maxEE 2 (F / R) |", "|---|---|---|---|---|"]
    for x in res.ee_table:
        label = _ee_text(x["ee"]) + (" (used)" if x["current"] else "")
        label = f"**{label}**" if x["current"] else label
        out.append(f"| {label} | {x['pass_rate']:.1%} | {x['yield_est']:.1%} | {_load_text(x)} | {_added_text(x)} |")
    out += ["", *ee_advice(res.ee_table, res.cfg.max_ee, res.relaxed), ""]
    return out


def _settings(res: Result) -> dict:
    r = res.recommendation
    ee_f, ee_r = ee_pair(res.cfg.max_ee)                      # maxEE is one number when both reads share it and a pair otherwise; maxEE_f and maxEE_r are always numbers
    out = {"truncLenF": r.trim_f, "truncLenR": r.trim_r, "maxEE": res.cfg.max_ee, "maxEE_f": ee_f, "maxEE_r": ee_r,
           "minOverlap": res.cfg.dada2_min_overlap, "maxMismatch": res.cfg.dada2_max_mismatch}
    if res.marker == "its":
        out["minLen"] = res.min_len
    return out


def screen_summary(res: Result, report_path: Path) -> List[str]:
    """The result as a few plain lines for the terminal."""
    r = res.recommendation
    if r is None:
        lines = ["", "No recommendation: nothing tested could merge your longest amplicons."]
    elif res.marker == "its":
        s = _settings(res)
        lines = ["", "Recommended DADA2 settings",
                 f"  no fixed truncation (truncLen = 0)   maxEE = {_ee_text(s['maxEE'])}   minLen = {s['minLen']}   "
                 f"minOverlap = {s['minOverlap']}", "",
                 f"  {r.yield_est:.1%} of the {res.sub.n_pairs:,} read pairs analysed should pass filtering and merge "
                 f"({r.ee_pass_rate:.1%} pass the filter; {r.merge_rate_tolerant:.1%} of those merge)."]
    else:
        s = _settings(res)
        lines = ["", "Recommended DADA2 settings",
                 f"  truncLenF = {s['truncLenF']}   truncLenR = {s['truncLenR']}   "
                 f"maxEE = {_ee_text(s['maxEE'])}   minOverlap = {s['minOverlap']}", "",
                 f"  {r.yield_est:.1%} of the {res.sub.n_pairs:,} read pairs analysed should pass filtering "
                 f"and merge.",
                 f"  Covers amplicons up to {r.max_mergeable_insert} nt (your longest: {res.longest_insert} nt)."]
    if r is not None and res.relaxed:
        lines += ["", "  " + _maxee_note(res, short=True)]
    if res.warnings:
        lines += [""] + [f"Warning: {w}" for w in res.warnings]
    lines += ["", f"Report: {report_path}"]
    return lines


def _figure_its(res: Result, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    if res.inserts:
        lo, hi = res.inserts[0], res.inserts[-1]
        bins = min(60, hi - lo + 1)
        ax.hist(res.inserts, bins=bins, color="0.55")
    ax.set_title("Amplicon lengths in your reads")
    ax.set_xlabel("Amplicon length (nt), pairs that merge at full read length")
    ax.set_ylabel("Read pairs")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _figure(res: Result, path: Path) -> None:
    if res.marker == "its":
        return _figure_its(res, path)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))
    r = res.recommendation

    if res.inserts:
        lo, hi = res.inserts[0], res.inserts[-1]
        counts = [0] * (hi - lo + 1)
        for v in res.inserts:
            counts[v - lo] += 1
        ax1.bar(range(lo, hi + 1), counts, color="0.55", width=0.85)
    ax1.axvline(res.longest_insert, color="tab:blue", ls="--",
                label=f"longest amplicon: {res.longest_insert} nt")
    if r is not None:
        ax1.axvline(r.max_mergeable_insert, color="tab:red",
                    label=f"longest these settings can merge: {r.max_mergeable_insert} nt")
    if res.inserts:
        ax1.set_ylim(0, max(counts) * 1.3)      # headroom for the legend
    ax1.set_title("Amplicon lengths in your reads")
    ax1.set_xlabel("Amplicon length (nt)")
    ax1.set_ylabel("Read pairs")
    ax1.legend(fontsize=8, loc="upper left")

    cands = res.candidates
    if cands:
        ax2.scatter([c.total_length for c in cands], [100 * c.yield_est for c in cands],
                    color="tab:blue", s=22, label="setting tested")
        best = {}
        for c in cands:
            best[c.total_length] = max(best.get(c.total_length, 0), c.yield_est)
        xs = sorted(best)
        ax2.plot(xs, [100 * best[x] for x in xs], color="tab:blue", lw=1)
    if cands:
        ys = [100 * c.yield_est for c in cands]
        low, high = max(0.0, min(ys) - 8.0), min(101.0, max(ys) + 3.0)
        low = min(low, max(0.0, high - 10.0))
    else:
        low, high = 0.0, 105.0
    left = res.floor - 6
    ax2.axvspan(left, res.floor, color="tab:red", alpha=0.12)
    ax2.text(left + 0.4, low + 0.05 * (high - low), "too short to\nmerge your longest\namplicons", fontsize=8,
             color="tab:red")
    if r is not None:
        ax2.scatter([r.total_length], [100 * r.yield_est], marker="*", s=220, color="tab:red",
                    zorder=3, label=f"recommended: {r.trim_f} + {r.trim_r}")
    ax2.set_xlim(left, max([c.total_length for c in cands] + [res.floor]) + 4)
    ax2.set_ylim(low, high)
    ax2.set_title("Usable read pairs by total trimmed length")
    ax2.set_xlabel("truncLenF + truncLenR (nt)")
    ax2.set_ylabel("Read pairs usable (%)")
    ax2.legend(fontsize=8, loc="lower right")

    prof = res.error_profile
    if prof:
        ef, er = ee_pair(res.cfg.max_ee)
        for key, color, name in (("f", "tab:blue", "forward"), ("r", "tab:orange", "reverse")):
            c = prof[key]
            x = range(1, c.shape[1] + 1)
            ax3.fill_between(x, c[0], c[2], color=color, alpha=0.15)
            ax3.plot(x, c[1], color=color, lw=1.4, label=f"{name} read: median, 10th-90th percentile")
        if ef == er:
            ax3.axhline(ef, color="0.3", ls=":", lw=1, label=f"maxEE = {_fmt(ef)}")
        else:
            ax3.axhline(ef, color="tab:blue", ls=":", lw=1, label=f"maxEE forward = {_fmt(ef)}")
            ax3.axhline(er, color="tab:orange", ls=":", lw=1, label=f"maxEE reverse = {_fmt(er)}")
        if r is not None:
            ax3.axvline(r.trim_f, color="tab:blue", ls="--", lw=1, label=f"cut: {r.trim_f} (F)")
            ax3.axvline(r.trim_r, color="tab:orange", ls="--", lw=1, label=f"cut: {r.trim_r} (R)")
        ax3.set_ylim(0, max(ef, er) * 2.2 + 0.5)
        ax3.set_title("Expected errors along the reads (dotted: maxEE)")
        ax3.set_xlabel("Position in read (nt)")
        ax3.set_ylabel("Expected errors so far")
        ax3.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _share_shorter(res: Result) -> float:
    """Share of measured amplicons shorter than the reads (those reads can run past the amplicon)."""
    if not res.inserts:
        return 0.0
    longest_read = max(res.sub.read_len_f, res.sub.read_len_r)
    return sum(1 for v in res.inserts if v < longest_read) / len(res.inserts)


def _central_range(inserts: List[int]):
    """The range holding 99% of the measured amplicon lengths (the extremes are often a few outliers)."""
    n = len(inserts)
    return inserts[int(0.005 * (n - 1))], inserts[int(0.995 * (n - 1))]


def _data_lines(res: Result) -> List[str]:
    sub, r = res.sub, res.recommendation
    lengths = (f"{sub.read_len_f} + {sub.read_len_r} nt" if sub.uniform_length
               else f"up to {sub.read_len_f} + {sub.read_len_r} nt")
    lines = [f"- Reads: {sub.n_pairs:,} pairs analysed ({lengths})."]
    if res.inserts:
        lo, hi = _central_range(res.inserts)
        median = res.inserts[len(res.inserts) // 2]
        span = (f"{res.inserts[0]}-{res.inserts[-1]} nt (median {median})" if (lo, hi) == (res.inserts[0], res.inserts[-1])
                else f"median {median} nt; 99% are {lo}-{hi} nt (all: {res.inserts[0]}-{res.inserts[-1]})")
        lines.append(f"- Amplicons: {span}, measured from the {res.frac_overlapping:.0%} of pairs "
                     f"that merge at full read length.")
        shorter = _share_shorter(res)
        if shorter >= 0.005 and sub.uniform_length:      # with reads of different lengths this cannot be told apart from trimming
            lines.append(f"- {shorter:.0%} of amplicons are shorter than the reads. If primers or adapters have "
                         f"not been removed, those reads run past the amplicon into them; remove them first.")
    else:
        lines.append(f"- Longest amplicon: {res.longest_insert} nt (as given).")
    qv = getattr(sub, "quality_values", 0)
    if 3 <= qv <= 5:
        lines.append(f"- Quality scores take only {qv} distinct values, as with binned scores on newer Illumina "
                     f"instruments. DADA2's error model can need adjusting for such data.")
    if r is not None and r.mean_overlap == r.mean_overlap:
        lines.append(f"- Overlap at the recommended settings: about {r.mean_overlap:.0f} nt on average.")
    return lines


def _data_summary(res: Result) -> dict:
    sub = res.sub
    has = bool(res.inserts)
    return {"pairs_analysed": sub.n_pairs, "read_length_f": sub.read_len_f, "read_length_r": sub.read_len_r,
            "uniform_read_length": bool(sub.uniform_length),
            "amplicon_min": res.inserts[0] if has else None,
            "amplicon_median": res.inserts[len(res.inserts) // 2] if has else None,
            "amplicon_max": res.inserts[-1] if has else None,
            "amplicon_99pct_min": _central_range(res.inserts)[0] if has else None,
            "amplicon_99pct_max": _central_range(res.inserts)[1] if has else None,
            "share_merged_at_full_length": round(res.frac_overlapping, 4),
            "share_amplicons_shorter_than_reads": round(_share_shorter(res), 4),
            "quality_values": getattr(sub, "quality_values", 0)}


def _split_sources(res: Result) -> str:
    spread = "a spread of forward/reverse splits among those with the highest filter pass rate"
    return spread


# Limits that both reports state, written once so that the two cannot drift apart.
_LIMIT_ACCURACY = "- Only read quality and mergeability are assessed, not taxonomic accuracy."
_LIMIT_DENOISING = ("- The yield is measured on the raw reads, before DADA2 denoises them. DADA2 cannot confirm a sequence it sees "
                    "about once, so in diverse samples (soil, for example) at modest depth the share of reads DADA2 keeps can be "
                    "far lower than this estimate. Use it to compare settings, not as the final retained fraction.")
_LIMIT_INPUT = "- The reads are assumed to have adapters and primers removed and no quality trimming."


def _len_text(x) -> str:
    """Mean cut length kept, forward / reverse (a dash when no pair is kept)."""
    return "-" if x.get("len_f") is None or x["len_f"] != x["len_f"] else f"{x['len_f']:.0f} / {x['len_r']:.0f}"


def _its_ee_section(res: Result) -> List[str]:
    rows = res.ee_table
    if not rows:
        return []
    out = ["## maxEE: retention against read quality", "",
           "ITS reads are not given a fixed truncation length, so maxEE is the one setting that trades reads for quality: each read is still cut at its first "
           "base of quality 2 or below, but a read with a noisy tail is not shortened, it is dropped when its expected errors exceed maxEE. For each value: the share of pairs that pass DADA2's filter, the estimated yield, how long the kept reads are, and how "
           "noisy they are.", "",
           "| maxEE | pairs passing the filter | estimated yield | mean length kept (F / R) | expected errors per kept read (F / R) | "
           "expected errors in the reads added beyond maxEE 2 (F / R) |", "|---|---|---|---|---|---|"]
    for x in rows:
        label = _ee_text(x["ee"]) + (" (used)" if x["current"] else "")
        label = f"**{label}**" if x["current"] else label
        out.append(f"| {label} | {x['pass_rate']:.1%} | {x['yield_est']:.1%} | {_len_text(x)} | {_load_text(x)} | {_added_text(x)} |")
    out += ["", *ee_advice(rows, res.cfg.max_ee, res.relaxed), "",
            "The yield is an upper bound: it stops before DADA2's denoising and chimera removal, so DADA2 keeps less, and the less so at low depth or high diversity "
            "(on a spring-wheat ITS study DADA2 merged 0.65 of the estimate at 5,000 pairs per sample and 0.87 at full depth). The gain from a looser maxEE is real but "
            "smaller than the estimate promises, and AmpliTrim does not measure what a looser maxEE does to the accuracy of your results.", ""]
    return out


def _its_loss_section(res: Result) -> List[str]:
    """ITS: where the pairs go, and what a forward-reads-only analysis would keep."""
    x = res.its_loss
    if not x:
        return []
    ee = _ee_text(res.cfg.max_ee)
    out = ["## Where the pairs go", "",
           f"Of the {x['pairs']:,} read pairs analysed, with maxEE {ee} and a minimum length of {res.min_len} nt (estimates; the merging is vsearch's):", "",
           "| step | share of the pairs |", "|---|---|",
           f"| too short after the cut at the first base of quality 2 or below (under {res.min_len} nt) | {x['too_short']:.1%} |",
           f"| contain an N | {x['has_n']:.1%} |",
           f"| more expected errors than maxEE {ee} | {x['over_maxEE']:.1%} |",
           f"| **pass the filter** | **{x['pass']:.1%}** |",
           f"| pass the filter but do not merge | {x['pass_not_merged']:.1%} |",
           f"| **pass the filter and merge** | **{x['pass_and_merged']:.1%}** |", "",
           f"- **{1 - x['overlap_at_full_length']:.0%} of the pairs do not merge even at full read length:** the amplicon is longer than your two reads can span, or the overlap is too poor. "
           f"No setting recovers these as pairs."]
    if x["forward_only_pass"] > x["pass_and_merged"]:
        out.append(f"- **Using the forward reads alone,** {x['forward_only_pass']:.1%} of the reads would pass the same filter, against {x['pass_and_merged']:.1%} of the pairs that pass and merge. "
                   f"DADA2's author says that using only the forward reads is \"totally reasonable\", and that when the reverse reads are bad enough they \"will cost more in sensitivity to "
                   f"low-frequency stuff than they add in a lower FP rate\" (DADA2 issue #236, 2017). AmpliTrim does not measure what this does to your ASVs.")
    return out + [""]


def _report_text_its(res: Result) -> str:
    cfg, r = res.cfg, res.recommendation
    out = ["# AmpliTrim report", ""]
    if r is None:
        out += ["## No recommendation", "", "No read pairs survived filtering and merging.", ""]
    else:
        s = _settings(res)
        out += [
            "## Recommended settings", "",
            "Do not truncate. ITS amplicons differ in length, so a fixed `truncLen` would discard every read "
            "shorter than it.", "",
            "| truncLenF | truncLenR | maxEE | minLen | minOverlap |", "|---|---|---|---|---|",
            f"| **0 (none)** | **0 (none)** | {_ee_text(s['maxEE'])} | {s['minLen']} | {s['minOverlap']} |", "",
            "In DADA2 (R):", "", "```r",
            "out <- filterAndTrim(fwd, filtF, rev, filtR, maxN = 0,",
            f"                     maxEE = c({_ee_f(s['maxEE'])}, {_ee_r(s['maxEE'])}), truncQ = 2, minLen = {s['minLen']})",
            "# ... learnErrors(), dada() ...",
            f"merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = {s['minOverlap']}, "
            f"maxMismatch = {s['maxMismatch']})", "```", "",
            "In QIIME 2:", "", "```",
            f"qiime dada2 denoise-paired ... --p-trunc-len-f 0 --p-trunc-len-r 0 \\\n    --p-max-ee-f {_ee_f(s['maxEE'])} "
            f"--p-max-ee-r {_ee_r(s['maxEE'])}", "```", "",
            "## What to expect", "",
            f"- About **{r.yield_est:.1%}** of read pairs should pass quality filtering and merge: "
            f"{r.ee_pass_rate:.1%} pass the filter, and {r.merge_rate_tolerant:.1%} of those merge.",
            f"- {_maxee_note(res)}", ""]
    if res.warnings:
        out += ["## Things to check", ""] + [f"- {w}" for w in res.warnings] + [""]
    out += _its_loss_section(res)
    out += _its_ee_section(res)
    out += ["## Your data", ""] + _data_lines(res) + [""]
    out += ["![summary](summary.png)", "",
            "## How these settings were chosen", "",
            f"- {res.sub.n_pairs:,} of {res.sub.n_input_pairs:,} read pairs were analysed.",
            f"- No fixed truncation length is tested, as in DADA2's ITS workflow: each read is cut at its first base of quality 2 or "
            f"below (`truncQ = 2`), then the pair is kept if both reads are at least {res.min_len} nt, free of Ns "
            f"and within maxEE {_ee_text(cfg.max_ee)}. The kept pairs are merged with vsearch.", "",
            "## Limits", "",
            "- The estimates come from merging raw reads with vsearch, not from DADA2. Run DADA2 once with these "
            "settings to confirm.",
            "- Pairs whose amplicon is longer than the reads can span cannot merge, with any setting.",
            _LIMIT_ACCURACY, _LIMIT_DENOISING, _LIMIT_INPUT, "",
            f"_AmpliTrim {__version__}_", ""]
    return "\n".join(out)


def _report_text(res: Result) -> str:
    if res.marker == "its":
        return _report_text_its(res)
    cfg, r, sub = res.cfg, res.recommendation, res.sub
    out = ["# AmpliTrim report", ""]
    if r is None:
        out += ["## No recommendation", "",
                "Nothing tested could merge your longest amplicons with usable yield.", ""]
    else:
        s = _settings(res)
        out += [
            "## Recommended settings", "",
            "| truncLenF | truncLenR | maxEE | minOverlap |", "|---|---|---|---|",
            f"| **{s['truncLenF']}** | **{s['truncLenR']}** | {_ee_text(s['maxEE'])} | {s['minOverlap']} |", "",
            "In DADA2 (R):", "", "```r",
            "out <- filterAndTrim(fwd, filtF, rev, filtR,",
            f"                     truncLen = c({s['truncLenF']}, {s['truncLenR']}), "
            f"maxEE = c({_ee_f(s['maxEE'])}, {_ee_r(s['maxEE'])}), truncQ = 2, maxN = 0)",
            "# ... learnErrors(), dada() ...",
            f"merged <- mergePairs(dadaF, filtF, dadaR, filtR, minOverlap = {s['minOverlap']}, "
            f"maxMismatch = {s['maxMismatch']})", "```", "",
            "In QIIME 2:", "", "```",
            f"qiime dada2 denoise-paired ... --p-trunc-len-f {s['truncLenF']} --p-trunc-len-r "
            f"{s['truncLenR']} \\\n    --p-max-ee-f {_ee_f(s['maxEE'])} --p-max-ee-r {_ee_r(s['maxEE'])}", "```", "",
            "## What to expect", "",
            f"- About **{r.yield_est:.1%}** of read pairs should pass quality filtering and merge.",
            f"- Your longest amplicon is {res.longest_insert} nt. These settings can merge amplicons up to "
            f"{r.max_mergeable_insert} nt, so your whole length range is covered.",
            f"- {_maxee_note(res)}", ""]
    if res.warnings:
        out += ["## Things to check", ""] + [f"- {w}" for w in res.warnings] + [""]
    out += _coverage_section(res)
    out += _ee_section(res)
    out += ["## Your data", ""] + _data_lines(res) + [""]
    out += ["![summary](summary.png)", "",
            "## How these settings were chosen", "",
            f"- {sub.n_pairs:,} of {sub.n_input_pairs:,} read pairs were analysed.",
            f"- Your longest amplicon is {res.longest_insert} nt, so the settings must keep at least "
            f"{res.floor} nt in total ({res.longest_insert} plus the minimum overlap of {cfg.dada2_min_overlap}). "
            f"Shorter totals are not tested.",
            f"- Totals were tested from {res.floor} nt in steps of {res.step} nt"
            + (" up to the read limit" if res.exhaustive else ", stopping once the share of usable pairs stopped "
                                                           "improving")
            + f". For each total, {_split_sources(res)} were quality-filtered (maxEE {_ee_text(cfg.max_ee)}) and merged "
            f"with vsearch.",
            f"- Settings within {100 * res.tolerance:g} percentage point{'' if res.tolerance == 0.01 else 's'} of "
            f"the best yield count as tied; the one with the shortest overlap was chosen.",
            "",
            "## Limits", "",
            "- The estimates come from merging raw reads with vsearch, not from DADA2. With long overlaps (well "
            "above 50 nt) DADA2's zero-mismatch rule can reject pairs that vsearch accepts, so the yield can "
            "be several points optimistic. Run DADA2 once with these settings to confirm.",
            "- Chimeras, and pairs that merge at the wrong offset, count as usable here, so the yield can be "
            "optimistic. This is most likely when amplicons have repetitive overlaps.",
            _LIMIT_ACCURACY, _LIMIT_DENOISING, _LIMIT_INPUT, "",
            f"_AmpliTrim {__version__}_", ""]
    return "\n".join(out)


def write_outputs(res: Result, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _figure(res, out_dir / "summary.png")
    r = res.recommendation
    (out_dir / "summary.json").write_text(encoding="utf-8", data=json.dumps({
        "schema_version": 1,                                  # keys may be added within a version, never removed or renamed
        **({"its": res.its_loss} if res.its_loss else {}),
        "maxEE_rule": {"relaxed": res.relaxed, "pass_rate_at_maxEE_2": res.pass_rate_default, "percentile": RULE_Q, "cap": RULE_CAP,
                       "cut_measured_at": None if res.rule is None else res.rule["cut"], "percentile_expected_errors": None if res.rule is None else res.rule["percentile"]},
        "recommended": None if r is None else {**_settings(res), "estimated_yield": round(r.yield_est, 4),
                                               "longest_mergeable_insert": r.max_mergeable_insert},
        "longest_amplicon": res.longest_insert, "pairs_analysed": res.sub.n_pairs,
        "data": _data_summary(res), "warnings": res.warnings,
        **({"maxEE_exploration": {"minRetention": res.min_retention,
                                  "suggested": None if _suggestion(res)[0] is None else
                                  {"maxEE": list(ee_pair(_suggestion(res)[0]["ee"])), "truncLenF": _suggestion(res)[0]["trim_f"],
                                   "truncLenR": _suggestion(res)[0]["trim_r"], "estimated_yield": round(_suggestion(res)[0]["yield_est"], 4),
                                   "pass_rate": round(_suggestion(res)[0]["pass_rate"], 4), "reachesFloor": _suggestion(res)[1]},
                                  "rows": [{"maxEE": list(ee_pair(x["ee"])), "truncLenF": x["trim_f"], "truncLenR": x["trim_r"],
                                            "pass_rate": None if x["pass_rate"] is None else round(x["pass_rate"], 4),
                                            "estimated_yield": None if x["yield_est"] is None else round(x["yield_est"], 4),
                                            "expected_errors_kept": None if x["ee_f"] is None else [round(x["ee_f"], 3), round(x["ee_r"], 3)],
                                            "meanLengthKept": None if x.get("len_f") is None or x["len_f"] != x["len_f"] else [round(x["len_f"], 1), round(x["len_r"], 1)],
                                            "beatenBy": [list(ee_pair(o["ee"])) for o in better]}
                                           for x, better in zip(res.ee_search, __import__("amplitrim.ranking", fromlist=["dominators"]).dominators(res.ee_search))]}}
                                 if res.ee_search else {}),
        **({"coverage_exploration": {"rows": [{"allowUnmergeable": x["level"], "longestCovered": x["longest"], "totalLength": x["total"],
                                               "truncLenF": x["trim_f"], "truncLenR": x["trim_r"],
                                               "leftOut": None if x["lost"] is None else round(x["lost"], 4),
                                               "pass_rate": None if x["pass_rate"] is None else round(x["pass_rate"], 4),
                                               "estimated_yield": None if x["yield_est"] is None else round(x["yield_est"], 4),
                                               "expected_errors_kept": None if x["ee_f"] is None else [round(x["ee_f"], 3), round(x["ee_r"], 3)]}
                                              for x in res.coverage_search]}} if res.coverage_search else {}),
        **({"maxEE_table": [{"maxEE": list(ee_pair(x["ee"])), "pass_rate": round(x["pass_rate"], 4), "estimated_yield": round(x["yield_est"], 4),
                             "expected_errors_kept": [round(x["ee_f"], 3), round(x["ee_r"], 3)],
                             "expected_errors_added": None if x["add_f"] != x["add_f"] else [round(x["add_f"], 3), round(x["add_r"], 3)]}
                            for x in res.ee_table]} if res.ee_table else {})}, indent=2))
    path = out_dir / "REPORT.md"
    path.write_text(_report_text(res), encoding="utf-8")
    return path
