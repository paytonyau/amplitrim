"""The project page (docs/index.html; the how-it-works page is on the validation archive's site) is self-contained and accessible, and says what the README says: every sentence, command and table row of the README is on the page, and the page adds no numbers of its own. No network."""
import re
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAGE_FILE, README_FILE, GUIDE_FILE = ROOT / "docs" / "index.html", ROOT / "README.md", ROOT / "docs" / "usage.md"
BLOB = "https://github.com/paytonyau/amplitrim/blob/main/"
VOID = {"meta", "link", "img", "br", "hr", "input"}
IMAGE_HOSTS = ("https://img.shields.io/", "https://github.com/paytonyau/amplitrim/")
BLOCKS = {"p", "li", "ul", "ol", "h1", "h2", "h3", "h4", "tr", "th", "td", "div", "section", "header", "footer", "nav", "main", "pre", "summary", "details", "table", "svg", "text", "title", "figure", "figcaption"}      # picture labels are separate pieces of text
norm = lambda s: " ".join(s.split())


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.stack, self.errors, self.headings, self.ids, self.links, self.imgs = [], [], [], [], [], []
        self.scripts, self.sheets, self.pres, self.tables, self.flow = [], [], [], [], []
        self.text, self.attrs_html, self.html_lang, self.title = [], 0, None, ""
        self._skip, self._cell, self._row, self._table, self._pre, self._flow = 0, None, None, None, None, False
        self._in_body, self._sr, self._link, self._step = False, 0, None, False
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag not in VOID: self.stack.append(tag)
        if "id" in a: self.ids.append(a["id"])
        if tag == "html": self.html_lang = a.get("lang")
        if tag == "body": self._in_body = True
        if tag in BLOCKS and self._in_body: self.text.append(" ")                                          # block boundaries are spaces; inline tags add nothing
        if tag in ("script", "style"): self._skip += 1
        if tag == "script": self.scripts.append(a)
        if tag == "link": self.sheets.append(a)
        if tag == "img": self.imgs.append(a)
        if re.fullmatch(r"h[1-6]", tag): self.headings.append(int(tag[1]))
        if tag == "a":
            self.links.append({"href": a.get("href", ""), "text": ""}); self._link = self.links[-1]
        if tag == "img" and self._link is not None: self._link["text"] += a.get("alt", "")
        if {"sr", "num"} & set(a.get("class", "").split()): self._sr += 1                                 # screen-reader-only text and step badges are not prose
        if tag == "ol" and "flow" in a.get("class", "").split(): self._flow = True
        if tag == "span" and "step" in a.get("class", "").split() and self._flow: self.flow.append(""); self._step = True
        if tag == "table": self._table = []
        if tag == "tr" and self._table is not None: self._row = []
        if tag in ("th", "td") and self._row is not None: self._cell = []
        if tag == "br" and self._cell is not None: self._cell.append("\n")
        if tag == "pre": self._pre = []

    def handle_endtag(self, tag):
        if tag in VOID: return
        if tag in BLOCKS and self._in_body: self.text.append(" ")
        if not self.stack or self.stack[-1] != tag: self.errors.append(f"</{tag}> closes {self.stack[-1:] or 'nothing'}")
        if self.stack and tag in self.stack:
            while self.stack and self.stack.pop() != tag: pass
        if tag in ("script", "style"): self._skip -= 1
        if tag == "a": self._link = None
        if tag == "ol": self._flow = False
        if tag == "span" and self._sr: self._sr -= 1
        if tag == "span": self._step = False
        if tag in ("th", "td") and self._cell is not None:
            lines = tuple(x for x in (norm(p) for p in "".join(self._cell).split("\n")) if x); self._row.append(lines); self._cell = None
        if tag == "tr" and self._row is not None: self._table.append(tuple(self._row)); self._row = None
        if tag == "table": self.tables.extend(self._table); self._table = None
        if tag == "pre": self.pres.append("".join(self._pre)); self._pre = None

    def handle_data(self, data):
        if self._skip: return
        if not self._in_body:
            if self.stack and self.stack[-1] == "title": self.title += data
            return
        if self._sr: return
        self.text.append(data)
        if self._cell is not None: self._cell.append(data)
        if self._pre is not None: self._pre.append(data)
        if self._link is not None: self._link["text"] += data
        if self._step: self.flow[-1] += data

    @property
    def body_text(self): return norm("".join(self.text))


def plain(md): return re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", re.sub(r"\*\*([^*]+)\*\*", r"\1", re.sub(r"`([^`]*)`", r"\1", md.replace("\\*", "\u0001")))).replace("\u0001", "*")


def readme_parts():
    text = README_FILE.read_text(encoding="utf-8")
    blocks = re.findall(r"```[a-z]*\n(.*?)```", text, flags=re.S)
    summary = norm(re.search(r"<summary>(.*?)</summary>", text, flags=re.S).group(1))
    body = re.sub(r"```.*?```", "", text, flags=re.S)
    body = re.sub(r"<details>.*?</details>", lambda m: "\n", body, flags=re.S)
    sentences, rows = [], []
    for line in body.splitlines():
        if line.startswith("|"):
            if not re.match(r"^\|[-| ]+\|$", line):
                rows.append(tuple(tuple(x for x in (norm(plain(p)) for p in c.strip().split("<br>")) if x) for c in line.strip("|").split("|")))
        elif line.strip() and not line.startswith(("#", "[![")):
            sentences += [norm(s) for s in re.split(r"(?<=[.!?])\s+(?=[A-Z])", plain(re.sub(r"^- ", "", line)))]
    for line in re.search(r"<details>.*?</details>", text, flags=re.S).group(0).splitlines():          # the collapsed table is part of the README too
        if line.startswith("|") and not re.match(r"^\|[-| ]+\|$", line):
            rows.append(tuple(tuple(x for x in (norm(plain(p)) for p in c.strip().split("<br>")) if x) for c in line.strip("|").split("|")))
    return text, blocks, summary, sentences, rows


PAGE = Page(PAGE_FILE.read_text(encoding="utf-8"))
ALL = [("home", PAGE_FILE, PAGE)]
HOW_URL = "https://paytonyau.github.io/amplitrim-validation/how-it-works/"        # the how-it-works page lives on the validation archive's site


def test_the_pages_are_self_contained_and_well_formed():
    for name, f, p in ALL:
        raw = f.read_text(encoding="utf-8")
        assert not p.errors and not p.stack, (name, p.errors, p.stack)                                  # every tag opened is closed
        assert p.html_lang == "en" and p.title.strip() and '<meta name="viewport"' in raw, name
        assert all("src" not in x for x in p.scripts), name                                             # only inline scripts
        assert all(l.get("rel") == "icon" and l["href"].startswith("data:") for l in p.sheets), name    # no external stylesheet or font
        assert all(i["src"].startswith(IMAGE_HOSTS) for i in p.imgs), (name, [i["src"] for i in p.imgs])   # the badges are the only images
        assert not re.search(r"<(iframe|form|object|embed)\b", raw) and "http://" not in raw.replace("http://www.w3.org/2000/svg", ""), name


def test_the_pages_are_accessible_in_the_basics():
    for name, f, p in ALL:
        assert p.headings.count(1) == 1, name                                                           # one h1
        assert all(b - a <= 1 for a, b in zip(p.headings, p.headings[1:])), (name, p.headings)          # no skipped heading level
        assert all(i.get("alt", "").strip() for i in p.imgs), name
        assert all(norm(l["text"]) for l in p.links), (name, [l["href"] for l in p.links if not norm(l["text"])])
        assert len(p.ids) == len(set(p.ids)) and p.links[0]["href"] == "#main", name                    # unique ids; the skip link comes first
        assert all(l["href"][1:] in p.ids for l in p.links if l["href"].startswith("#")), name          # every in-page link has a target


def test_every_link_goes_somewhere_that_exists():
    guide = GUIDE_FILE.read_text(encoding="utf-8")
    slug = lambda h: re.sub(r"\s+", "-", re.sub(r"[^\w\s-]", "", h.strip().lower()))
    guide_anchors = {slug(re.sub(r"^#+\s*", "", l)) for l in guide.splitlines() if l.startswith("#")}
    for name, f, p in ALL:
        for l in p.links:
            href = l["href"]
            if href.startswith("#"):
                continue
            if not href.startswith("https://"):                                                         # a link to the other page of the site
                target = (f.parent / href.split("#")[0]).resolve()
                if href.endswith("/") or target.is_dir():
                    target = target / "index.html"
                assert target.exists() and (ROOT / "docs").resolve() in target.parents, (name, href)
                continue
            if href.startswith(BLOB):
                path, _, anchor = href[len(BLOB):].partition("#")
                assert (ROOT / path).exists(), (name, href)
                if anchor:
                    assert path == "docs/usage.md" and anchor in guide_anchors, (name, href)


def test_everything_the_readme_says_is_on_the_page():
    _, blocks, summary, sentences, rows = readme_parts()
    page_text = PAGE.body_text
    missing = [s for s in sentences if s not in page_text]
    assert not missing, missing[:3]                                                                     # every sentence of the README, caveats included
    assert [norm(b) for b in blocks if "→" not in b] and all(any(norm(b) == norm(p) for p in PAGE.pres) for b in blocks if "→" not in b)   # every command and the example output, as written
    assert PAGE.flow == [s.strip() for s in [b for b in blocks if "→" in b][0].split("→")]               # the same steps, in the same order
    assert summary in page_text
    assert Counter(rows) == Counter(PAGE.tables), (set(rows) ^ set(PAGE.tables))                       # the same table rows, no more and no fewer


def test_the_home_page_adds_no_numbers_of_its_own():
    text, _, _, _, _ = readme_parts()
    number = re.compile(r"\d+(?:\.\d+)*")
    extra = set(number.findall(PAGE.body_text)) - set(number.findall(text))
    assert not extra, extra                                                                             # a number on the page must be in the README


def test_no_page_makes_a_claim_we_decided_not_to_make():
    for name, f, p in ALL:
        raw = f.read_text(encoding="utf-8").lower().split("<style>")[0]
        for word in ("FIGARO", "guarantee", "resolution", "fastest", "faster", "state-of-the-art", "outperform", "best-in-class"):
            assert word.lower() not in p.body_text.lower() and word.lower() not in raw, (name, word)


def test_the_why_cards_are_two_by_two():
    raw = PAGE_FILE.read_text(encoding="utf-8"); css = raw[raw.index("<style>"):raw.index("</style>")]
    assert raw.count('<li class="card">') == 4
    assert re.search(r"\.cards\{[^}]*grid-template-columns:repeat\(2,minmax\(0,1fr\)\)", css)               # two columns, so four cards make two rows
    assert "@media (max-width:40rem){.cards{grid-template-columns:1fr}}" in css                            # one column on a phone


def test_the_title_has_a_coloured_trim_and_the_wastewater_footnote_is_explained():
    raw = PAGE_FILE.read_text(encoding="utf-8"); css = raw[raw.index("<style>"):raw.index("</style>")]
    assert '<h1>Ampli<span class="trim">Trim</span></h1>' in raw and re.search(r"\.trim\{[^}]*color:var\(--accent\)", css)
    assert any("wastewater (1 dataset): behind by 1.0 point*" in line for row in PAGE.tables for cell in row for line in cell)     # the asterisk on the result
    assert re.search(r'<p class="note">\* Wastewater is shown on its own\.', raw) and "which cost about a point" in PAGE.body_text      # and the footnote that explains it


def theme_vars(raw):
    """The colour variables of the light theme, of the dark theme chosen by the system, and of the dark theme chosen with the switch."""
    css = raw[raw.index("<style>"):raw.index("</style>")]
    pairs = lambda block: dict(re.findall(r"--([a-z-]+):(#[0-9a-fA-F]{3,6})", block))
    light = {}
    for b in re.findall(r"(?m)^:root\{([^}]*)\}", css): light.update(pairs(b))
    auto, chosen = {}, {}
    for b in re.findall(r'@media \(prefers-color-scheme:dark\)\{:root:not\(\[data-theme="light"\]\)\{([^}]*)\}\}', css): auto.update(pairs(b))
    for b in re.findall(r':root\[data-theme="dark"\]\{([^}]*)\}', css): chosen.update(pairs(b))
    return light, auto, chosen


def luminance(h):
    h = h.lstrip("#"); h = h if len(h) == 6 else "".join(c * 2 for c in h)
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = sorted((luminance(a), luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def test_the_themes_are_soft_readable_and_can_be_switched():
    for name, f, p in ALL:
        raw = f.read_text(encoding="utf-8")
        light, auto, chosen = theme_vars(raw)
        assert auto and auto == chosen, name                                                             # the system's dark theme and the chosen one cannot differ
        assert luminance(light["bg"]) < 0.97 and luminance(chosen["bg"]) > 0.012, name                    # neither pure white nor near-black
        assert light["bg"].lower() not in ("#fff", "#ffffff", "#000", "#000000") and chosen["bg"].lower() not in ("#000", "#000000")
        for label, theme in (("light", light), ("dark", chosen)):
            for fg in ("fg", "muted", "accent", "warn", "ok", "rev"):
                if fg not in theme: continue
                for bg in ("bg", "card"):
                    assert contrast(theme[fg], theme[bg]) >= 4.5, (name, label, fg, bg, round(contrast(theme[fg], theme[bg]), 2))      # readable text and marks
            assert contrast(theme["accent-fg"], theme["accent"]) >= 4.5, (name, label)                   # the text on a filled button
        assert 'localStorage.getItem(\'amplitrim-theme\')' in raw and "theme-toggle" in raw and "prefers-color-scheme" in raw, name     # remembered choice and the switch
        head_scripts = re.findall(r"<script>(.*?)</script>", raw[:raw.index("<style>")], flags=re.S)
        assert any("localStorage.getItem('amplitrim-theme')" in x and "setAttribute('data-theme'" in x and x.count("{") == x.count("}") for x in head_scripts), name     # a whole script before the first paint: no flash


def test_the_home_page_links_to_the_how_it_works_page_on_the_archive_site():
    assert [l["href"] for l in PAGE.links].count(HOW_URL) == 2                                          # the menu and the button
    assert not (ROOT / "docs" / "how-it-works").exists()                                                # the page is not in this repository


def test_the_readme_and_the_home_page_say_that_the_maxee_rule_has_no_upper_limit_and_notes_a_limit_above_5():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "never below 2, with a note if it goes above 5" in readme and "never above 5" not in readme
    assert "never above" not in PAGE.body_text                                                          # the home page is made from the README
