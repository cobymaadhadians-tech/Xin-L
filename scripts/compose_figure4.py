#!/usr/bin/env python3
"""Compose the quantitative and biological-summary panels of Figure 4."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import re

from lxml import etree
import cairosvg


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
INKSCAPE_NS = "http://www.inkscape.org/namespaces/inkscape"
NS = {"s": SVG_NS}

FIGURE_STEM = "Figure_4_SCZ_STRATIFIED_PI_BRAIN_CELL"
PANEL_A_WIDTH = 418.986421
PANEL_A_HEIGHT = 334.118556
PANEL_B_SCALE = 0.225
PANEL_B_X = 431.0
PANEL_B_Y = 13.0 - (90.0 * PANEL_B_SCALE)
FIGURE_WIDTH = 662.0
FIGURE_HEIGHT = PANEL_A_HEIGHT

# The source panels remain SVG text elements and CairoSVG preserves them as
# editable text in the vector PDF. These settings document the equivalent
# publication-safe vector-text contract used by the figure workflow.
EDITABLE_TEXT_SETTINGS = {
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
}
FONT_SIZE_PT = 8.0


def qname(local: str) -> str:
    return f"{{{SVG_NS}}}{local}"


def parse_svg(path: Path) -> etree._Element:
    parser = etree.XMLParser(remove_blank_text=False)
    return etree.parse(path, parser).getroot()


def normalize_panel_b_styles(root: etree._Element) -> None:
    for text_node in root.xpath(".//s:text", namespaces=NS):
        if text_node.get("class") == "st31" and "".join(text_node.itertext()).strip() == "b":
            tspans = text_node.xpath("./s:tspan", namespaces=NS)
            if tspans:
                tspans[0].text = "(b)"
            else:
                text_node.text = "(b)"
        if text_node.get("class") == "st18" and text_node.get("transform") == "translate(105 89)":
            text_node.set("transform", "translate(145 89)")
    for style in root.xpath(".//s:style", namespaces=NS):
        text = style.text or ""
        text = text.replace("TimesNewRomanPS-BoldMT, 'Times New Roman'", "Arial, Helvetica, sans-serif")
        text = text.replace("TimesNewRomanPSMT, 'Times New Roman'", "Arial, Helvetica, sans-serif")
        text = text.replace("Arial-BoldMT, Arial", "Arial, Helvetica, sans-serif")
        text = text.replace("#111824", "#262626")
        text = text.replace("#364351", "#262626")
        # Panel a establishes the final-size hierarchy: 8 px panel labels,
        # 7 px headings, 5.8 px bold row labels and 5.7 px supporting text.
        # Panel b is placed at 0.225 scale, so its source SVG sizes are
        # normalized to the same final sizes before composition.
        size_map = {
            ".st18": 31.1,
            ".st19": 23.1,
            ".st20": 25.3,
            ".st22": 25.3,
            ".st25": 25.3,
            ".st28": 31.1,
            ".st30": 25.8,
            ".st31": 35.6,
        }
        for selector, size in size_map.items():
            text = re.sub(
                rf"({re.escape(selector)}\s*\{{\s*font-size:\s*)[0-9.]+px",
                rf"\g<1>{size:g}px",
                text,
            )
        # The remaining source rules place font-size after other properties
        # or group several selectors in one rule, so normalize those exact
        # source values as well.
        for old, new in {
            "17.5": "23.1",
            "23": "25.3",
            "27": "25.3",
            "28": "25.8",
            "33.4": "31.1",
            "35": "31.1",
            "43": "35.6",
        }.items():
            text = text.replace(f"font-size: {old}px;", f"font-size: {new}px;")
        style.text = text


def append_defs(parent: etree._Element, source: etree._Element) -> None:
    source_defs = source.find(qname("defs"))
    if source_defs is None:
        return
    target_defs = parent.find(qname("defs"))
    if target_defs is None:
        target_defs = etree.SubElement(parent, qname("defs"))
    for child in source_defs:
        target_defs.append(deepcopy(child))


def append_panel_content(parent: etree._Element, source: etree._Element, transform: str) -> None:
    group = etree.SubElement(parent, qname("g"), transform=transform)
    for child in source:
        if not isinstance(child.tag, str):
            continue
        local = etree.QName(child).localname
        if local in {"metadata", "defs"}:
            continue
        group.append(deepcopy(child))


def build_composite(panel_a: Path, panel_b: Path) -> bytes:
    a = parse_svg(panel_a)
    b = parse_svg(panel_b)
    normalize_panel_b_styles(b)

    root = etree.Element(
        qname("svg"),
        nsmap={None: SVG_NS, "xlink": XLINK_NS, "inkscape": INKSCAPE_NS},
        version="1.1",
        width=f"{FIGURE_WIDTH}pt",
        height=f"{FIGURE_HEIGHT}pt",
        viewBox=f"0 0 {FIGURE_WIDTH} {FIGURE_HEIGHT}",
    )
    defs = etree.SubElement(root, qname("defs"))
    append_defs(root, a)
    append_defs(root, b)
    style = etree.SubElement(defs, qname("style"), type="text/css")
    style.text = "text { font-family: Arial, Helvetica, 'DejaVu Sans', sans-serif !important; }"
    etree.SubElement(
        root,
        qname("rect"),
        x="0",
        y="0",
        width=str(FIGURE_WIDTH),
        height=str(FIGURE_HEIGHT),
        fill="white",
    )
    append_panel_content(root, a, "translate(0 0) scale(1)")
    append_panel_content(root, b, f"translate({PANEL_B_X} {PANEL_B_Y}) scale({PANEL_B_SCALE})")
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-a-svg", type=Path, required=True)
    parser.add_argument("--panel-b-svg", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    svg_bytes = build_composite(args.panel_a_svg, args.panel_b_svg)
    svg_path = args.outdir / f"{FIGURE_STEM}.svg"
    pdf_path = args.outdir / f"{FIGURE_STEM}.pdf"
    png_path = args.outdir / f"{FIGURE_STEM}.png"
    svg_path.write_bytes(svg_bytes)
    cairosvg.svg2pdf(bytestring=svg_bytes, write_to=str(pdf_path))
    cairosvg.svg2png(bytestring=svg_bytes, write_to=str(png_path), output_width=2648, output_height=1336)
    print(f"wrote {svg_path}")
    print(f"wrote {pdf_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
