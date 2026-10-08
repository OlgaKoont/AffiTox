#!/usr/bin/env python3
"""Convert existing Matplotlib SVG figures to BMC-sized vector PDFs.

Only figures independent of PoseBusters are exported: correlations, enrichment,
and inferential statistics. PDFs are written next to their source SVG files.
"""

from __future__ import annotations

import argparse
import base64
import io
import re
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import cairosvg
from PIL import Image


MM_PER_INCH = 25.4
POINTS_PER_INCH = 72.0
FULL_WIDTH_MM = 170.0
HALF_WIDTH_MM = 85.0
MAX_HEIGHT_MM = 225.0
MIN_FINAL_LINE_WIDTH_PT = 0.25
INCLUDED_SECTIONS = ("correlations", "enrichment", "inferential")
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"


def _number(value: str) -> float:
    match = re.match(r"\s*([0-9.+\-eE]+)", value)
    if not match:
        raise ValueError(f"Cannot parse SVG length: {value!r}")
    return float(match.group(1))


def _replace_style_stroke_widths(style: str, minimum: float) -> str:
    pattern = re.compile(r"(stroke-width\s*:\s*)([0-9.+\-eE]+)")

    def replace(match: re.Match[str]) -> str:
        width = float(match.group(2))
        if width <= 0:
            return match.group(0)
        return f"{match.group(1)}{max(width, minimum):.6g}"

    return pattern.sub(replace, style)


def _rgba_css(rgba: tuple[int, int, int, int]) -> tuple[str, str | None]:
    red, green, blue, alpha = rgba
    color = f"#{red:02x}{green:02x}{blue:02x}"
    opacity = None if alpha == 255 else f"{alpha / 255.0:.6g}"
    return color, opacity


def _vectorize_embedded_colorbars(root: ET.Element) -> None:
    """Replace Matplotlib's raster colorbar gradients with vector rectangles."""
    parent_by_child = {child: parent for parent in root.iter() for child in parent}
    images = list(root.iter(f"{{{SVG_NS}}}image"))
    for image_element in images:
        href = image_element.get(f"{{{XLINK_NS}}}href", "")
        if not href.startswith("data:image/png;base64,"):
            raise ValueError("SVG contains a non-PNG or externally linked raster image")

        png = Image.open(io.BytesIO(base64.b64decode(href.split(",", 1)[1]))).convert("RGBA")
        pixel_width, pixel_height = png.size
        long_side = max(pixel_width, pixel_height)
        short_side = min(pixel_width, pixel_height)
        if short_side == 0 or long_side / short_side < 2.0:
            raise ValueError(
                f"Embedded raster is not a colorbar ({pixel_width}x{pixel_height})"
            )

        x = _number(image_element.get("x", "0"))
        y = _number(image_element.get("y", "0"))
        width = _number(image_element.get("width", str(pixel_width)))
        height = _number(image_element.get("height", str(pixel_height)))
        strips = min(256, long_side)

        group = ET.Element(f"{{{SVG_NS}}}g")
        transform = image_element.get("transform")
        if transform:
            group.set("transform", transform)

        vertical = pixel_height >= pixel_width
        for index in range(strips):
            if vertical:
                pixel_y = min(int((index + 0.5) * pixel_height / strips), pixel_height - 1)
                rgba = png.getpixel((pixel_width // 2, pixel_y))
                rect_x = x
                rect_y = y + height * index / strips
                rect_width = width
                rect_height = height / strips * 1.002
            else:
                pixel_x = min(int((index + 0.5) * pixel_width / strips), pixel_width - 1)
                rgba = png.getpixel((pixel_x, pixel_height // 2))
                rect_x = x + width * index / strips
                rect_y = y
                rect_width = width / strips * 1.002
                rect_height = height

            color, opacity = _rgba_css(rgba)
            rect = ET.SubElement(group, f"{{{SVG_NS}}}rect")
            rect.set("x", f"{rect_x:.8g}")
            rect.set("y", f"{rect_y:.8g}")
            rect.set("width", f"{rect_width:.8g}")
            rect.set("height", f"{rect_height:.8g}")
            rect.set("fill", color)
            if opacity is not None:
                rect.set("fill-opacity", opacity)

        parent = parent_by_child[image_element]
        position = list(parent).index(image_element)
        parent.remove(image_element)
        parent.insert(position, group)


def _prepare_svg(svg_path: Path) -> tuple[bytes, float, float]:
    tree = ET.parse(svg_path)
    root = tree.getroot()
    _vectorize_embedded_colorbars(root)
    view_box = root.get("viewBox")
    if not view_box:
        raise ValueError(f"SVG has no viewBox: {svg_path}")
    _, _, view_width, view_height = (float(x) for x in view_box.split())
    aspect = view_height / view_width

    target_width_mm = FULL_WIDTH_MM
    target_height_mm = target_width_mm * aspect
    if target_height_mm > MAX_HEIGHT_MM:
        target_width_mm = HALF_WIDTH_MM
        target_height_mm = target_width_mm * aspect
    if target_height_mm > MAX_HEIGHT_MM:
        raise ValueError(
            f"Figure is too tall for BMC dimensions: {svg_path} "
            f"({target_width_mm:.1f} x {target_height_mm:.1f} mm)"
        )

    target_width_pt = target_width_mm / MM_PER_INCH * POINTS_PER_INCH
    scale = target_width_pt / view_width
    minimum_source_width = MIN_FINAL_LINE_WIDTH_PT / scale

    root.set("width", f"{target_width_mm:.6g}mm")
    root.set("height", f"{target_height_mm:.6g}mm")
    for element in root.iter():
        style = element.get("style")
        if style:
            element.set(
                "style", _replace_style_stroke_widths(style, minimum_source_width)
            )
        stroke_width = element.get("stroke-width")
        if stroke_width is not None:
            width = _number(stroke_width)
            if 0 < width < minimum_source_width:
                element.set("stroke-width", f"{minimum_source_width:.6g}")

    return ET.tostring(root, encoding="utf-8", xml_declaration=True), target_width_mm, target_height_mm


def _convert(svg_path: Path) -> tuple[Path, float, float]:
    svg_bytes, width_mm, height_mm = _prepare_svg(svg_path)
    pdf_path = svg_path.with_suffix(".pdf")
    cairosvg.svg2pdf(bytestring=svg_bytes, write_to=str(pdf_path))
    return pdf_path, width_mm, height_mm


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("/mnt/tank/scratch/okonovalova/AffiTox/analysis/figures"),
    )
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()

    svg_files = sorted(
        path
        for section in INCLUDED_SECTIONS
        for path in (args.figures_dir / section).rglob("*.svg")
    )
    if not svg_files:
        raise SystemExit(f"No eligible SVG figures under {args.figures_dir}")

    print(f"Converting {len(svg_files)} vector figures with {args.workers} workers")
    failures: list[tuple[Path, Exception]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_convert, path): path for path in svg_files}
        for future in as_completed(futures):
            source = futures[future]
            try:
                pdf_path, width_mm, height_mm = future.result()
                print(
                    f"{pdf_path}: {width_mm:.1f} x {height_mm:.1f} mm",
                    flush=True,
                )
            except Exception as exc:
                failures.append((source, exc))
                print(f"ERROR {source}: {exc}", flush=True)

    if failures:
        raise SystemExit(f"{len(failures)} figure(s) failed conversion")
    print(f"Completed: {len(svg_files)} PDFs")


if __name__ == "__main__":
    main()
