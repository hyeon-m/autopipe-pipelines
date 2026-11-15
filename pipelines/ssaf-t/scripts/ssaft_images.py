"""The image downsampling uses (image_key), apart from the map's dimension, and the element names.
Everything the pipeline does about this is in this file.

SSAF-T takes the two apart: `downsample_vector(z=2|3, image_key=...)`, where z is
the dimension of the map (2D or 3D) and image_key names ONE image element, either a 2D image (one
tissue mask; in a 3D map SSAF-T places it at the middle of the transcripts' z range) or a z-stack
built with Analysis.zstack_image (one mask per plane). So:

  * a 2D map uses a 2D image only;
  * a 3D map uses a 2D image or a z-stack. A z-stack has its planes chosen (z_planes, from the
    molecules per plane of Analysis.zstack_planes) and is built from those planes only. A 2D image
    in a 3D map needs no planes, which makes 3D possible without a z-stack (e.g. CosMx composite).

The candidates a raw folder offers: the 2D image its reader loads (Xenium morphology_focus, MERSCOPE
the z_layers mosaic <folder>_z<N>, CosMx composite) and its z-stack (Xenium morphology.ome.tif,
MERSCOPE the <stain> mosaics of every plane, or zstack_source), named as SSAF-T's zstack_image names
them (`morphology_zstack`, `DAPI_zstack`). A SpatialData zarr offers its own 2D and 3D image elements.
The recommendation: a 2D map, the platform's 2D image; a 3D map, the z-stack with the most planes
when there is one with several, else the 2D image.

MERSCOPE element names: <data folder>_z<N> and <data folder>_transcripts, from the data folder's own
name, whatever path the folder is reached by in the container (a link's target, or the folder
mounted as /input itself, whose name is read from the mount table). spatialdata_io names them after
the parent folder too (<parent>_<folder>_z<N>); they are renamed right after loading.

Functions marked [SSAF-T] mirror what SSAF-T does or names; plan.py and ssaft_load.py only call them
through this module.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

# the 2D image each platform's reader loads
ONE_LAYER = {"xenium": "morphology_focus", "cosmx": "composite"}

# 2D images preferred per platform, and why
PREFERRED_2D = {
    "xenium": (("morphology_focus", "morphology_mip"),
               "the Xenium morphology image projected over z (focus), the one Xenium makes for the 2D tissue"),
    "cosmx": (("composite",), "the CosMx composite image"),
}

MODE = {2: "2D image", 3: "z-stack"}


class ImageError(ValueError):
    pass


# ============================================================================ element names
def _mount_root(p: Path) -> str | None:
    """The source path of the mount at p (a folder bind-mounted into the container), if p is one."""
    try:
        if not os.path.ismount(p):
            return None
        root = None
        with open("/proc/self/mountinfo") as fh:
            for line in fh:
                f = line.split()
                if len(f) > 4 and f[4].encode().decode("unicode_escape") == str(p):
                    root = f[3].encode().decode("unicode_escape")       # the last mount there wins
        return root
    except (OSError, UnicodeDecodeError):
        return None


def folder_name(dataset_path) -> str:
    """The data folder's own name: the path's last part, or for a folder mounted as /input itself the
    name of the folder mounted there."""
    p = Path(dataset_path)
    root = _mount_root(p)
    return Path(root).name if root and Path(root).name else p.name


def merscope_reader_id(dataset_path) -> str:
    """The id spatialdata_io.merscope (and ssaft_load.merscope_load) names the elements with."""
    p = Path(dataset_path)
    return f"{p.parent.name}_{p.name}"


def rename_merscope(sdata, cfg):
    """[SSAF-T] After spatialdata_io.merscope: <parent>_<folder>_* images and points renamed
    <folder>_* (cfg['merscope_images'], set by plan.py)."""
    from ssaft_core import log
    M = cfg.get("merscope_images") or {}
    rid, name = M.get("reader_id"), M.get("name")
    if not rid or not name or rid == name:
        return sdata
    for kind in ("images", "points"):
        els = getattr(sdata, kind)
        for k in [k for k in els if k.startswith(rid + "_")]:
            new = name + k[len(rid):]
            els[new] = els.pop(k)
            log(f"{kind}[{k!r}] -> {kind}[{new!r}] (named after the data folder)")
    return sdata


def zstack_key(source, stain=None) -> str:
    """[SSAF-T] The key Analysis.zstack_image stores a z-stack under: the stain's (or the file's) name
    + _zstack (morphology.ome.tif -> morphology_zstack; stain DAPI -> DAPI_zstack)."""
    if stain:
        base = str(stain)
    else:
        src = source[0] if isinstance(source, (list, tuple)) and len(source) else source
        base = Path(str(src)).name
        base = re.sub(r"\[Z=\d+\]$", "", base)
        for suf in (".ome.tif", ".ome.tiff", ".tif", ".tiff"):
            if base.lower().endswith(suf):
                base = base[:-len(suf)]
                break
        base = re.sub(r"[_-]?z[\d*?]+$", "", base, flags=re.I)
        base = re.sub(r"[*?\[\]]", "", base).strip("_-")
    return f"{base or 'image'}_zstack"


# ======================================================================= plan time (dry_run)
def _tif_yx(f: Path):
    """(height, width) of a TIFF from its header, without reading the pixels."""
    try:
        import tifffile
        with tifffile.TiffFile(str(f)) as t:
            return list(t.pages[0].shape[-2:])
    except Exception:  # noqa: BLE001
        pass
    try:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None
        with Image.open(f) as im:
            return [im.size[1], im.size[0]]
    except Exception:  # noqa: BLE001
        return None


def raw_image_rows(plat, dataset_path, insp, planes, sel, *, name=None, stains=(), stain=None,
                   mask_z=None, source=None, stack_name=None) -> list[dict]:
    """[SSAF-T] The images a raw Xenium / MERSCOPE / CosMx folder offers for downsampling: the 2D
    image its reader loads, and the z-stack (when it has several planes)."""
    p = Path(dataset_path)
    rows = []
    yx = None
    if plat == "merscope":
        yx = _tif_yx(p / "images" / f"mosaic_{stain}_z{mask_z}.tif")
        rows.append({"name": f"{name}_z{mask_z}", "dims": 2, "axes": "cyx",
                     "shape": [len(stains) or 1] + yx if yx else None, "channels": list(stains),
                     "origin": f"2D image: the z{mask_z} mosaic spatialdata_io.merscope loads (merscope_mask_z), "
                               f"its {stain} channel makes the mask"})
        chans = [stain]
    elif plat == "cosmx":
        rows.append({"name": "composite", "dims": 2, "axes": "cyx", "shape": None, "channels": [],
                     "origin": "2D image: the CellComposite FOV images joined into one slide by the CosMx slide loader"})
        chans = []
    else:
        focus = (insp.get("files") or {}).get("morphology_focus") or []
        yx = _tif_yx(p / focus[0]) if focus else None
        rows.append({"name": "morphology_focus", "dims": 2, "axes": "cyx",
                     "shape": [len(focus)] + yx if yx and len(focus) > 1 else ([1] + yx if yx else None),
                     "channels": [], "origin": "2D image: the morphology image spatialdata_io.xenium loads"})
        chans = []
    if planes and len(planes) > 1:
        use = sel or None
        n = len(use) if use else len(planes)
        if use:
            origin = f"z-stack: Analysis.zstack_image(z={use}) of {source}, built from the chosen planes only"
        else:
            origin = f"z-stack of {source}, built with Analysis.zstack_image from the planes chosen in z_planes"
        rows.append({"name": stack_name or "zstack", "dims": 3, "axes": "czyx",
                     "shape": ([1, n] + yx if yx else None), "channels": chans, "planes": list(planes),
                     "origin": origin})
    return rows


def zarr_images(z: dict, own) -> list[dict]:
    """The images of an input zarr (SSAF-T outputs of an earlier run left out)."""
    rows = []
    for k, v in z["images"].items():
        if k in own:
            continue
        rows.append({"name": k, "dims": 3 if v["is_3d"] else 2, "axes": "".join(map(str, v["axes"])),
                     "shape": v["shape"], "channels": v.get("channels") or [],
                     "origin": ("z-stack: this 3D image cut to the planes chosen in z_planes" if v["is_3d"]
                                else "2D image: this image element")})
    return rows


def _area(r):
    s = r.get("shape") or []
    return s[-1] * s[-2] if len(s) >= 2 else 0


def _planes(r):
    s = r.get("shape") or []
    return s[-3] if len(s) >= 3 else len(r.get("planes") or [])


def recommend(plat, dim, rows, mask_z=None, stain=None):
    """[SSAF-T] (image name, why) recommended for downsampling a dim-D map, or (None, why not).
    2D map: the platform's 2D image. 3D map: the z-stack with the most planes (several planes
    sample the anchors through the tissue depth), else the 2D image."""
    if dim == 3:
        stacks = [r for r in rows if r["dims"] == 3 and _planes(r) > 1]
        if stacks:
            best = max(stacks, key=lambda r: (_planes(r), _area(r)))
            return best["name"], (f"the z-stack with {_planes(best)} planes: one tissue mask per plane samples "
                                  f"the anchors through the tissue depth"
                                  + ("" if len(stacks) == 1 else " (the one with the most planes)"))
        two, why = recommend(plat, 2, rows, mask_z, stain)
        if two is None:
            return None, "there is no image"
        return two, (f"the data has no z-stack: the 2D image ({why}), placed at the middle of the transcripts' "
                     f"z range")
    cand = [r for r in rows if r["dims"] == 2]
    if not cand:
        return None, "there is no 2D image"
    names = [r["name"] for r in cand]
    pref, why = PREFERRED_2D.get(plat, ((), ""))
    pick = next((k for k in pref if k in names), None)
    if pick:
        return pick, why
    if plat == "merscope" and mask_z is not None:
        pick = next((k for k in names if k.endswith(f"_z{mask_z}")), None)
        if pick:
            return pick, (f"the z{mask_z} mosaic (merscope_mask_z, the layer spatialdata_io loads); its "
                          f"{stain} channel makes the mask (merscope_stain)")
    if len(cand) == 1:
        return names[0], "the only 2D image"
    best = max(cand, key=_area)
    return best["name"], "the largest 2D image"


def choose(asked, dim, rows, rec):
    """(image used, how it was chosen). A 2D map takes a 2D image; a 3D map a 2D image or a z-stack."""
    if asked is None or not str(asked).strip():
        return rec, "image_key blank: the recommendation"
    asked = str(asked).strip()
    r = next((r for r in rows if r["name"] == asked), None)
    if r is None:
        offered = ", ".join(x["name"] + " (" + MODE[x["dims"]] + ")" for x in rows)
        raise ImageError(f"image_key {asked!r} is not one of the images this data offers: {offered}")
    if dim == 2 and r["dims"] == 3:
        other = [x["name"] for x in rows if x["dims"] == 2]
        raise ImageError(f"image_key {asked!r} is a z-stack, but this map is 2D (dimension 2), which "
                         f"downsamples on a 2D image only. Set image_key to one of the 2D images "
                         f"{other or '(none)'}, or make the map 3D (dimension: 3)")
    return asked, "image_key"


def lines(rows, used, dim, rec, why, how) -> list[str]:
    """The `images :` block of the dry_run plan."""
    if dim == 3:
        kinds = ("a 2D image (one tissue mask, placed at the middle of the transcripts' z range) or a "
                 "z-stack (one mask per plane, from the planes chosen in z_planes)")
    else:
        kinds = "a 2D image (one tissue mask); a 2D map takes a 2D image only"
    out = [f"images     : the image downsampling uses (image_key), chosen apart from the map's dimension: "
           f"{kinds}; recommend one and ask the user"]
    w = max([len(r["name"]) for r in rows] + [4])
    for r in rows:
        mark = "*" if r["name"] == used else " "
        off = "  [not for a 2D map]" if dim == 2 and r["dims"] == 3 else ""
        ch = f"; channels {', '.join(map(str, r['channels']))}" if r.get("channels") else ""
        pl = f"; planes {r['planes']}" if r.get("planes") else ""
        out.append(f"  {mark} {r['name']:<{w}}  {MODE[r['dims']]:<8} {r['axes']} {r.get('shape') or ''}{ch}{pl}  "
                   f"({r['origin']}){off}")
    out.append(f"  recommended: {rec} -- {why}")
    out.append(f"  this run   : {used} ({how}); marked *")
    return out


# ============================================================================ run time (stages)
def keep_planes(el, planes):
    """Only these planes (positions along z) of a 3D image element, each at the z it had.

    New plane j is old plane planes[0] + j * stride. In intrinsic coordinates (a pixel's centre at
    index + 0.5) that is u_old = stride * u_new + planes[0] + 0.5 - 0.5 * stride, put in front of
    the element's own transformations, so its x/y placement and z mapping are kept as they were."""
    from ssaft_core import scale0
    from spatialdata.models import Image3DModel
    from spatialdata.transformations import Scale, Sequence, Translation, get_transformation
    node = scale0(el)
    stride = planes[1] - planes[0] if len(planes) > 1 else 1
    data = node.data[:, planes, :, :]
    shift = Sequence([Scale([float(stride)], axes=("z",)),
                      Translation([planes[0] + 0.5 - 0.5 * stride], axes=("z",))])
    tfs = {cs: Sequence([shift, t]) for cs, t in get_transformation(el, get_all=True).items()}
    kw = {}
    if "c" in node.coords:
        kw["c_coords"] = [str(c) for c in node.coords["c"].values]
    out = Image3DModel.parse(data, dims=("c", "z", "y", "x"), scale_factors=[{"x": 2, "y": 2, "z": 1}] * 3,
                             chunks=(1, 1, 2048, 2048), transformations=tfs, **kw)
    srcs = (getattr(el, "attrs", None) or {}).get("ssaft_plane_sources")
    if srcs is not None and len(srcs) == node.sizes["z"]:
        out.attrs["ssaft_plane_sources"] = [srcs[k] for k in planes]
    return out
