# Setting up an SSAF-T run (chat)

This pipeline runs the SSAF-T workflow, on one sample or on a batch of samples that
share their settings. Which values matter depends on the input (platform, the images the data has,
what a SpatialData zarr contains), so they are decided in conversation, in this order: the input,
then a dry_run that reads the data, then the questions, each answered from what the dry_run showed.
**Every key of config.yaml that applies to this run is explained to the user and set with them**:
none is left at its default without the user having seen it. Each question gives what the key does
(from its config.yaml comment, in plain words), the values it takes, its default, the value the
dry_run read from the data or recommends (if any) and when another value would be chosen. A value the
data decides (the CosMx pixel size, a zarr's bin size, the Stereo-seq blur raster) is shown with
where it comes from and confirmed, or overridden when the user wants. Keys that do not apply to this
run (another platform's keys, the z-stack keys on a 2D image, the batch keys for one sample) are
named in one line with the reason, not asked one by one.

## 1. Ask for the input first

Ask which dataset(s) to run. Accepted:

- a raw output folder: Xenium (`experiment.xenium`), MERSCOPE (`detected_transcripts.csv` +
  `images/`), CosMx flat files (`*_tx_file.csv`, `*fov_positions_file.csv`), Visium HD (the `outs`
  folder with `binned_outputs/square_002um`), Stereo-seq SAW output, either the numbered layout
  (`03.register/`, `04.tissuecut/`, ...) or every file in one folder (`<chip>.tissue.gem.gz`,
  `<chip>_regist.tif`, `<chip>.cellbin.gef`, `<chip>.tissue.gef`, ...). The load stage reads it with
  the platform's reader (spatialdata_io xenium, merscope, visium_hd, stereoseq for the numbered SAW
  layout; the pipeline's CosMx slide loader, which joins the FOVs into one slide; for a SAW output in one folder the pipeline's flat SAW
  loader, which sums the tissue GEM's DNBs into bins) into the run's
  SpatialData zarr; dry_run names it on its `loader :` line. Only SAW output is read for Stereo-seq:
  a folder with a GEM table but no SAW `.gef` file is refused;
- a SpatialData `.zarr` (any of these platforms, or an SSAF-T count zarr with `uns['cellman']`).

Stage it with prepare_input: it makes a link to the dataset in the workspace input folder; pass that
folder as `input_dir` and set `input` to the link's name. The link may point anywhere on the server:
AutoPipe 0.0.36 or later mounts each link's target at its real path in the container, for dry_run and
for the run, so the links prepare_input makes work as they are. For a batch, link every sample into
the same input folder.

## 2. Read the data with a first dry_run (before asking anything else)

- `list_files` on the dataset folder (and on a zarr's `images/`, `points/`, `shapes/`, `tables/`).
- Set only `input` (or the batch's `samples:`) and leave `species` and `organ` blank, then run
  `dry_run`. It reads the data and prints the whole plan made from it, then stops with
  `the plan made from the data, before the values it still needs are set` and a numbered `still
  needed:` list (here "species and organ are required"; for a 3D map on MERSCOPE or CosMx also the z
  plane spacing, and `z_planes` when a z-stack used lacks some planes): this stop is expected. The
  plan shows the detected platform and its candidates, the reader, the files the load opens (all
  present), the zarr elements or the raw folder's files and hints (stains and planes, FOVs, bin
  sizes, the Xenium panel name, full-resolution images, nucleus boundaries, a GEM header), the CosMx
  pixel size, the `map :` line (2D or 3D), the `images :` block (the images downsampling can use,
  each marked 2D image or z-stack, the recommendation and the one used) and, for a z-stack, the z
  plane table, the grid and bin size, the nuclei, GPUs and cores, and the summary table.
  Show it to the user. The items under `still needed` are asked in the order of section 3, each when
  its turn comes, not all at once.
- A dry_run that stops earlier, without the plan, says what is wrong in the data (a file the reader
  opens, a layout it cannot read) or in a value that is set; most of those stops are questions for
  the user below.
- Do not read raw data files. Small metadata (e.g. `experiment.xenium`, a zarr's `zarr.json`) is under
  the input directory, so reading it needs the user's consent; the dry_run plan usually makes it
  unnecessary.

## 3. Ask, in this order (from what the dry_run showed)

1. **Platform.** Show the detected platform and its candidates (the `platform :` line) and ask the user
   to confirm it. Set `platform` only if the user names another candidate.
2. **Segmentation** (`segment`, `nuclei`, `nuclei_file`): explain how it works on this platform before
   anything else is decided, and ask whether to use it and which nuclei. The map can be cut into cells seeded by nuclei; where the nuclei come from:
   - Xenium: `nucleus_boundaries` (shipped with the output) — segmentation on by default;
   - CosMx: `cell_centroids` (the vendor's cell centres) — on by default;
   - Stereo-seq: `cells_polygons`, the SAW cell segmentation (`041.cellcut/` of the numbered layout,
     `<chip>.cellbin.gef` of a SAW output in one folder) — on by default. A folder without the
     cellbin.gef stops at dry_run listing it: ask whether to run without segmentation (`segment:
     false`) or to add the file;
   - MERSCOPE and Visium HD: the data has no nuclei, so segmentation is off unless the user gives
     nucleus polygons as `nuclei_file` (parquet / geojson / csv with a WKT geometry column), which
     turns it on.
   In a zarr, the default element may be missing: dry_run then stops and lists its shapes and points;
   let the user pick one, give a `nuclei_file`, or set `segment: false` (map only). Never switch
   segmentation off without asking; without nuclei the result is the cell-type map only.
3. **`species`, `organ`** (**required, never guess**): propose them from what the dry_run showed (the
   panel or region name, the dataset or folder name, a zarr's `uns['cellman']`) and ask the user to
   confirm. `organ` must be a Cell-Ontology module name; the stop message lists the modules served.
4. **The map: 2D or 3D** (`dimension`). Ask this first, then the image (item 5); the two are chosen
   apart.
   - Visium HD and Stereo-seq: always a 2D map; do not ask. Tell the user in one line that these are
     count platforms (bins on a 2D lattice, no z).
   - Xenium, MERSCOPE, CosMx: ask **"a 2D map or a 3D map?"**. A 3D map assigns cell types through the
     tissue depth from the transcripts' z; it needs no z-stack image (item 5). Blank = 3D.
5. **The image downsampling uses** (`image_key`). Show the dry_run's `images :` block: the candidates
   found in the data, each marked **2D image** or **z-stack**, with the `recommended:` line and its
   reason, and ask which to use. Write the name the user picks as `image_key` (blank = the
   recommendation).
   - A **2D map** uses a 2D image only; the z-stack rows are marked `[not for a 2D map]`. Usually
     there is one 2D image (Xenium `morphology_focus`, MERSCOPE `<folder>_z3` in the `merscope_stain`
     channel, CosMx `composite`): say so and confirm it.
   - A **3D map** can use either:
     - a **z-stack** (Xenium `morphology_zstack` from `morphology.ome.tif`, MERSCOPE `<stain>_zstack`
       from the mosaics of every plane, or a zarr's 3D image): one tissue mask per plane. Recommended
       when the data has one with several planes. Then show the `z planes` table (ssaft
       `Analysis.zstack_planes`: per plane its number, its z in µm once the spacing is known,
       molecules, fraction and cumulative fraction, the molecules on no plane, and the recommended
       range: the shortest range with at least 95% of the molecules on the planes) and **ask which
       planes to use**: the recommended range, all planes, or another range. Only a **continuous
       range** (`3-8`) or **a range at a constant step** (`0-3:3`, `2-10:2`) can be used. Write
       `z_planes`; the z-stack is then built from those planes only (`Analysis.zstack_image(z=...)`).
       The planes decide only which planes' images the anchors are sampled on; the cell-type map still
       covers the whole tissue height the transcripts span. Run `dry_run` again: the chosen planes are
       marked `*`.
       - All planes: `z_planes` may stay blank when the z-stack has every plane from 0 (the table's
         `this run :` line then says "all N planes"). When it has only some planes (e.g. z4 to z6), a
         blank `z_planes` is listed under `still needed` with the table: the planes have to be named.
     - a **2D image**: one tissue mask, which SSAF-T places at the middle of the transcripts' z range.
       No z-stack step and **no question about planes**; the map is still 3D. Recommended when the
       data has no z-stack: CosMx flat files, a Xenium bundle without `morphology.ome.tif`, a zarr
       without a 3D image. A 3D map on a 2D image is not tiled.
   - `z_planes` set with a 2D image on a 3D map, or on a MERSCOPE 2D map, stops the dry_run: clear it,
     or pick the z-stack. On a Xenium or CosMx 2D map it is only noted and not used.
   - In a 3D map on a z-stack the summary table still shows a 2D image on its `image` row; the
     row says it is the 2D reference grid for tiling and pixel coordinates and names the z-stack used.
     On a 2D image the row says so.
6. **Platform-specific values.**
   - **The z plane spacing of a 3D map**: the transcripts' z is a plane number on MERSCOPE and CosMx.
     **Settle it after the map and the image (and the planes) are decided**, never before, and never
     for a 2D map.
     - **Xenium** reads its plane spacing from the data: show it (`z_step` below).
     - **MERSCOPE**: the dry_run reads it from the OME headers of the mosaics (the `z spacing :` line
       says so) and shows each plane's z in the `z planes` table: show the table, confirm the spacing,
       and ask which planes to use (item 5). If the headers have no z, the table is by plane number
       and `merscope_z_plane_um` is under `still needed`: offer the spacing of their experiment or the
       usual MERSCOPE 1.5 µm, say it is not read from their data, and write the one they choose.
       `merscope_z_plane_um` overrides the headers. `z_step` follows it.
     - **CosMx** has no z-stack here: a 3D map uses the 2D composite image, and only the cell-type map
       is 3D, so there are no planes to choose. `cosmx_z_plane_um` is not in the flat files: ask for
       the experiment's value and do not propose one of your own; without a value, offer a 2D map.
   - **MERSCOPE**:
     - `merscope_stain`: every stain in the folder is loaded; ask which one makes the tissue mask and
       the z-stack, offering only the stains the dataset has (the `stains` hint; e.g. DAPI, PolyT,
       Cellbound3). Blank = DAPI.
     - The 2D mosaic's z layer (`merscope_mask_z`): asked in the analysis details (item 10), with the
       layers the dataset has a mosaic of (default 3). If the dataset has no mosaic of the chosen
       stain at layer 3, dry_run stops and lists the layers there are: ask it then.
     - The elements are named after the data folder (`<folder>_z3`, `<folder>_transcripts`, the
       `element names` line), whatever path the folder is reached by.
   - **CosMx: the pixel size (µm per pixel) — show and confirm.** The load takes it from the FOV image
     size: CellComposite images of 5472 x 3648 px are 0.18 µm/px (the earlier camera), 4256 x 4256 px
     0.12028 µm/px (the newer one). Show the value and the camera from the `pixel size :` line and ask
     the user to confirm it (`cosmx_px_um` stays blank) or give the instrument's own value. When the
     dry_run stops because the images have another size (not a known camera), the value has to be
     given. A value in `cosmx_px_um` overrides the image size (the line then says so); a zarr input,
     which has no FOV image file, uses `cosmx_px_um` or 0.12028.
   - **CosMx: the FOVs** (`cosmx_fovs`), a per-FOV zarr (spatialdata_io.cosmx) only: show its FOVs
     and ask whether to run all of them (blank) or only some (comma-separated ids). A raw CosMx folder
     always joins every FOV into one slide: say so; the key does not apply.
   - **Other image sources** (`zstack_source`, Xenium, a zarr): say which z-stack the dry_run
     found (Xenium `morphology.ome.tif`, a zarr's 3D images) and ask whether another z-stack file
     under `/input` should be offered as well (blank = none). MERSCOPE, a CosMx folder: not used.
   - **The z-stack plane spacing** (`z_step`, a 3D map on a z-stack only): show the value the run
     uses (Xenium: the OME metadata; MERSCOPE: the spacing above) and confirm it (blank).
7. **Sampling distance** (`sampling_distance`, the grid spacing of the vector field and the cell-type
   map). Show the minimum (the `grid :` line), say that blank uses exactly that minimum and that a
   larger value gives a coarser, faster map; any value from the minimum up may be chosen, and a
   smaller one stops at dry_run. Ask which to use.
   - Visium HD: first `visiumhd_bin_um`, 2, 8 or 16 µm (default 2), offering only the sizes the
     dataset has (an SSAF-T count zarr records its bin size: show it and confirm). The minimum is that
     bin size. Then `visiumhd_fullres_image`: show the full-resolution H&E the dry_run found (blank =
     that one) and ask whether to use it, another file under `/input`, or `none` (Space Ranger images
     only).
   - Stereo-seq: the minimum is the bin size the counts are blurred on (`vf_bin_um`), but at least
     1 µm: max(bin size, 1). E.g. 2 µm for a zarr with 2 µm bins, 1 µm for a raw SAW folder.
   - Xenium, MERSCOPE, CosMx: 1 µm.
8. **The count platforms' image and blur raster** (Visium HD, Stereo-seq): show and confirm.
   - Visium HD: the blur raster `vf_bin_um` is the bin size (`visiumhd_bin_um`); keep it blank.
   - SSAF-T reads the count data's own image (an SSAF-T count zarr: the staining image it names; a
     SAW output in one folder: `<chip>_regist.tif`, ssDNA or H&E); the `mask image :` line says which.
   - Stereo-seq: the blur raster `vf_bin_um` and where it comes from (the `note :` line: the bin size
     an SSAF-T zarr records in `uns['cellman']['meta']['bin_um']`, else 0.5 µm, the bin1 table or the
     GEM's own bin). For a SAW output in one folder the load sums the GEM's DNBs into bins of this
     size (the `bin size :` line). Ask
     whether to keep it (blank) or set another (a multiple of 0.5 µm, at least the GEM's own bin).
9. **Annotation.**
   1. `annotate_host`: claude or codex. **Propose the CLI you (the AI in this chat) run in**: codex when
      you are Codex, claude when you are Claude Code, and ask the user to confirm. The `claude` in
      `config.yaml` is only the value for a run nobody set up; do not propose it because it is there.
   2. `llm_model`: **recommend your own model**, the model you (the AI in this chat) run as, by the
      name the chosen CLI accepts (e.g. `claude-opus-5-5`), and ask the user to confirm it. Write it
      when they confirm; write the model they name instead if they want another; leave it blank (the
      CLI's default model) if they prefer that. If your model cannot run on the chosen host (e.g. you
      are a Claude model and the host is codex), say so and ask for a model of that host, or blank.
   3. Tell the user that when the run starts, `check_env` calls the CLI once with this model on the
      server user's account and stops the run with the reason if the model cannot be used there;
      the model actually used is recorded (`00_check/env.json` and `05_annotation/summary.json`,
      `model_used`), also when `llm_model` is blank.
   4. `condition` (e.g. normal, tumor) and `disease` (an optional hint, e.g. glioblastoma; blank =
      none): ask both; propose from the dataset name when it says.
   5. **Reference** (CeLLMan annotate-reference-mode): explain that the annotation is open by default
      (any cell type of the organ module) and that a reference limits it to given cell types; ask
      whether to use one (default: none). If yes, take one of:
      - an annotated h5ad: `reference_h5ad` + `reference_column` (the label column);
      - a list of cell types: `reference_celltypes`, comma-separated in config (`"hepatocyte, Kupffer
        cell, cholangiocyte"`), or a label-list file under `/input` (.json / .csv / .tsv);
      and ask `reference_scope`: strict (default; labels only from the list) or refine (their ontology
      descendants too). Cell types with no Cell-Ontology term in the module end up Unknown.
   6. `replicates` (without a reference): explain that CeLLMan merges this many independent judgments
      into one annotation and reports where they disagreed (`consensus_report.md`); it uses about that
      many times the annotation quota. Blank = 3; 1 = a single judgment; an even number is raised by
      one. A reference run is always one judgment.
10. **Analysis details** (the keys marked `[detail]` in config.yaml). After items 1–9, go through
    every detail key that applies to this run, grouped by step as below. For each: what it does, its
    default (and the value it gives here, e.g. `n_samples` 5 for a 3D map), the values it takes and
    when another value helps. Ask per group: keep the defaults, or change which. The user may answer
    "defaults" for a whole group; record each value they set.
    - Tissue mask: every platform `otsu_mult` (e.g. 0.6 for a wider mask); Xenium, CosMx `image_channel`,
      `revert_image` (true for dark tissue on a bright background); MERSCOPE `revert_image` and
      `merscope_mask_z` (the layers with a mosaic of `merscope_stain`; its stain decides the channel).
      Visium HD and Stereo-seq read neither `image_channel` nor `revert_image`.
    - Tiling (Xenium, MERSCOPE): `tile`, `tile_size`.
    - Anchors and vector field: transcript platforms `n_samples`, `bandwidth`, `kernel`; count
      platforms `n_samples_per_block`, `vf_bandwidth_um`, `vf_kernel`, `hvg_n`.
    - Clustering: `resolution`, `n_neighbors`.
    - Mapping and filtering: `map_method`, `norm_f`, `show_unknown`.
    - Segmentation (when on): `nuclei_frame` (only with a `nuclei_file`); `seed_rule` (3D map only:
      majority, or all, which gives more cells than nuclei).
    - Annotation: `annotate_note`; `codex_sandbox` (codex host only).
    - Run: `memory_gate` (wait before a memory-heavy stage until the host has room; false = never
      wait).
11. **Compute** (section 5).

Keys that do not apply to the run (another platform's keys, z-stack keys on a 2D image, `z_step`
without a z-stack, `nuclei_frame` without a `nuclei_file`, `codex_sandbox` with claude, the batch
keys for one sample) are listed in one line with the reason. In the final summary (section 6) list
every key with the value set, marking the ones left at their default.

## 4. A batch (several samples)

Ask the shared decisions once (section 3 applies to all samples) and
put only what differs per sample into the sample list.

- `samples:` in config.yaml, one item per sample in run order: `name`, `input`, and the config keys
  that differ (typically `nuclei` / `segment`, sometimes `organ`, `dimension`, `image_key`,
  `z_planes`, `condition`, `merscope_stain`, `merscope_z_plane_um`, `cosmx_z_plane_um`,
  `visiumhd_bin_um`, `sampling_distance`). Or a `manifest:` TSV (columns order, name, path, ..., overrides) if
  the user already has one.
- The first dry_run (species / organ blank) shows every sample's plan; each sample then stops with
  its `still needed` list (`FAIL`), which is expected at that point.
- Leave `input` blank in a batch; `threads` and `gpu` cannot differ per sample.
- The map in a batch: ask 2D or 3D once per transcript platform present (Xenium / MERSCOPE / CosMx),
  never for Visium HD or Stereo-seq (always 2D, also when the shared `dimension` is blank). A shared
  `dimension: 3` would make the count-platform samples fail; set `dimension: 2` on their items.
- The image per sample: the dry_run prints each sample's `images :` block; image names differ between
  samples (a MERSCOPE `<folder>_z3`, a zarr's element names), so set `image_key` per sample where
  the recommendation is not wanted, and `z_planes` per sample for a z-stack. A sample whose z-stack
  lacks some planes is `FAIL` until `z_planes` is set.
- CosMx samples: leave `cosmx_px_um` blank (each sample's pixel size comes from its own images); set
  it per sample only for a sample whose camera is unknown.
- A shared `sampling_distance` must be at least every sample's minimum; leave it blank to use each
  sample's own minimum.
- Per-sample nuclei: dry_run lists each zarr's shapes and points; ask per sample, never switch a
  sample's segmentation off without asking.
- Ask once, each with its default: on a failure, continue with the next sample (default,
  `batch_on_fail: continue`) or stop the batch; `batch_retry_failed` (re-run samples that failed in an
  earlier run of this batch; default false); `batch_prepare_only` (stop every sample after the
  clustering, to run all annotations later; default false); `batch_min_free_gb` (stop before the next
  sample when the output disk has less free space; default 200, 0 = no check). Ask whether the
  samples come as a list (`samples:`) or an existing TSV (`manifest:`).
- Tell the user that each sample gets its own AI annotation (`replicates` judgments merged), so the labels of two samples are not
  directly comparable: the same input can get different labels on different runs. Anchors and maps
  are computed, the labels are judged. Ask whether a shared list of expected cell types
  (`reference_celltypes`) should restrict every sample to the same labels.
- Tell the user the run time and the annotation quota scale with the number of samples, and that the
  batch can be stopped before the next sample by creating `STOP` in the output directory (write_file),
  and resumed by running again with the same run_name.

## 5. Compute: ask, do not assume

AutoPipe reports the server's GPUs and CPU cores. Ask the user:

- **which GPU** to use (execute_pipeline `gpus`, e.g. `"1"`; prefer the least used), or none (`gpu:
  cpu`, very slow on full sections);
- **the device setting** (`gpu`): auto (default: the attached GPU if usable, else the CPU), cuda (stop
  without a usable GPU) or cpu;
- **how many cores** (execute_pipeline `cores`, or `cores` in config.yaml, which AutoPipe reads when
  the run starts; 0 = every core of the server), and `threads`, the cores per stage: blank (default)
  uses the cores the run is given (dry_run shows `cores per stage`); a number is capped at those.

`gpu: auto` uses the attached GPU when torch can use it (CUDA 13 stack, driver >= 580) and the CPU
otherwise; the device used is recorded in `00_check/env.json` and in every stage's `summary.json`.
Inside the container the attached GPUs are numbered from 0 (`cuda:0`, nvidia-smi's "GPU 0"); the
server GPU it is (the `gpus` given to execute_pipeline) is recorded in `00_check/env.json`
(`server_gpus`) and shown on the `device` line of `report.md` (e.g. "cuda:0 = server GPU 1").

## 6. Write config.yaml, dry_run, run

Write the chosen values, run `dry_run`, show the plan block to the user with a summary of every key
that applies (grouped as in section 3: input and tissue, the map and image, platform values, sampling,
annotation and reference, the analysis details, compute), each with its value and "default" where it
was kept, and one line for the keys that do not apply; get a go-ahead, then
execute. The annotation step sends the anchors' expression profiles to the AI provider through the
server user's subscription, and takes tens of minutes to hours per sample. It needs the server user's
Claude / Codex login inside the container, which AutoPipe 0.0.36 or later provides. config.yaml is
read when a run starts: do not change it while a run is still being planned, and give separate
datasets that need different settings their own runs or one batch.

A batch runs its samples one after another. To run more at the same time, start a second run with
another run_name on another GPU, and split the cores between the runs. Once a run has started it no
longer reads config.yaml (a batch keeps each sample's values in `batch_plan.json`), so write the second
run's values after the first run has started, and keep a copy of each run's config.yaml outside the
pipeline folder for later re-runs. Do not change `scripts/` or the Snakefile while a run is going. The
runs share the server user's AI subscription quota and the server's memory.

Results come in numbered stage folders (`00_check` ... `09_final`, see README.md); `report.md` sums
them up and `09_final/spatialdata.zarr` holds the final SpatialData. A batch puts each sample's folders
in `runs/<name>/` and sums the batch up in `batch_report.md`; `batch_state.json` shows the sample and
step running now. The cell type count in the reports leaves out the `unknown1..N` groups, which the map
keeps as separate classes (one per cluster of Unknown anchors); the report gives their number apart.

When the user asks how a run is going, read the end of `pipeline.log` (a batch: `batch.log` and
`waiting` in `batch_state.json`), not the GPU load, which may come from another process. A stage that
needs GPU or host memory waits for it without a time limit (`<stage> waits for 10 GiB free on cuda:0`):
it starts by itself when there is room. Tell the user it is waiting, for what, and what holds the GPU
(nvidia-smi); they can wait, free the GPU, or stop the run and start it again with the same run_name on
another GPU (finished stages are kept).

Propose a run_name and confirm it with the user: it names the output folder, and the same name resumes
that run.

If `build_image` says the image name already belongs to another pipeline, keep the folder and its name
and pass another `image_name` (for example `autopipe-ssaf-t-2`); use the same name for the runs. Do not
use `replace_image`, which takes the image away from the other pipeline, unless the user asks for it.

## Re-running with changed values

Use the SAME run_name, with that run's values in config.yaml. Each stage re-runs only when a value it
reads changed, and then everything after it: e.g. `norm_f` re-runs filter and
later only; `resolution` re-runs cluster, annotation and later; `organ` re-runs the annotation (count
platforms: from panel_vf); `sampling_distance` re-runs from panel_vf (Visium HD, Stereo-seq) or
downsample (transcript platforms); `image_key` re-runs from downsample (from load for a zarr, or
when it switches between a 2D image and a z-stack); `visiumhd_bin_um`, `merscope_stain`,
`merscope_mask_z`, `cosmx_px_um`, the z plane spacing, `dimension` (MERSCOPE, CosMx or a zarr; others
from the next stage that reads it) and, for a SAW output in one folder, `vf_bin_um` re-run from load; `z_planes` re-runs from zstack (a zarr's z-stack element: from
load); `otsu_mult` re-runs from downsample; `llm_model` re-runs check_env and the annotation; `replicates` re-runs the annotation and later;
`seed_rule` re-runs the segmentation only. Tell the
user beforehand when the annotation will run again. A previous annotation is kept as
`05_annotation/annotation_out.prev-<time>`. In a batch the same holds per sample; failed samples are
retried only with `batch_retry_failed: true`.

To re-run one sample of a batch with another value: same run_name and the same `samples:`, the changed
key on that sample's item, and `skip: true` on the other items. The skipped samples' folders are not
touched, and `batch_report.md` keeps their rows. A re-run replaces that sample's earlier map and cells;
to keep them, use a new run_name.

A run that was interrupted is resumed the same way, with unchanged values: finished stages are kept,
and when the CeLLMan annotation had already been written for the same values and clustering, it is
used again instead of starting a new annotation (`05_annotation/summary.json` says `reused`).

## dry_run errors

- `the plan made from the data, before the values it still needs are set: ... still needed: 1. ...`:
  the whole plan, then the values that have to come from the user, numbered. Expected on the first
  dry_run. Ask them in the order of section 3, each at its turn: species and organ (item 3);
  `z_planes` when a z-stack used lacks some planes (item 5: show the table and the recommended
  range); the z plane spacing of a 3D map on MERSCOPE / CosMx once the map and the image are decided
  (item 6). Set them and run dry_run again.
- `<reader> cannot load <folder>: these files it opens are missing: ...`: the folder lacks files the
  load stage opens (each listed with what it is for). The data has to be completed (a copy still in
  progress, a partial export); nothing in config.yaml replaces them, except where the line offers a
  setting (a missing Stereo-seq cellbin.gef: `segment: false`, asked first).
- `<folder> looks like a <platform> output folder (it has ...), but these files ... are missing`: the
  folder shows a platform's files but not the ones needed to recognise and load it; complete it.
- `<folder> holds Stereo-seq files (...) but is not a SAW output: a GEM table alone is not read`: ask
  for the SAW output folder (in one folder or the numbered layout) or a SpatialData zarr.
- `cosmx_px_um: the CellComposite images are W x H px, not a known CosMx camera`: ask the user for the
  instrument's µm per pixel and set `cosmx_px_um`.
- `input ... does not exist`: wrong path; the message lists what `/input` holds.
- `could not tell the platform`: set `platform` to one of the listed candidates.
- `no ... Stereo-seq layout recognised`: a raw Stereo-seq input must be a SAW output (the numbered
  layout, or every file in one folder with its `.gef` files) or a SpatialData zarr.
- `dimension ... is not possible for visiumhd / stereoseq`: these map in 2D only; leave `dimension`
  blank.
- `a 3D map needs the transcripts' z, and the zarr's transcripts ... have no z column`: a 2D map.
- `image_key ... is not one of the images this data offers`: pick from the listed images.
- `image_key ... is a z-stack, but this map is 2D`: pick a 2D image, or make the map 3D.
- `z_planes is set, which chooses the planes of a z-stack, but the image this 3D map downsamples on is
  the 2D image ...` (or `... but the map is 2D`): ask the user whether to use the z-stack on those
  planes (`image_key` the z-stack, `dimension: 3`) or the 2D image (clear `z_planes`).
- `tile: a tiled 3D map carries a z-stack in every tile`: leave `tile` blank with a 2D image.
- `image_key ... and zstack_source ... name different images`: keep only `image_key`.
- `zstack_source is set for a MERSCOPE folder`: leave it blank.
- `merscope_stain ... is not in this dataset` / `choose the stain`: pick one of the stains listed.
- `merscope_mask_z ...: ... has no mosaic_<stain>_z<N>.tif`: the message lists the z layers that have
  a mosaic; show them and set `merscope_mask_z` to one.
- `visiumhd_bin_um must be one of [2, 8, 16]` / `this dataset has no square_###um bins`: pick one of
  the bin sizes the message lists.
- `vf_bin_um ...: the load sums the GEM's DNBs ... into bins of vf_bin_um` (SAW output in one folder):
  a multiple of 0.5 µm, at least the GEM's own bin; or leave it blank.
- `nuclei ...: a flat SAW folder's cells are loaded as cells_polygons`: leave `nuclei` blank.
- `sampling_distance ... is too small`: the message gives the minimum (Visium HD: the bin size;
  Stereo-seq: max(vf_bin_um, 1); the others: 1 µm); set it to that or more, or leave it blank.
- `z_planes ... not evenly spaced` / `write a range`: use a continuous range (`3-8`) or a range at a
  constant step (`2-10:2`). `z_planes: no plane numbered ...`: the message lists the planes the
  z-stack has.
- `CosMx 3D on a zarr's z-stack (or a per-FOV zarr) is not supported`: use the zarr's 2D composite.
- `give the reference either as reference_h5ad ... or as reference_celltypes`: keep one of the two.
- `not a ... Cell-Ontology module`: pick from the listed modules.
- `there is no image to make the tissue mask from`: the zarr has no image; check the input.
- `segmentation needs nuclei`: pick a listed shapes/points element, a nuclei_file, or `segment: false`.
- `the configuration cannot run`: the plan's validation (z-step, reference column,
  nuclei file, ...); each line names the key to fix.
- batch: `FAIL [n] name: ...` is one of the above for that sample (it would be skipped); `not
  per-sample config keys` names a key that does not exist or cannot differ per sample; `no sample of
  the batch can run` when every sample fails (nothing would run, so the dry_run stops).
- At run time, `check_env` fails when the container holds no login for the chosen CLI (the message
  asks for AutoPipe 0.0.36 or later, which mounts it, and a logged-in server user); when the CLI cannot
  answer one short prompt with `llm_model` (`could not answer with the model ...`: the model name is
  wrong or the server user's account cannot use it; the message has the CLI's reason; fix the name,
  pick another model with the user, or leave it blank for the CLI's default); or when `gpu: cuda` has
  no usable GPU.
