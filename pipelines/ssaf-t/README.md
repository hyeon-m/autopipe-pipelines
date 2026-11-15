# ssaf-t

Segmentation-free spatial cell-type mapping with SSAF-T: transcripts (or count bins) become a
kernel-density vector field, anchors sampled from the tissue are clustered and annotated by CeLLMan,
and every voxel of the section is assigned a cell type; optionally the map is cut into cells seeded
by nuclei. The pipeline runs the SSAF-T workflow stage by stage, with the CeLLMan annotation done by
Claude Code (`claude -p`) or Codex (`codex exec`) headless. One sample, or a batch of samples with
shared settings.

Platforms: Xenium, MERSCOPE, CosMx (experimental), Visium HD, Stereo-seq, as raw output folders or
as SpatialData `.zarr`.

## Inputs

| key | what |
|---|---|
| `input` | one sample: dataset under `/input`: Xenium / MERSCOPE / CosMx output folder, Visium HD `outs` folder, Stereo-seq SAW output (the numbered layout, or every file in one folder), or a SpatialData `.zarr` |
| `samples` / `manifest` | a batch: the samples and what differs per sample (see Batch) |
| `species`, `organ` | Cell-Ontology module, e.g. `mouse` / `kidney` (an SSAF-T count zarr may record them) |
| `dimension`, `image_key` | the map (2D / 3D) and the image downsampling uses (a 2D image or a z-stack), chosen apart (see below) |
| `reference_h5ad` + `reference_column`, or `reference_celltypes` | optional: limit the annotation to given cell types (reference mode); chat setup asks whether to use one |
| `nuclei_file` | optional nucleus polygons for segmentation (needed on MERSCOPE and Visium HD) |

A raw folder is loaded by the load stage (`scripts/ssaft_load.py`) into the run's SpatialData zarr,
and the stages work on that: with spatialdata_io's `xenium`, `merscope` (plus the pipeline's own
reading of the transcripts and mosaics), `visium_hd` and `stereoseq` readers, and for CosMx with the
slide loader, which reads the flat files directly and joins the FOVs into one slide. A SpatialData zarr is copied in. dry_run names the loader
(`loader :`).

Stereo-seq takes SAW output only:

- the numbered SAW layout (`03.register/`, `04.tissuecut/`, `041.cellcut/`, `051.cellcluster/`) is
  read with spatialdata_io.stereoseq;
- a SAW output with every file in one folder (`<chip>.tissue.gem.gz`, `<chip>_regist.tif`,
  `<chip>.cellbin.gef`, `<chip>.tissue.gef`, ...), which spatialdata_io does not read, is
  read by the pipeline's flat SAW loader (`scripts/saw_flat.py`): the DNBs of the tissue GEM
  summed into bins of `vf_bin_um` (blank: the GEM's own bin, 0.5 µm), each placed at its centre in
  the chip's µm frame; the registered image (1 px = 1 DNB = 0.5 µm) as the staining image; and, added
  to the script, the cells of the cellbin.gef as `shapes['cells_polygons']` (each cell's border
  offsets added to its centre, as spatialdata_io builds them, in µm), which segmentation seeds on.
  Where the script differs, the data decides: the stain polarity follows the image (H&E: dark tissue
  on bright; ssDNA: bright on dark, where the script, made for an H&E chip, always takes H&E), and
  with several registered images the ssDNA one is used;
- a folder with a GEM table but no SAW `.gef` file is not read (dry_run says so).

The input may be a link (e.g. made by AutoPipe's prepare_input in the workspace input folder) to a
dataset anywhere on the server: AutoPipe 0.0.36 or later mounts each link's target at its real path,
for dry_run and for the run.

A zarr's platform is recognised from its element names (e.g. `morphology_focus` / `nucleus_boundaries`
= Xenium, `<id>_transcripts` + `<id>_z<N>` = MERSCOPE, `<fov>_points` = CosMx, `square_###um` /
`*_full_image` = Visium HD, `*_regist` / `bin1_table` = Stereo-seq, `uns['cellman']` = SSAF-T count
zarr), the reader attribute only breaking ties; `platform` overrides it. dry_run shows the platform
detected and its candidates, to be confirmed with the user.

## Checks at dry_run

- **The files the load opens** (raw folders, `scripts/reader_files.py`): taken from the readers as
  installed in the image (spatialdata_io 0.7.0) called the way the load calls them, and from what the
  load (or the flat SAW loader) reads itself. A missing file stops the dry_run with its name and what
  it is for (`files :` lists them when all are there):

  | platform | files |
  |---|---|
  | Xenium | `experiment.xenium`, `cells.zarr.zip`, `cells.parquet`, `cell_feature_matrix.h5`, `nucleus_boundaries.parquet` (the reader reads it in every run, segmentation or not), `transcripts.parquet`, the morphology focus image (`morphology_focus.ome.tif` before Xenium Analyzer 2.0, else `morphology_focus/morphology_focus_0000.ome.tif` or `ch0000_dapi.ome.tif`); downsampling on the z-stack: `morphology.ome.tif` |
  | MERSCOPE | `detected_transcripts.csv`, `images/micron_to_mosaic_pixel_transform.csv`, `images/mosaic_<stain>_z<N>.tif` of every stain loaded at the mask layer; downsampling on the z-stack: the `merscope_stain` mosaic of every plane chosen. The cell table (`cell_by_gene.csv`, `cell_metadata.csv`) is optional |
  | CosMx | `*tx_file.csv`, `*fov_positions_file.csv`, `CellComposite/` images; with segmentation `*metadata_file.csv` |
  | Visium HD | `[<id>_]feature_slice.h5`, `binned_outputs/square_<NNN>um/` `filtered_feature_bc_matrix.h5`, `spatial/scalefactors_json.json`, `spatial/tissue_positions.parquet` (the bin size chosen) |
  | Stereo-seq, numbered SAW | `03.register/*_regist.tif`, `04.tissuecut/*.gef`, `041.cellcut/*.cellbin.gef`, `051.cellcluster/*cell.cluster.h5ad` |
  | Stereo-seq, SAW in one folder | `<chip>.tissue.gem.gz` (else another `<chip>*.gem.gz`), `<chip>_regist.tif`; with segmentation `<chip>.cellbin.gef` (else `<chip>.adjusted.cellbin.gef`) |

  A folder that shows a platform's files but lacks the ones needed to recognise it (e.g. MERSCOPE
  mosaics without `detected_transcripts.csv`, a Stereo-seq GEM table without the SAW `.gef` files) is
  reported as that platform with the files missing.
- **CosMx pixel size** (`pixel size :`): the load takes the µm per pixel from the size of the
  CellComposite images (5472 x 3648 px: 0.18 µm/px, the earlier camera; 4256 x 4256 px: 0.12028
  µm/px, the newer one). dry_run reads the first image's size and shows the value and the camera; on
  any other size it stops and asks for `cosmx_px_um`. A `cosmx_px_um` value overrides the image size
  (passed to the load as `load.px_um_override` when it differs from what the size gives). The
  summary table notes on its `load` row the value this run uses.
- **Values that have to come from the user** (species and organ; the z plane spacing of a 3D map on
  CosMx `cosmx_z_plane_um`, or MERSCOPE `merscope_z_plane_um` when the mosaics' OME headers have no z; `z_planes` when a z-stack used lacks
  some planes): they do not stop the plan where they are found. The dry_run reads the data and prints
  the whole plan made from it (platform, loader, files, the map, the images and the z plane table,
  by plane number while the spacing is not set, bin size, grid, nuclei, the stage table), then stops
  with a numbered `still needed:` list, so that they can be decided from what the data shows. Chat
  setup starts with this dry_run.

## Stages

The plan (input inspection, run.json from `scripts/run_config.py`) is made when Snakemake reads the
workflow, so dry_run prints it; a run writes it to `.ssaft/run.json` and `00_check/`. Then:
`check_env` (AI CLI login, one short prompt with the chosen model, GPU/CPU and the server GPUs
recorded; runs first) → `load` → `zstack` (a 3D map on a z-stack) → `tile` (Xenium/MERSCOPE,
optional) → `panel_vf` (count platforms) → `downsample` → `cluster` → `annotate` (CeLLMan) → `map` →
`filter` → `segment` (optional) → `final` → `report`.

Each stage takes the run.json values it reads as Snakemake params: re-running with a changed value
re-runs that stage and everything after it; earlier stages are kept (e.g. `resolution` re-runs from `cluster`). No stage depends
on a job that can run again (run.json is not written by a rule, and `check_env` is not an input of
any stage): when such a job reruns, Snakemake does not examine the jobs that take its output as
`ancient()` input, so their changed params, and even their missing outputs, would go unnoticed.
Every stage uses `threads` cores, or, with `threads` blank, the cores the run is given.

A re-annotation keeps the previous one as `05_annotation/annotation_out.prev-<time>`. When the
annotate stage runs again with the same values and the same clustering (a run interrupted after, or
during, the annotation), the annotation already written is used and no new CeLLMan session starts
(`05_annotation/cellman/annotation_inputs.json` records what an annotation was made from).

## The map (2D / 3D) and the image downsampling uses

SSAF-T takes these apart: `downsample_vector(z=2|3, image_key=...)`, `z` the
dimension of the map, `image_key` the ONE image element the tissue mask is made from. The pipeline
does the same (`scripts/ssaft_images.py`):

- **`dimension`**: a 2D or a 3D cell-type map. Visium HD and Stereo-seq map in 2D only (counts on a
  2D bin lattice). Xenium, MERSCOPE and CosMx: blank = 3D, built on the transcripts' z; a 3D map
  needs no z-stack image.
- **`image_key`**: the image downsampling uses, one of the candidates dry_run finds in the data,
  each marked **2D image** or **z-stack**:

  | input | 2D image | z-stack |
  |---|---|---|
  | Xenium | `morphology_focus` | `morphology_zstack` (`morphology.ome.tif`, or `zstack_source`) |
  | MERSCOPE | `<folder>_z<N>` (the `merscope_mask_z` mosaic, its `merscope_stain` channel) | `<stain>_zstack` (that stain's mosaics of every plane) |
  | CosMx | `composite` | none: only the cell-type map is 3D |
  | zarr | its 2D image elements | its 3D image elements (or `zstack_source`) |

  The z-stack names follow SSAF-T's `zstack_image` keys. A 2D map takes a 2D image only. A 3D map
  takes either: a **z-stack** gives one tissue mask per plane, from the planes chosen in `z_planes`
  (below), built with `Analysis.zstack_image(z=...)`; a **2D image** gives one mask, which SSAF-T
  places at the middle of the transcripts' z range, so no z-stack step and no planes. That makes a 3D
  map possible on CosMx flat files (`composite`) and on any data without a z-stack.
- **Recommendation** (blank `image_key`): a 2D map, the platform's 2D image; a 3D map, the z-stack
  with the most planes when there is one with several, else the 2D image. dry_run shows the
  candidates, the recommendation with its reason, and the one this run uses (`*`) in its `images :`
  block, kept in `00_check/images.json`.
- **The z plane spacing** of a 3D map on MERSCOPE (`merscope_z_plane_um`) or CosMx
  (`cosmx_z_plane_um`): their transcripts' z is a plane number and the µm between planes is not in
  the data, so a 3D map needs it whichever image is used; asked after the map and the image. MERSCOPE
  `z_step` follows it. Xenium's z is in µm already.
- **Tiling** (Xenium, MERSCOPE): blank = tiled when downsampling on a z-stack (every tile carries it),
  not tiled otherwise; a 3D map on a 2D image is never tiled.
- **Anchors** per 10 µm block of each imaging plane (`n_samples`): 10 for a 2D map, 5 for a 3D map
  (SSAF-T's defaults); Visium HD and Stereo-seq (`n_samples_per_block`): 10.

A 3D map on a z-stack still keeps a 2D image as `image.key`: the reference grid for tiling and pixel
coordinates; the summary table says so on its `image` row and names the z-stack.

**MERSCOPE element names** are `<folder>_z<N>` and `<folder>_transcripts`, after the data folder's
own name (e.g. `<data folder>_z3`), whatever path the folder is reached by in the container: a
link's target, or the folder mounted as `/input` itself (its name is then read from the container's
mount table). spatialdata_io names them `<parent>_<folder>_*`; they are renamed right after loading.

Visium HD and Stereo-seq use the count data's own image: raw folders the full-resolution H&E / the
registered ssDNA image (`*_regist`; a SAW output in one folder: `<chip>_regist.tif`, ssDNA or H&E),
an SSAF-T count zarr the staining image it names in `uns['cellman']` (`staining`). The tissue mask
uses the Otsu threshold as computed (`otsu_mult` blank), on every platform, the count platforms
included (SSAF-T 23dd889 takes `otsu_mult` there too), and on every plane of a z-stack.

## Sampling distance

`sampling_distance` is the grid spacing of the vector field and the cell-type map. It has a minimum,
and blank uses exactly that minimum; any value from the minimum up can be chosen, a smaller one stops
the dry_run:

- Visium HD: the bin size (`visiumhd_bin_um`); the counts are measured per bin;
- Stereo-seq: the bin size the counts are blurred on (`vf_bin_um`: the bin size an SSAF-T zarr
  records, else 0.5 µm: the bin1 table of a numbered SAW folder, the GEM's own bin of a SAW output in
  one folder), but at least 1 µm: max(bin size, 1); e.g. 2 µm for a zarr with 2 µm bins, 1 µm for a
  raw SAW folder. For a SAW output in one folder `vf_bin_um` is also the bin size the load sums the
  DNBs into (a multiple of 0.5 µm);
- Xenium, MERSCOPE, CosMx: 1 µm.

dry_run prints the value, the minimum and how it was derived on its `grid :` line.

## Choosing the z planes (a 3D map on a z-stack, Xenium and MERSCOPE)

A z-stack usually reaches past the section on both sides; a plane outside the tissue still gets a
tissue mask and a full plane of sampling, and its anchors carry ambient signal. The dry_run plan
therefore counts the molecules at each plane of the z-stack (ssaft's `Analysis.zstack_planes`; for a
z-stack image element of a zarr input, the same count over the element's planes) and prints, per
plane, its number, its z (once the plane spacing is known), the molecules, their fraction and the
cumulative fraction, plus the molecules on no plane and a **recommended range**: the shortest run of
planes holding at least 95% of the molecules on the planes.

`z_planes` keeps only some planes: a continuous range (`3-8`) or a range at a constant step (`2-10:2`),
because the z-stack places its planes at one constant spacing. The z-stack is built from those planes
only (`Analysis.zstack_image(z=...)`, which places each plane at its own number); a zarr's z-stack
element is cut down to them; every kept plane stays at the z it had. The count is cached in
`.ssaft/zstack_census.json` (reading a MERSCOPE CSV's z column takes minutes) and kept in
`00_check/zstack_planes.json`.

The planes chosen only decide on which planes' images the anchors are sampled (downsample). The
cell-type map is made over the whole tissue height the transcripts span, whichever planes are chosen.

`z_planes` may be left blank (every plane) when the z-stack has every plane from 0 (0, 1, ..., n-1),
on Xenium and MERSCOPE alike (for MERSCOPE every plane number is then passed to `zstack_image` as
`z=`). Given no `z=`, ssaft stacks the planes one `z_step` apart from 0 in list order, which is right
only for such a stack: a bundle shipping some planes only (z4 to z6) would have them placed at 0, 1,
2 × z_step instead of their own numbers, away from their molecules. For such a stack a blank
`z_planes` is listed under `still needed`, after the molecules per plane and the recommended range.
With a 2D image `z_planes` stays blank.

## Outputs (one sample: in the run's output directory)

One numbered folder per stage; a stage the platform does not run has no folder. Every folder has the
stage's `summary.json` (status, time, device, parameters, the stage's summary) and its log.

| folder | contents |
|---|---|
| `00_check/` | `env.json` (AI CLI check with the model that answered, GPU, the server GPUs the run was given, CPU), `gpu_info.txt`, `plan.txt` (the dry-run plan), `inspect.json`, `answers.json`, `run.json`, `images.json` (the candidate images, the recommended and the chosen one), `zstack_planes.json` (a z-stack, Xenium / MERSCOPE) |
| `01_load/` | load summary (the reader's elements; a SAW output in one folder: GEM rows, bins, genes, image, cells) |
| `02_zstack/` | z-stack summary (a 3D map on a z-stack; planes kept) |
| `02a_tile/` | `crops/` tiles (tiled Xenium / MERSCOPE) |
| `02_panel_vf/` | `panel.txt` gene panel (Visium HD, Stereo-seq) |
| `03_downsample/` | `anchors.h5ad` |
| `04_cluster/` | `<species>_<organ>_cls.h5ad`, `clusters_spatial.png` |
| `05_annotation/` | `annotation_out/` (CeLLMan: `*_annotated.h5ad`, `celltype_markers.yaml`, rationale, `consensus_report.md` where the judges disagreed), `cellman/` (prompt, session events, annotation summary, `annotation_inputs.json`) |
| `06_map/` | mapping summary |
| `07_filter/` | `celltype_map.png` (map preview with legend) |
| `08_segment/` | `ssaft_cells.h5ad` (when segmentation runs) |
| `09_final/` | `spatialdata.zarr`: the input elements (MERSCOPE: `<folder>_z<N>`, `<folder>_transcripts`; a SAW output in one folder: `tables['table']`, `images['staining']`, `shapes['cells_polygons']`) plus `labels['ssaft_celltype_map']`, `images['ssaft_celltype_map_rgb']`, `tables['ssaft_downsampled_vectors']`, and with segmentation `shapes['ssaft_cell_boundaries_2d']` (+3D), `tables['ssaft_cells']`; a tiled run on a z-stack also `images['ssaft_zstack']` (the chosen planes) |
| `ssaft_store/` | the SSAF-T zarr store (vector field, maps), shared by several stages |
| `report.md`, `report.json` | per unit: cell types (the `unknown1..N` groups counted apart), kept voxel fraction, cells, stage times, output paths; the device, with the server GPU number |
| `.ssaft/` | the run directory: `run.json`, `state.json`, `units/<unit>/` (links to the files above), unit logs, `devices.tsv`, `reference_labels.json` (a reference given as a list) |

A per-FOV CosMx run has one subfolder per FOV inside each stage folder.

**Anchor positions**: `obsm['spatial']` of `anchors.h5ad`, the cluster h5ad and the annotated h5ad
is what SSAF-T's `downsample_vector` returns. On Xenium, MERSCOPE and CosMx these are vector-field
grid indices, not µm (µm = `coordinate_offset` + index × `sampling_distance`, both in
`ssaft_store/`); they equal µm only when `sampling_distance` is 1 and the offset is 0. The anchors
in µm, in the frame of the images and the map, are `tables['ssaft_downsampled_vectors']` in
`09_final/spatialdata.zarr`. `ssaft_cells.h5ad` is in µm.

**The GPU**: inside the container the attached GPUs are numbered from 0 (`cuda:0`, nvidia-smi's
"GPU 0"). The server GPU that is comes from AutoPipe's run record (`.autopipe-run.json`, the `gpus`
the run was given), kept as `server_gpus` in `00_check/env.json` and shown on the `device` line of
`report.md`, e.g. `cuda:0 = server GPU 1, NVIDIA GeForce RTX 4090 (container GPU 0: ...)`.

## Batch (several samples, shared settings)

config.yaml holds the shared values and the sample list what differs per sample:

- every key of config.yaml is the value shared by all samples;
- `samples:` lists the samples; an item is a path, or a mapping of `name`, `input` and the config keys
  that differ for that sample (`nuclei`, `segment`, `organ`, `dimension`, `image_key`, `z_planes`,
  ...); `skip: true` leaves an item out. Alternatively `manifest:` names a TSV with the columns
  `order`, `name`, `path`, `platform`, `species`, `organ`, `dimension`,
  `condition`, `overrides` as a JSON object of config keys, `check`; a further column named after a
  config key overrides that key, and order `-` leaves a row out. `threads` and `gpu` are shared by
  the whole batch.

```yaml
samples:
  - name: patient_01
    input: /input/patient_01.zarr
    segment: false
  - name: patient_02
    input: /input/patient_02.zarr
    nuclei: nucleus_boundaries
```

At dry_run every sample is planned on its own (detection, defaults, validation) and its plan is
printed, with one line per sample (`ok` / `FAIL`). Then (`scripts/batch.py`):

- samples run strictly one after another, in list order, each with the single-sample workflow in
  `runs/<name>/` (the folders above, its own `report.md` and `pipeline.log`);
- each sample gets its own AI annotation (`replicates` judgments merged, 3 by default), and the same
  input can get different labels on different runs, so compare labels across samples with care; `reference_celltypes` restricts every sample to
  one list of expected cell types;
- a sample whose plan did not validate is `skipped`; a failed sample is recorded and the batch goes on
  (`batch_on_fail: stop` ends it there);
- a sample that failed in an earlier run is skipped unless `batch_retry_failed: true`;
- re-running the batch with a changed value re-runs each sample only from the stage the change
  affects (each sample is its own Snakemake run); with nothing changed nothing runs;
- `batch_prepare_only: true` stops every sample before the annotation (`prepared`, with the AI CLI
  checked); running again with `false` continues them from the annotation;
- a file named `STOP` in the output directory stops the batch before the next sample, and so does a
  disk with less than `batch_min_free_gb` free (default 200); running again resumes.

Batch outputs, next to `runs/`: `batch_report.md` / `batch_report.json` (per sample: status, platform,
hours, cell types, kept voxels, cells, link to its report, errors), `batch_state.json` (status, step, what it waits for,
start / end per sample), `batch_plan.json` and `batch_manifest.tsv` (the resolved samples and their
validation), `batch.log`.

## Running

Run through AutoPipe **0.0.36 or later** (it mounts the server user's Claude / Codex login into the
container, and the targets of links in the input folder; without a login there, `check_env` stops and
says so); input values are set in chat (`SETUP.md`). GPU is **optional**: the run uses the GPU
AutoPipe attaches when torch can use it (CUDA 13 stack, NVIDIA driver >= 580), else the CPU, which is
very slow on full sections. The device used is recorded.

Without AutoPipe: `docker build -t ssaf-t .`, then run `snakemake --cores 16` in `/pipeline` with the
data mounted at `/input` (read-only), an output directory at `/output` as the working directory
(`--directory /output`), and the server user's `~/.claude` (or `~/.codex`) mounted with
`CLAUDE_CONFIG_DIR` (`CODEX_HOME`) pointing at it. There is no `.autopipe-run.json` then, so the report
names the GPU by its container number only.

## The annotation step (AI CLI)

- It uses the Claude or ChatGPT login of the server user (no API key; never put one in the
  pipeline) and that subscription's quota; a batch uses it once per sample, and `replicates` (3 by
  default) multiplies the use of each annotation.
- The model is `llm_model` (chat setup proposes the model of the AI doing the setup, and the user
  confirms it), or the CLI's default when blank. When the run starts, `check_env` sends the CLI one
  short prompt with that model on the server user's account; if the account cannot use it (an unknown
  name, a model the plan does not include, no quota), the run stops there with the CLI's reason,
  before any stage. The model that answered is recorded as `model_used` in `00_check/env.json`, and
  the model the annotation session ran on in `05_annotation/summary.json` (`model_used`, plus
  `models_called` when a subagent used another model). Claude Code reports the model in its JSON
  output; Codex does not, so for Codex `model_used` is `llm_model` or the `model` of the mounted
  `config.toml`, marked `reported_by_cli: false`.
- It sends the anchors' expression profiles (gene symbols and values) and the tissue description to
  the AI provider, and CeLLMan's helper reads the Cell-Ontology server (cellman.pnucolab.com); the
  gene panel step (count platforms) reads the same server.
- Its labels are not exactly reproducible between runs. The CLI, its version and the model are
  recorded in `05_annotation/summary.json`.
- The prompt (`annotate.py`) runs CeLLMan fully autonomous: no questions are asked. The session's
  working directory is `.ssaft/`, prepared with CeLLMan's skill and agent; AutoPipe hides the server
  user's own skills, agents and settings in the run.
- The pipeline runs CeLLMan's standard annotation only, finalized once: the interactive free annotation
  is not part of it, and its outputs fail the stage. Without a reference, `replicates` independent
  judges (annotate-expert) are merged by an annotate-lead; a reference run is one annotate-expert.
- Without a reference the annotation is open (the dry_run's `annotation :` line says "no
  reference"); chat setup asks whether to use one (default none). Reference mode (CeLLMan annotate-reference-mode),
  when the user wants to limit the annotation to given cell types: an annotated h5ad
  (`reference_h5ad` + `reference_column`), or a list of cell types (`reference_celltypes`: inline, or a
  .json / .csv / .tsv label list under `/input`; an inline list is written to
  `.ssaft/reference_labels.json` and the prompt's `- reference :` line points CeLLMan at it).
- An annotation is made again only when a value it reads (species, organ, condition, the annotation
  keys, including `llm_model`) or the clustering changed. To force a new one with unchanged values,
  delete `05_annotation/cellman/annotation_inputs.json` and re-run the annotate stage.

## Configuration

Every key of `config.yaml` is commented, and chat setup explains and sets every key that applies to
the run (SETUP.md). Keys marked `[detail]` are analysis details, asked last, each with its default (`image_channel`, `revert_image`, `otsu_mult`, `n_samples`, `bandwidth`, `kernel`,
`n_samples_per_block`, `vf_bandwidth_um`, `vf_kernel`, `hvg_n`, `resolution`, `n_neighbors`,
`map_method`, `norm_f`, `show_unknown`, `tile`, `tile_size`, `memory_gate`, `codex_sandbox`,
`annotate_note`, `nuclei_frame`, `seed_rule`, `merscope_mask_z`).

Values the pipeline decides from the data or the map (the others: `scripts/platform_defaults.json`):

| key | blank means | checked |
|---|---|---|
| `dimension` | the map: 2D for Visium HD / Stereo-seq; 3D for Xenium / MERSCOPE / CosMx (on any image) | dry_run: Visium HD / Stereo-seq 2D only; a zarr's transcripts have z for 3D |
| `image_key` | the recommended image: a 2D map, the 2D image; a 3D map, the z-stack with the most planes, else the 2D image | dry_run: one of the candidates; a 2D map takes a 2D image only |
| `threads` | the cores the run is given | |
| `cosmx_px_um` | from the FOV image size as the load decides it (5472 x 3648 px: 0.18; 4256 x 4256 px: 0.12028); a value overrides the image size; a zarr input: 0.12028 | dry_run: shows the value and the camera; stops on an unknown image size |
| `cosmx_z_plane_um` | required for a 3D map on CosMx | dry_run: listed under `still needed` after the plan |
| `merscope_stain` | DAPI (every stain is loaded) | dry_run: one of the stains the data has |
| `merscope_mask_z` | layer 3 (spatialdata_io's `z_layers`), the MERSCOPE 2D image | dry_run: a mosaic of `merscope_stain` exists at that layer; else the layers there are are listed |
| `merscope_z_plane_um` | required for a 3D map on MERSCOPE (either image) | dry_run: listed under `still needed` after the plan (z planes by number) |
| `zstack_source` | Xenium `morphology.ome.tif`; a zarr's own 3D images | a z-stack candidate of `image_key` |
| `z_planes` (a z-stack, Xenium / MERSCOPE) | every plane, allowed only when the z-stack has every plane from 0; a 2D image: blank | dry_run: a continuous range or constant step, of planes the z-stack has; blank with some planes only: listed under `still needed`; set with a 2D image: stops |
| `tile` | tiled when downsampling on a z-stack, not tiled otherwise | dry_run: a 3D map on a 2D image is not tiled |
| `otsu_mult` | the Otsu threshold as computed (also Xenium 3D, and the count platforms) | |
| `n_samples` (Xenium, MERSCOPE, CosMx) | anchors per 10 µm block of each imaging plane: 10 for a 2D map, 5 for a 3D map (SSAF-T's defaults) | |
| `n_samples_per_block` (Visium HD, Stereo-seq) | anchors per 10 µm block: 10 | |
| `visiumhd_bin_um` | 2 µm (choices 2, 8, 16; an SSAF-T count zarr: the bin size it records) | dry_run: the dataset has that `square_###um` bin size |
| `sampling_distance` | the minimum: Visium HD the bin size; Stereo-seq max(vf_bin_um, 1) (2 µm for a 2 µm-bin zarr, 1 µm for a raw SAW folder); Xenium, MERSCOPE, CosMx 1 µm | dry_run: not below the minimum |
| `vf_bin_um` | Visium HD: the bin size; Stereo-seq: the SSAF-T zarr's recorded bin size, else 0.5 µm (a numbered SAW folder's bin1 table; a SAW output in one folder: the GEM's own bin, which the load also sums the DNBs into) | dry_run (SAW in one folder): a multiple of 0.5 µm, at least the GEM's bin |
| `z_step` (MERSCOPE) | the MERSCOPE z plane spacing | |
| `segment` (Stereo-seq) | on, seeded by `cells_polygons` (numbered SAW: 041.cellcut; SAW in one folder: `<chip>.cellbin.gef`) | dry_run: the cellbin.gef is there |
| `reference_celltypes` | no reference | dry_run: not together with `reference_h5ad` |
| `llm_model` | the CLI's default model (its name recorded) | check_env: the CLI answers with it on the server user's account |

## Implementation (`scripts/`)

- `plan.py` (at parse time, so dry_run prints it): inspects the input (`inspect_raw.py` for raw
  folders; a zarr from its metadata and element names), turns config.yaml into answers and builds
  run.json with `run_config.py` over `platform_defaults.json`; checks the files the load opens and the
  CosMx pixel size (`reader_files.py`), the images (`ssaft_images.py`) and the z planes (`zplanes.py`),
  and shows the whole plan before the values still needed are asked for.
- `ssaft_stage.py` runs one stage for every unit: `ssaft_load.py` (load) or `ssaft_steps.py` (zstack,
  tile, panel_vf, downsample, cluster, map, filter, segment), through `ssaft_core.run_one`, which tees
  the unit log and records the stage in `.ssaft/state.json`. After each stage its files are moved
  into the stage's numbered folder with a relative link left behind (`layout.py`).
- `annotate.py`: the CeLLMan caller session (prompt, workdir, `cellman-kg` MCP server, the wait for
  its end, reuse of an annotation made from the same values and clustering); `check_env.py` before
  any stage; `report.py` after.
- `saw_flat.py`: the Stereo-seq loader for a SAW output in one folder.
- `batch_plan.py` and `batch.py`: the batch, each sample a Snakemake run of the single-sample workflow.
- The memory and GPU gates (`ssaft_core.py`): a heavy stage waits until the host (and its card) has
  room; a full-section estimate above what the host can ever free starts without waiting;
  `memory_gate: false` turns the memory gate off.
- Downsampling takes the one image chosen as `downsample_vector(image_key=...)` (SSAF-T 23dd889):
  the 2D image, the z-stack (a tiled run: the `ssaft_zstack` element every tile carries), or, for a 3D
  map on a 2D image, that 2D image; `otsu_mult` on every platform.
- PCA (cluster) on the GPU: when the full PCA does not fit, an exact covariance PCA in row chunks,
  halving the chunk on out-of-memory. Each retry starts after the failed attempt's exception has been
  handled, so its GPU memory is released first.
- Workaround for SSAF-T 23dd889, to remove once SSAF-T is fixed:
  - the tiled (chunked) 3D downsample reads each tile without its transcripts before asking their z
    range, and its messages name a variable `images` it does not have (a NameError under verbose).
    The pipeline lets that z-range call read the tile's transcripts and gives the messages the name
    (`ssaft_core.analysis`).
- No Python bytecode caches are written into the pipeline folder (`PYTHONDONTWRITEBYTECODE`); the
  Snakefile removes any left there.

## Vendored code

| folder | source | commit | license |
|---|---|---|---|
| `vendor/SSAF-T` | https://github.com/pnucolab/SSAF-T | 23dd889c0b4e618313bd0d1c078a0a9f7bbec8e4 ("minor update") | BSD-3-Clause |
| `vendor/CeLLMan` | https://github.com/pnucolab/CeLLMan | 4495e309c7c735d90075ffd3fe488eea73ebbbe5 (the submodule commit of SSAF-T 23dd889) | BSD-3-Clause |

None of them is modified. SSAF-T and CeLLMan are public repositories of the PNU Computational
Biology Lab under BSD-3-Clause; their `LICENSE` files are kept in `vendor/`.
