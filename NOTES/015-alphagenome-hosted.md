# 015 — The hosted AlphaGenome MCP, probed over HTTP

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Method:** `curl` only. **No install, no model download, no clone, no new venv.** The existing
`src/` Paper2Agent checkout was not touched.
**Result:** **The server is live and the handshake works — and this laptop still cannot run a
single AlphaGenome prediction with it.** Both halves matter; the second is the point.

## Short version

Note 001 probed this space nine days ago and recorded "handshake OK, 22 tools". That reproduces
exactly. What note 001 flagged in a single line — the API key, and artifacts staying on the
server — is measured here, and it is worse than a footnote: **the two walls together mean the
hosted space can demonstrate tool *shape* and nothing else.** Every call either needs a key this
project does not have, or returns a file path that cannot be fetched.

## URLs tried and HTTP status

| # | Request | Status | What came back |
|---|---|---|---|
| 1 | `GET https://paper2agent-alphagenome-mcp.hf.space/` | **200** | 58 bytes of plain text: `MCP is on https://Paper2Agent-alphagenome-mcp.hf.space/mcp`. No UI, no docs. |
| 2 | `GET …/mcp` | **406** | Not Acceptable. The endpoint is POST-only and wants `Accept: application/json, text/event-stream`. |
| 3 | `POST …/mcp` → `initialize` | **200** | `text/event-stream`, `mcp-session-id: df0df134…`, `server: uvicorn`. |
| 4 | `POST …/mcp` → `notifications/initialized` | **202** | Accepted, no body. |
| 5 | `POST …/mcp` → `tools/list` | **200** | 33,534 bytes, **22 tools**. |
| 6 | `POST …/mcp` → `tools/call create_genomic_interval` | **200** | `isError: false`, plus a server-side CSV path. |
| 7 | `POST …/mcp` → `tools/call predict_dna_sequence` (no key) | **200** | `isError: **true**` — see below. |
| 8 | `GET …/data/tmp_outputs/genomic_interval_20261001_030204.csv` | **404** | The artifact from step 6 is not retrievable. |
| 9 | `GET …/file=/data/tmp_outputs/…csv` (Gradio-style) | **404** | Nor this way. |
| 10 | `GET …/health` | **200** | — |
| 11 | `GET …/sse` | **404** | No legacy SSE transport; streamable HTTP only. |
| 12 | `GET https://huggingface.co/spaces/Paper2Agent/alphagenome_mcp` | **200** | The canonical Space page, as advertised in the `link:` header. |

Note the trap in rows 6–7: **HTTP 200 does not mean the call worked.** MCP reports tool failures
inside a 200 body with `isError: true`. Anything that judges this server by status code alone
will record a refusal as a success.

## What the remote actually returned

**`initialize`** — the server identifies itself, and the version is worth recording:

```json
{"protocolVersion":"2025-06-18",
 "capabilities":{"experimental":{},"prompts":{"listChanged":true},
                 "resources":{"subscribe":false,"listChanged":true},
                 "tools":{"listChanged":true}},
 "serverInfo":{"name":"AlphaGenome","version":"1.13.1"}}
```

**`tools/list`** — 22 tools, identical to note 001's count:

`score_variants_batch`, `filter_variant_scores`, `create_genomic_interval`,
`create_genomic_variant`, `create_track_data`, `create_variant_scores`,
`genomic_interval_operations`, `variant_interval_operations`, `track_data_operations`,
`track_data_resolution_conversion`, `predict_dna_sequence`, `predict_genome_interval`,
`ism_analysis`, `explore_output_metadata`, `count_tracks_by_output_type`,
`visualize_variant_effects`, `visualize_gene_expression`, `visualize_chromatin_accessibility`,
`visualize_splicing_effects`, `visualize_histone_modifications`, `visualize_tf_binding`,
`visualize_contact_maps`

**The split that decides what is usable.** Reading the input schemas rather than the prose:

| | Count | Tools |
|---|---|---|
| `api_key` **required** | **3** | `predict_dna_sequence`, `predict_genome_interval`, `ism_analysis` |
| `api_key` optional | 10 | `score_variants_batch`, `explore_output_metadata`, `count_tracks_by_output_type`, and all 7 `visualize_*` |
| No `api_key` at all | 9 | the `create_*`, `*_operations`, `filter_variant_scores`, `track_data_resolution_conversion` |

**The three that require a key are the three that run the model.** The nine that need no key are
constructors and coordinate arithmetic — they build interval and variant objects, convert track
resolutions, filter score tables. Useful plumbing; no genomics.

**A working call** (`create_genomic_interval`, `chr11:116837600-116837700`, the same coordinates
note 001 used):

```
isError : false
message : Genomic interval created: chr11:116837600-116837700
artifact: /data/tmp_outputs/genomic_interval_20261001_030204.csv
```

**A blocked call** (`predict_dna_sequence`, `sequence: "ACGT"`, no key):

```
HTTP 200, isError: true
Input validation error: 'api_key' is a required property
```

The refusal is clean and arrives before any computation — the same quality the Scanpy tools were
praised for in note 003. Credit where due: the server validates rather than crashing.

## The artifact is stranded

Step 6 succeeded and returned `/data/tmp_outputs/genomic_interval_20261001_030204.csv`. That path
is **on the Hugging Face container**, and both plausible retrieval routes return **404** (rows
8–9). So the full extent of what a successful, key-free call delivers to this laptop is *a
sentence and a filename*. The CSV itself is unreachable.

This is the same design limit note 001 recorded for hosted Scanpy and TISSUE — those take file
paths on the server, so local data cannot be sent in. AlphaGenome's constructors take inline
arguments, so data can go *in*; it just cannot come back *out*.

## What this laptop cannot do without a local install

Stated plainly, because the live 200s make it tempting to overclaim:

1. **No prediction, of any kind.** The 3 model tools require a Google DeepMind AlphaGenome API
   key. This project has none, and obtaining one is outside the shop order. Nothing was attempted
   beyond the single unauthenticated call above, which was made to record the refusal.
2. **No results retrieved, even from calls that work.** Artifacts are 404. A hosted call cannot
   put a file on this disk.
3. **No GPU, and no model weights.** AlphaGenome is a large genomic model. Note 001 already
   classified it as "avoid for now: needs an API key and targets a huge genomic model", and
   nothing here changes that. Running it locally is not a `pip install` away — it is an API key
   *plus* a model this laptop was never going to host.
4. **Even a local MCP would not help.** Running `alphagenome_mcp` locally would fix the artifact
   404 — files would land on this disk — but the 3 prediction tools would still demand the key.
   Local installation solves the wrong wall.
5. **No conversion of the AlphaGenome paper was attempted**, and none should be: a `/paper2agent`
   run would produce tools that cannot be executed here.

So the honest summary of the hosted space, for this machine: **22 tools reachable, 9 runnable,
0 scientifically useful, 0 files retrievable.**

## What failed in this probe

Keeping the failures visible rather than folded into prose:

- `GET /mcp` → **406**. My first attempt used the wrong method; recorded rather than hidden.
- `GET /sse` → **404**. I guessed a legacy SSE endpoint. It does not exist; the server is
  streamable-HTTP only.
- Both artifact fetches → **404**. Two guesses at the retrieval URL, both wrong. There may be a
  route I did not find — this is a negative result from two attempts, not proof none exists.
- `predict_dna_sequence` → refused, as expected. Counted as a successful *probe* and a failed
  *call*; they are not the same thing.

## Caveats

- **One space, one day, one probe.** HF spaces sleep and restart; a 200 today is not an uptime
  claim. The `x-proxied-replica` header shows this is one replica behind a proxy.
- **I did not verify the server runs the real AlphaGenome.** `serverInfo` says
  `AlphaGenome 1.13.1`, which is a self-report. Without a key there is no way to check that the
  prediction tools do what they claim.
- **The 10 "optional key" tools were not called.** They may work keyless on pre-computed inputs,
  or may fail asking for a key once they reach the model. Untested, and I am not guessing which.
- No rate limits were probed, and nothing was called more than once.

## Next

Nothing here changes the Phase 1 open list: **Stage 5** (`requirements.txt` has never been
installed), the **TPM/CPM hole** in Scanpy QC, and the Scanpy **ZIP + relocation test**.

AlphaGenome stays where note 001 put it — **out of scope on this machine**, now with measurements
behind the judgement instead of a one-line guess. If a key ever appears, the first honest test is
`predict_dna_sequence` on a 4-base sequence against the hosted space, before anything local is
contemplated.
