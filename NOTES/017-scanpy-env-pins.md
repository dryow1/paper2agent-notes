# 017 — The scanpy-env freeze, as a rollback map

**Date:** 2026-10-01
**Machine:** Dell laptop, Omarchy Linux, 15 GiB RAM, 8 CPU cores, **no GPU**
**Interpreter:** `/home/kee/Work/paper2agent/Scanpy_Agent/scanpy-env/bin/python` — **Python 3.12.14**
(`sys.executable`, printed by the environment itself)
**Method:** read-only inventory. **Nothing installed, nothing downloaded, no new environment.**
**Result:** **177 packages recorded**, in this note and in `NOTES/017-scanpy-env-pins.txt`.

> **This freeze is the rollback map before any future pip into `scanpy-env`.**

## How it was taken, and why not with `pip`

The ticket asked for `python -m pip freeze`. That **does not work here**:

```
/home/kee/Work/paper2agent/Scanpy_Agent/scanpy-env/bin/python: No module named pip
```

`scanpy-env` was built by `uv`, which does not place `pip` in the virtual environment. Installing
pip to inventory the environment would have been a package install — the exact thing this ticket
forbids — so the freeze was taken two other ways instead, and cross-checked:

1. `uv pip freeze --python Scanpy_Agent/scanpy-env/bin/python` (uv 0.12.17, the tool that built it)
2. `importlib.metadata.distributions()` from the stdlib, inside that interpreter

Both report **177 distributions with identical versions**. The only differences are cosmetic — uv
normalizes names to PEP 503 form (`Authlib` → `authlib`, `jaraco.classes` → `jaraco-classes`,
`docstring_parser` → `docstring-parser`) while the stdlib reports raw metadata names. The list
below is uv's normalized output, because that is the form `uv pip install -r` expects.

No activation script was sourced. The interpreter was invoked by absolute path, the same way the
`harmonypy` check in note 016 was run.

## Three things this freeze does *not* let you rebuild

Recording these because a rollback map that is trusted too far is worse than none:

1. **`scanpy` is not a PyPI pin — it is a local path install.** uv reports it as
   `scanpy @ file:///home/kee/Work/paper2agent/Scanpy_Agent/repo/scanpy`, and the stdlib reports the
   version as `1.13.0a3.dev0+g0d5fd1623`. That `g0d5fd162` is a git hash, and it matches the pinned
   checkout from note 002 (`repo/scanpy` at HEAD `0d5fd16`). **Restoring from this freeze requires
   that checkout to be present at that commit.** It is not in git, and it is *not* in
   `tests/MANIFEST.md` either — so nothing currently detects if it changes.
2. **177 installed ≠ 122 pinned.** `tmp/build/requirements.txt` holds **122** package lines
   (note 002 called it 123). The gap is build-time and transitive packages the compile step did not
   record. So `requirements.txt` has never been proven to reproduce *this* set — which is Stage 5's
   open question, now with a number attached.
3. **This is a version list, not the bytes.** It pins no hashes. Two installs of
   `numpy==2.3.5` from different wheels would both satisfy it.

## What it is good for

Exactly one thing, and it is the thing that was wanted: **if a future ticket allows a `pip`/`uv`
install into `scanpy-env` — the `harmonypy` question left open in note 016 — this is the state to
compare against afterwards, and the list to reinstall from if the install drags in a version
change nobody asked for.** Take a second freeze after any install and diff it against this file.

## The freeze — 177 packages

```
aiofile==3.12.3
aiohappyeyeballs==2.7.1
aiohttp==3.14.3
aiosignal==1.4.0
anndata==0.13.4
annotated-types==0.8.0
anyio==4.15.1
array-api-compat==1.15.0
asttokens==3.0.2
attrs==26.1.0
authlib==1.8.0
beartype==0.22.9
beautifulsoup4==4.15.0
bleach==6.4.0
cachetools==7.2.0
caio==0.12.4
certifi==2026.7.22
cffi==2.1.1
charset-normalizer==3.5.1
click==8.5.0
cloudpickle==3.1.2
comm==0.2.3
contourpy==1.4.0
cryptography==50.0.1
cycler==0.12.1
cyclopts==4.25.3
debugpy==1.8.22
defusedxml==0.7.1
dnspython==2.8.0
docstring-parser==0.18.0
donfig==0.8.1.post1
email-validator==2.3.0
entrypoints==0.4
exceptiongroup==1.3.1
executing==2.2.1
fast-array-utils==1.5.1
fastjsonschema==2.22.2
fastmcp==4.0.3
fastmcp-slim==4.0.3
fonttools==4.65.0
formulaic==1.2.2
frozenlist==1.8.0
google-crc32c==1.8.0
griffelib==2.3.0
h11==0.16.0
h5py==3.16.0
httpcore2==2.13.0
httpx2==2.13.0
idna==3.20
igraph==1.0.0
imageio==2.37.4
iniconfig==2.3.0
interface-meta==2.0.1
ipykernel==7.3.0
ipython==9.17.1
ipython-pygments-lexers==1.1.1
jaraco-classes==3.4.0
jaraco-context==6.1.2
jaraco-functools==4.6.0
jedi==0.20.0
jeepney==0.9.0
jinja2==3.1.6
joblib==1.6.0
joserfc==1.7.5
jsonref==1.1.0
jsonschema==4.26.0
jsonschema-path==0.5.0
jsonschema-specifications==2025.9.1
jupyter-client==8.10.0
jupyter-core==5.9.1
jupyterlab-pygments==0.3.0
jupytext==1.19.5
keyring==25.7.0
kiwisolver==1.5.1
lazy-loader==0.6
legacy-api-wrap==1.5
leidenalg==0.12.0
llvmlite==0.49.0
markdown-it-py==4.2.0
markupsafe==3.0.3
matplotlib==3.11.2
matplotlib-inline==0.2.2
mcp==2.2.0
mcp-types==2.2.0
mdit-py-plugins==0.6.1
mdurl==0.1.2
mistune==3.3.4
more-itertools==11.1.0
msgspec==0.21.1
multidict==6.9.1
narwhals==2.26.0
natsort==8.4.0
nbclient==0.11.0
nbconvert==7.17.1
nbformat==5.11.1
nest-asyncio2==1.7.2
networkx==3.7
numba==0.67.0
numcodecs==0.17.0
numpy==2.5.3
openapi-pydantic==0.5.1
opentelemetry-api==1.44.0
packaging==26.3
pandas==3.0.6
pandocfilters==1.5.1
papermill==2.7.0
parso==0.8.7
pathable==0.6.0
patsy==1.0.3
pexpect==4.9.0
pillow==12.3.0
platformdirs==4.11.12
pluggy==1.6.0
pooch==1.9.0
prompt-toolkit==3.0.53
propcache==0.5.4
psutil==7.2.2
ptyprocess==0.7.0
pure-eval==0.2.4
py-key-value-aio==0.4.6
pycparser==3.0
pydantic==2.13.5
pydantic-core==2.46.5
pydantic-settings==2.15.0
pygments==2.21.0
pyjwt==2.14.0
pynndescent==0.6.0
pyparsing==3.3.3
pyperclip==1.11.0
pytest==9.1.1
pytest-asyncio==1.4.0
python-dateutil==2.9.0.post0
python-dotenv==1.2.3
python-multipart==0.0.32
pyyaml==6.0.3
pyzmq==27.2.0
referencing==0.37.0
requests==2.34.2
rich==15.0.0
rich-rst==2.1.0
rpds-py==2026.6.3
scanpy @ file:///home/kee/Work/paper2agent/Scanpy_Agent/repo/scanpy
scikit-image==0.26.0
scikit-learn==1.9.1
scipy==1.18.1
scverse-misc==0.1.6
seaborn==0.13.2
secretstorage==3.5.0
session-info2==0.4.2
six==1.17.0
soupsieve==2.9.2
sse-starlette==3.4.11
stack-data==0.6.3
starlette==1.6.0
statsmodels==0.15.0
tenacity==9.1.4
texttable==1.7.0
threadpoolctl==3.7.0
tifffile==2026.9.20
tinycss2==1.5.1
tornado==6.5.10
tqdm==4.70.1
traitlets==5.16.1
truststore==0.10.4
typing-extensions==4.16.0
typing-inspection==0.4.4
umap-learn==0.5.12
uncalled-for==0.4.0
urllib3==2.8.0
uvicorn==0.53.0
watchfiles==1.3.0
wcwidth==0.8.4
webencodings==0.6.1
websockets==17.1
wrapt==2.4.1
yarl==1.25.1
zarr==3.4.0
```

Machine-readable copy: [`017-scanpy-env-pins.txt`](017-scanpy-env-pins.txt) — the same 177 lines,
no prose, suitable for `uv pip install -r` once the `repo/scanpy` caveat above is handled.

## Caveats

- **One moment in time.** Taken 2026-10-01, before any install. It describes the environment as it
  stands with `harmonypy` **absent** (note 016).
- **`pip freeze` itself was never run**, because pip is not installed. The two substitutes agree
  with each other, but neither was checked against pip's output.
- **Versions only, no hashes**, and `scanpy` is a path install — see above. This file alone cannot
  rebuild the environment.
- Nothing here was installed, upgraded, removed or downloaded. The environment is byte-identical to
  before this ticket.

## Next

Unchanged: **Stage 5** (`requirements.txt` has never been installed — and point 2 above now
quantifies the gap: 122 pins against 177 installed packages), the **TPM/CPM hole** in Scanpy QC, and
the Scanpy **ZIP + relocation test**. Harmony/paper 2 stays **refused** until a pip install into
this environment is explicitly allowed (note 016); when it is, this freeze is the before-picture.
