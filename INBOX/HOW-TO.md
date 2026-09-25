# How to use this tray

## The short version

1. Copy **one** file into this folder.
2. Run the two commands below.
3. Read what it says.

That's it. If it can work with your file, it runs a quality check. If it can't, it tells you why
in one sentence and stops.

## Copy your file here

```bash
cp /path/to/your_data.h5ad ~/Work/paper2agent/INBOX/
```

**One file only.** Not two, not a folder with several inside. If there's more than one thing in
here, it refuses rather than guessing which you meant.

The file must be either an `.h5ad` file or a 10x Genomics `.h5` file. Those are the two formats
the tools can read.

## Run it

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env

# Look at the file without doing anything to it
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py --inspect-only

# Then actually run the quality check
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py
```

Start with `--inspect-only`. It reads the file and tells you how big it is, what's inside, and
whether it's the right kind of data — without computing anything. Nothing is changed or written.

Add `--mt-prefix mt-` if your data is from **mouse**. The default assumes human.

**Every run leaves a receipt** — a small JSON file recording what it decided and why. It goes to
`RESULTS/inbox-last.json` by default; use `--report somewhere/else.json` to put it elsewhere, or
`--no-report` to skip it. You get one whether it accepted your file or refused it.

## What actually happens

It looks at your file first, every time, and prints:

- how many cells and genes
- the range of the numbers inside
- whether the file already has other things stored in it
- **whether the numbers look like raw counts** — the thing that matters most

Then it either runs the quality check or refuses. The quality check measures each cell (how many
genes were detected, how much mitochondrial content) and drops obviously bad cells and genes that
barely appear. Results go into `INBOX_OUT/`: a new data file plus two figures.

**Then it stops.** It does not group your cells, does not decide what cell types they are, and
does not name anything. Those are separate decisions that you should make deliberately.

## When it refuses

| It says | It means |
|---|---|
| *INBOX is empty* | Nothing here. It will **not** download a dataset to keep busy. |
| *INBOX holds 2 items* | Take one out. |
| *'x' is a folder* | It needs a single file. A 10x folder of `matrix.mtx` + `barcodes` + `features` has to be converted to `.h5ad` first — that's a choice for you, not for a script. |
| *has suffix '.txt'* | Wrong format. `.h5ad` or 10x `.h5` only. |
| *over the 80 MB door limit* | Too big for this laptop. Take a subset, or use a bigger machine. |
| *would peak near … over the budget* | The file is small but the data inside is enormous. Same answer: subset it. |
| *X is not raw counts* | Your numbers have already been processed — logged, scaled, or rescaled. The quality check only means something on raw counts. If the file has the original counts tucked away inside, the message tells you exactly how to put them back. |

A refusal is not a failure. It's the tool declining to give you a number that would be wrong.
Earlier versions didn't refuse, and quietly produced nonsense — that's why these checks exist.

## If something goes wrong

- *ERROR reading the file* — the file is damaged or isn't really the format its name claims.
- Anything else unexpected: the exit code is `0` if it ran, `2` if it refused on purpose, `3` if
  it couldn't read the file.

Nothing you put in this folder is ever uploaded, shared, or committed to the repository.
