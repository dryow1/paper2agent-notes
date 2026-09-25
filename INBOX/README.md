# INBOX

Put **one** file here, then run the door.

**Accepted:** a single `.h5ad`, or a 10x Genomics `.h5`.
**Not accepted:** a 10x mtx folder (`matrix.mtx` + `barcodes` + `features`) — the Scanpy tools
do not read that layout, and converting it is your decision, not the door's.

**Limits:** file ≤ 80 MB, and the estimated memory peak must stay under 8 GB on this 16 GB box.
`X` must be raw counts — logged, scaled or target-summed data is refused with the reason.

```bash
cd ~/Work/paper2agent && . ./env.sh && . ./Scanpy_Agent/project.env
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py --inspect-only   # look first
cap10g "$PROJECT_PYTHON" Scanpy_Agent/src/inbox.py --mt-prefix MT-  # then run QC
```

The door runs quality control and stops. It does not cluster, annotate or name anything, and it
never downloads a dataset — an empty INBOX is a refusal.
