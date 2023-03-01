SRAID="SRR14800536"
prefetch $SRAID --max-size u
fasterq-dump --split-files --include-technical "$SRAID"sra
