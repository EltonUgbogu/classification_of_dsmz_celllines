
echo "Protected file:"
ls -lh /work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_count.rds

echo
echo "Deleted paths status:"
for p in \
/work/ugbogu/pipeline/data/nbl/GSE189367/sra \
/work/ugbogu/pipeline/data/nbl/GSE189367/fastq \
/work/ugbogu/pipeline/data/nbl/GSE189367/results/star \
/work/ugbogu/pipeline/data/nbl/SRP409177/sra \
/work/ugbogu/pipeline/data/nbl/SRP409177/fastq \
/work/ugbogu/pipeline/data/nbl/SRP409177/results/star \
/work/ugbogu/pipeline/data/nbl/GSE100148/sra \
/work/ugbogu/pipeline/data/nbl/GSE100148/fastq \
/work/ugbogu/pipeline/data/nbl/GSE100148/results/star
do
if [ -e "$p" ]; then
    echo "STILL EXISTS: $p"
else
    echo "DELETED: $p"
fi
done

echo
echo "Remaining top-level sizes:"
du -sh /work/ugbogu/pipeline/data/nbl/* 2>/dev/null | sort -h



failed=0
for p in \
/work/ugbogu/pipeline/data/nbl/GSE189367/sra \
/work/ugbogu/pipeline/data/nbl/GSE189367/fastq \
/work/ugbogu/pipeline/data/nbl/GSE189367/results/star \
/work/ugbogu/pipeline/data/nbl/SRP409177/sra \
/work/ugbogu/pipeline/data/nbl/SRP409177/fastq \
/work/ugbogu/pipeline/data/nbl/SRP409177/results/star \
/work/ugbogu/pipeline/data/nbl/GSE100148/sra \
/work/ugbogu/pipeline/data/nbl/GSE100148/fastq \
/work/ugbogu/pipeline/data/nbl/GSE100148/results/star
do
if [ -e "$p" ]; then
    echo "FAIL: $p still exists"
    failed=1
else
    echo "OK: $p deleted"
fi
done

test -f /work/ugbogu/pipeline/data/nbl/count_data/nbl_tumour_count.rds || {
echo "FAIL: protected file missing"
failed=1
}

exit $failed
