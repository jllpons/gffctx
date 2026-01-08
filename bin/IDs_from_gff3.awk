#!/opt/homebrew/bin/gawk -f

BEGIN {
    FS = "\t";
}

{
    # skip header lines
    if ($0 ~ /^#/) {
        next;
    }
    # In GFF3, attributes are in the 9th column
    attribtues = split($9, a, ";");
    gene_id = "";
    for (i = 1; i <= attribtues; i++) {
        split(a[i], b, "=");
        if (b[1] == "ID") {
            gene_id = b[2];
        }
    }
    if (gene_id != "") {
        print gene_id
    }
}
