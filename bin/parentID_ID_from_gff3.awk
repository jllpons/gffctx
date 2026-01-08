#!/opt/homebrew/bin/gawk -f

BEGIN {
    FS = "\t"
    OFS = ","
}

# skip header lines
/^#/ { next }

{
    n = split($9, a, ";")
    id = ""
    parent_id = ""

    for (i = 1; i <= n; i++) {
        split(a[i], b, "=")
        if (b[1] == "ID")     id = b[2]
        else if (b[1] == "Parent") parent_id = b[2]
    }

    if (parent_id == "") {
        print "", id
    }
    else if (id != "" && parent_id != "") {
        np = split(parent_id, p, ",")
        for (j = 1; j <= np; j++) {
            if (p[j] != "")
                print p[j], id
        }
    }
}
