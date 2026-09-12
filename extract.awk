BEGIN {
    FS = ","
    OFS = ","
    getline cols < "data/wanted_columns.txt"
    n = split(cols, c, ",")
}
{
    out = $c[1]
    for (i = 2; i <= n; i++) out = out OFS $c[i]
    print out
}
