{
    unknown=homref=het=homalt=0
    for (i=10; i<=NF; i++) {
        split($i, a, ":")
        split(a[1], GT, "[/|]")
        if (GT[1]=="." && GT[2]==".") unknown++
        else if (GT[1]==0 && GT[2]==0) homref++
        else if (GT[1]==GT[2]) homalt++
        else het++
    }
    print $1, $2, $3":"$4"-"$5, $4, $5, unknown, homref, het, homalt, het+homalt
}