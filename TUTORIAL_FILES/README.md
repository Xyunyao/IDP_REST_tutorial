# System setup

```bash
awk '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top
```