#!/bin/bash

#set -euo pipefail
for ROW = 1,8
 PARAMS_FILE=jobs_param.tsv; OUTDIR="runs/job_${ROW}"; CFG="${OUTDIR}/input_${ROW}.nml";
mkdir -p "$OUTDIR"
gawk -v k="$ROW" -v cfg="$CFG" '
function trim(s){ sub(/^[[:space:]]+/, "", s); sub(/[[:space:]]+$/, "", s); return s }
function is_number(s){ s=trim(s); return s ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ }
function is_numlist(s, n,i){ s=trim(s); if (index(s,",")==0) return 0;
  n=split(s, A, /[[:space:]]*,[[:space:]]*/); for(i=1;i<=n;i++) if(!is_number(A[i])) return 0; return 1 }
function is_complex(s, t,n){ t=trim(s); if (t !~ /^\(.*\)$/) return 0;
  t=substr(t,2,length(t)-2); n=split(t, A, /[[:space:]]*,[[:space:]]*/); if(n!=2) return 0;
  return is_number(A[1]) && is_number(A[2]) }
function toval(v, t){ t=trim(v);
  if (is_complex(t)) return t;           # 复数：保留括号格式 (a,b)
  if (is_numlist(t)) return t;           # 数组：数字逗号列表 8,8,10
  if (is_number(t))  return t;           # 标量数字
  gsub(/"/, "\\\"", t); return "\"" t "\""}  # 其他当作字符串
BEGIN{ OFS="\t" }
/^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
{
  line = trim($0)
  if (++ln == 1) { N = split(line, H, /[[:space:]]+/); next }
  if (++d == k) {
    split(line, V, /[[:space:]]+/)
    print "&params" > cfg
    for (i=1; i<=N; i++) {
      val = (i in V) ? V[i] : ""
      printf "  %s = %s\n", H[i], toval(val) >> cfg
    }
    print "/" >> cfg; close(cfg); exit
  }
}
' "$PARAMS_FILE"
