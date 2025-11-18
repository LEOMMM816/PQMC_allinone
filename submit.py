#!/usr/bin/env python3
import sys
import shutil
import re
import os
from pathlib import Path
from datetime import date
import subprocess  
# ----------------- 工具函数 ----------------- #

def parse_value(raw):
    """
    将单元格字符串解析成描述字典：
      - scalar: 普通数字/字符串
      - tuple: "(8,8,8)"  -> inner = "8,8,8"
      - list:  "[1,2,3]"
      - range: "[1:0.1:1.7]"
    """
    raw = raw.strip()
    if raw.startswith("(") and raw.endswith(")"):
        inner = raw[1:-1].strip()
        return {"kind": "tuple", "inner": inner}

    if raw.startswith("[") and raw.endswith("]"):
        inner = raw[1:-1].strip()
        # 形如 start:step:stop 的 range
        if inner.count(":") == 2 and "," not in inner:
            start_s, step_s, stop_s = [s.strip() for s in inner.split(":")]
            return {
                "kind": "range",
                "start": float(start_s),
                "step": float(step_s),
                "stop": float(stop_s),
            }
        else:
            # 普通 list: 按逗号切
            parts = [p.strip() for p in inner.split(",") if p.strip()]
            return {"kind": "list", "values": parts}

    # 默认：标量
    return {"kind": "scalar", "value": raw}


def format_float(x):
    # 用通用格式，避免 1.7000000000002 这种
    return f"{x:g}"


def expand_for_blocks(desc, nblock, varname):
    """
    根据 nblock 展开到每个 block 的取值列表（长度 = nblock）。
    desc 为 parse_value 的结果。
    """
    kind = desc["kind"]

    if kind == "scalar":
        return [desc["value"]] * nblock

    if kind == "tuple":
        # 去掉括号后的内容整体写入
        return [desc["inner"]] * nblock

    if kind == "list":
        values = desc["values"]
        if len(values) != nblock:
            raise ValueError(
                f"变量 '{varname}' 的 list 元素数为 {len(values)}，"
                f"但 nblock = {nblock}，不一致，程序终止。"
            )
        return values

    if kind == "range":
        start = desc["start"]
        step = desc["step"]
        stop = desc["stop"]
        if step == 0:
            raise ValueError(f"变量 '{varname}' 的 range 步长为 0，非法。")

        # 先计算理论元素个数
        expected_n = int(round((stop - start) / step)) + 1
        values = []
        for i in range(expected_n):
            val = start + i * step
            # 防止数值误差轻微越界
            if step > 0 and val > stop + 1e-9:
                break
            if step < 0 and val < stop - 1e-9:
                break
            values.append(format_float(val))

        if len(values) != nblock:
            raise ValueError(
                f"变量 '{varname}' 的 range [{start}:{step}:{stop}] "
                f"生成了 {len(values)} 个元素，但 nblock = {nblock}，不一致，程序终止。"
            )
        return values

    raise ValueError(f"未知 kind: {kind}")


def replace_in_nml(text, varname, new_value):
    """
    在 namelist 文本中，把一行形如
      varname = 原值
    改成
      varname = new_value
    如果找不到该变量，则在末尾追加一行。
    """
    pattern = re.compile(rf"^(\s*{re.escape(varname)}\s*=\s*).*$", re.MULTILINE)
    if pattern.search(text):
        # 用函数替换，避免 \1{new_value} 变成 \18 之类的问题
        def _repl(m):
            return m.group(1) + str(new_value)
        new_text = pattern.sub(_repl, text, count=1)
        return new_text
    else:
        # 没找到变量，就追加到文件末尾
        print(f"警告: 在 nml 文件中找不到变量 '{varname}'，将在末尾追加该变量。", file=sys.stderr)
        append_line = f"\n{varname} = {new_value}\n"
        return text + append_line



# ----------------- 主逻辑 ----------------- #

def main():

    # 这里写你真正要 sbatch 的脚本，比如 run_sbatch.sh
    JOB_SCRIPT = "script/mpi_32c.sh"
    JOB_SCRIPT_local = "script/testlocal.sh"
    is_local = False
    if len(sys.argv) < 2:
        param_path = Path("jobs_param.tsv")
    else:
        param_path = Path(sys.argv[1])

   
    if not param_path.exists():
        print(f"找不到参数文件: {param_path}", file=sys.stderr)
        sys.exit(1)
    

    # 读取 param.tsv，忽略空行和以 # 开头的行
    with param_path.open("r") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]

    if not lines:
        print("param.tsv 为空。", file=sys.stderr)
        sys.exit(1)

    header_line = lines[0].strip()
    headers = header_line.split()
    if "model" not in headers:
        print("表头中找不到 'model' 列。", file=sys.stderr)
        sys.exit(1)

    idx_model = headers.index("model")
    sbatch_cols = headers[:idx_model]      # divide 之前
    nml_cols = headers[idx_model + 1:]     # divide 之后

    if "job-name" not in headers:
        print("表头中没有 'job-name' 列。", file=sys.stderr)
        sys.exit(1)
    if "nblock" not in headers:
        print("表头中没有 'nblock' 列。", file=sys.stderr)
        sys.exit(1)
    # 获取今天日期字符串
    today_str = date.today().strftime("%Y%m%d")
    # 模板 nml 文件路径
    template_dir = str('model/')
    suffix = str('.nml')
    print(f"# 读取参数文件: {param_path}")
    print()

    # 逐行（从第二行开始）处理
    for line_idx, line in enumerate(lines[1:], start=2):
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) != len(headers):
            print(
                f"第 {line_idx} 行的列数 ({len(parts)}) 与表头 ({len(headers)}) 不一致，程序终止。",
                file=sys.stderr,
            )
            sys.exit(1)

        row = dict(zip(headers, parts))

        jobname = row["job-name"]
        try:
            nblock = int(row["nblock"])
        except ValueError:
            print(f"第 {line_idx} 行 nblock = '{row['nblock']}' 不是整数。", file=sys.stderr)
            sys.exit(1)

        print(f"== 处理第 {line_idx - 1} 个任务: jobname={jobname}")
        # 模板文件路径
        base_name = row["model"]
        template_path = Path(f"{template_dir}/{base_name}{suffix}")
        if not template_path.exists():
            print(f"找不到模板文件: {template_path}", file=sys.stderr)
            sys.exit(1)
        # 1) 组装 sbatch 命令（只打印，不执行）
        sbatch_opts = [f"--{col}={row[col]}" for col in sbatch_cols]
        
        # 2) 解析 divide 后各变量，并展开到每个 block
        per_var_values = {}  # var -> [v_block1, v_block2, ...]
        for var in nml_cols:
            raw_val = row[var]
            desc = parse_value(raw_val)
            values = expand_for_blocks(desc, nblock, var)
            per_var_values[var] = values

        # 3) 为该任务复制 nblock 份 nml 并做替换
        block_files = []
        name_str = f"{base_name}_{jobname}_{today_str}"
        for ib in range(nblock):
            new_name = f"{name_str}_b{ib+1}"
            new_path = Path(f"{template_dir}temp/{new_name}{suffix}")
            shutil.copy2(template_path, new_path)
            block_files.append(new_path)

        # 4) 在每个 nml 文件中替换变量
        for ib, nml_path in enumerate(block_files):
            text = nml_path.read_text()
            for var, values in per_var_values.items():
                new_val = values[ib]
                text = replace_in_nml(text, var, new_val)
            nml_path.write_text(text)
            print(f"  已生成并修改: {nml_path}")

        # 5) 生成并执行 sbatch 命令（确保此时 nml 文件已经就绪）
        if is_local:
            env = os.environ.copy()
            env["NBLOCK"] = str(nblock)  # 环境变量一定要是字符串
            env["NMLFILE"] = f"{template_dir}temp/{name_str}"
            cmd = [f"./{JOB_SCRIPT_local}"]
            print("本地测试命令：", " ".join(cmd))
            subprocess.run(cmd, env=env, check=True)
            print("# 本地测试任务处理完成。")
            exit(0)
        else:
            export_opt = f"--export=ALL,NBLOCK={row['nblock']},NMLFILE={template_dir}temp/{name_str}"
            cmd = ["sbatch"]  + sbatch_opts + [export_opt, JOB_SCRIPT]
            #cmd = ["sbatch"] + sbatch_opts + [JOB_SCRIPT]
            print("提交 sbatch 命令：", " ".join(cmd))
            # 真正执行 sbatch，如果出错会抛异常
            subprocess.run(cmd, check=True)
        print("# 全部任务处理完成。")

        # 6) 清空 temp 目录
        #for nml_path in block_files:
        #    nml_path.unlink()
        #print("已清理临时文件。")
        print() 
if __name__ == "__main__":
    main()
