#!/usr/bin/env python3
"""
DEG和Pathway分析结果总结表格生成器
"""

import pandas as pd
import os
from pathlib import Path

# 定义路径
BASE_DIR = Path("/disk2/cai113/data/stateTrans")
DEG_DIR = BASE_DIR / "06.DEG_analysis/DEG_results"
PATHWAY_DIR = BASE_DIR / "07.pathway_analysis/GSEA_results/batch_all"
OUT_DIR = BASE_DIR / "07.pathway_analysis"

# 定义比较组
COMPARISONS = [
    "c1_vs_c1_star",
    "c2_vs_c1_star", 
    "c3_vs_c1_star",
    "c2_vs_c1",
    "c3_vs_c1",
    "c3_vs_c2"
]

# 定义pathway类型
PATHWAY_TYPES = ["GO_BP", "MSigDB_C2", "MSigDB_C7", "MSigDB_C8", "MSigDB_H"]


def parse_deg_summary(comparison):
    """解析DEG summary文件"""
    summary_file = DEG_DIR / comparison / "DEG_summary.txt"
    
    if not summary_file.exists():
        return {"all": 0, "up": 0, "down": 0}
    
    with open(summary_file) as f:
        content = f.read()
    
    # 提取数字
    lines = content.split('\n')
    deg_all = up = down = 0
    
    for line in lines:
        if "Significant DEGs:" in line:
            deg_all = int(line.split(':')[1].strip())
        elif "Up-regulated:" in line:
            up = int(line.split(':')[1].strip())
        elif "Down-regulated:" in line:
            down = int(line.split(':')[1].strip())
    
    return {"all": deg_all, "up": up, "down": down}


def parse_pathway_file(filepath):
    """解析pathway文件，统计all/up/down数量"""
    if not filepath.exists():
        return {"all": 0, "up": 0, "down": 0}
    
    try:
        df = pd.read_csv(filepath, sep='\t')
        
        if df.empty or 'enrichmentScore' not in df.columns:
            return {"all": 0, "up": 0, "down": 0}
        
        total = len(df)
        up = (df['enrichmentScore'] > 0).sum()
        down = (df['enrichmentScore'] < 0).sum()
        
        return {"all": total, "up": up, "down": down}
    
    except Exception as e:
        print(f"警告: 读取文件 {filepath} 时出错: {e}")
        return {"all": 0, "up": 0, "down": 0}


def generate_summary_table():
    """生成总结表格"""
    
    results = []
    
    for comp in COMPARISONS:
        row = {"Comparison": comp}
        
        # 添加DEG统计
        deg_stats = parse_deg_summary(comp)
        row["DEG_all"] = deg_stats["all"]
        row["DEG_up"] = deg_stats["up"]
        row["DEG_down"] = deg_stats["down"]
        
        # 添加各类pathway统计
        for ptype in PATHWAY_TYPES:
            filename = f"GSEA_{ptype}.txt"
            filepath = PATHWAY_DIR / comp / filename
            
            pathway_stats = parse_pathway_file(filepath)
            row[f"{ptype}_all"] = pathway_stats["all"]
            row[f"{ptype}_up"] = pathway_stats["up"]
            row[f"{ptype}_down"] = pathway_stats["down"]
        
        results.append(row)
    
    # 创建DataFrame
    df = pd.DataFrame(results)
    
    # 调整列顺序
    cols = ["Comparison", "DEG_all", "DEG_up", "DEG_down"]
    for ptype in PATHWAY_TYPES:
        cols.extend([f"{ptype}_all", f"{ptype}_up", f"{ptype}_down"])
    
    df = df[cols]
    
    return df

def main():
    """主函数"""
    print("开始生成DEG和Pathway分析总结表格...")
    
    # 生成表格
    summary_df = generate_summary_table()
    
    # 保存为CSV和Excel
    output_csv = OUT_DIR / "DEG_Pathway_Summary.csv"
    output_xlsx = OUT_DIR / "DEG_Pathway_Summary.xlsx"
    
    summary_df.to_csv(output_csv, index=False)
    summary_df.to_excel(output_xlsx, index=False)
    
    print(f"\n总结表格已保存:")
    print(f"  CSV: {output_csv}")
    print(f"  Excel: {output_xlsx}")

if __name__ == "__main__":
    main()