#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
画总体 OTR 图：

1. iFP_ratio
   FP = reads_with_any_iFP / total_reads
   X7 = reads_with_any_iFP_and_anchor / total_reads

2. iFP_spacer_ratio
   iFP_spacer = reads_with_any_iFP_spacer / total_reads
   X7 = reads_with_any_iFP_spacer_and_anchor / total_reads

输入：
    result/04_on_target_rates/overall_onTargetRate_ALL_samples.csv

输出：
    result/05_plots/
        iFP_ratio.png
        iFP_ratio.pdf
        iFP_spacer_ratio.png
        iFP_spacer_ratio.pdf
        overall_OTR_plot_data.csv

默认：
    如果图已经存在且非空，则跳过。
    --force 可以强制重画。
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 0. 默认路径
# ============================================================

DEFAULT_PROJECT_DIR = Path(
    "/data2/lindf/project/HPV_project/"
    "20260831_国妇婴组织样本_靶向建库_数据分析/"
    "On_target_rates"
)

DEFAULT_INPUT = (
    DEFAULT_PROJECT_DIR
    / "result"
    / "04_on_target_rates"
    / "overall_onTargetRate_ALL_samples.csv"
)

DEFAULT_OUTDIR = (
    DEFAULT_PROJECT_DIR
    / "result"
    / "05_plots"
)


# ============================================================
# 1. 参数
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Plot overall iFP and iFP_spacer on-target rates."
    )

    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
        help="overall_onTargetRate_ALL_samples.csv"
    )

    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help="Output directory"
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force redraw existing figures"
    )

    return parser.parse_args()


# ============================================================
# 2. 自然排序
#    Tissue1, Tissue15, Tissue26...
# ============================================================

def natural_sort_key(text):

    import re

    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", str(text))
    ]


# ============================================================
# 3. 通用双柱图函数
# ============================================================

def plot_grouped_bar(
    df,
    col1,
    col2,
    label1,
    label2,
    title,
    output_prefix,
    force=False
):

    png_file = output_prefix.with_suffix(".png")
    pdf_file = output_prefix.with_suffix(".pdf")

    # --------------------------------------------------------
    # 存在即跳过
    # --------------------------------------------------------

    if (
        not force
        and png_file.exists()
        and png_file.stat().st_size > 0
        and pdf_file.exists()
        and pdf_file.stat().st_size > 0
    ):
        print(f"[SKIP] Figure already exists: {png_file}")
        return


    samples = df["sample_id"].astype(str).tolist()

    value1 = df[col1].astype(float).to_numpy()
    value2 = df[col2].astype(float).to_numpy()

    x = np.arange(len(samples))

    width = 0.36


    # --------------------------------------------------------
    # 图尺寸
    # --------------------------------------------------------

    fig_width = max(9, len(samples) * 0.9)

    fig, ax = plt.subplots(
        figsize=(fig_width, 6.2)
    )


    # --------------------------------------------------------
    # 柱状图
    #
    # 为了和你原来的图接近：
    # FP / iFP_spacer = 浅蓝
    # X7              = 浅紫
    # --------------------------------------------------------

    bars1 = ax.bar(
        x - width / 2,
        value1,
        width,
        label=label1,
        color="#B7E1EA"
    )

    bars2 = ax.bar(
        x + width / 2,
        value2,
        width,
        label=label2,
        color="#DED1EE"
    )


    # --------------------------------------------------------
    # 百分比标签
    # --------------------------------------------------------

    def add_labels(bars):

        for bar in bars:

            h = bar.get_height()

            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 1.0,
                f"{h:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10
            )

    add_labels(bars1)
    add_labels(bars2)


    # --------------------------------------------------------
    # 坐标轴
    # --------------------------------------------------------

    max_value = max(
        float(np.max(value1)),
        float(np.max(value2))
    )

    # 至少显示到100附近
    upper = max(
        100,
        np.ceil((max_value + 10) / 10) * 10
    )

    ax.set_ylim(
        -2,
        upper
    )

    ax.set_ylabel(
        "% total reads",
        fontsize=15
    )

    ax.set_xlabel(
        "Library",
        fontsize=15
    )

    ax.set_title(
        title,
        fontsize=18,
        loc="left",
        pad=15
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        samples,
        rotation=45,
        ha="right",
        fontsize=11
    )

    ax.tick_params(
        axis="y",
        labelsize=11
    )


    # --------------------------------------------------------
    # 图例
    # --------------------------------------------------------

    ax.legend(
        title="type",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left"
    )


    # --------------------------------------------------------
    # 去掉上和右边框
    # --------------------------------------------------------

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)


    # 不画背景网格
    ax.grid(False)


    plt.tight_layout()


    # --------------------------------------------------------
    # 保存
    # --------------------------------------------------------

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(f"[DONE] {png_file}")
    print(f"[DONE] {pdf_file}")


# ============================================================
# 4. 主程序
# ============================================================

def main():

    args = parse_args()

    input_file = Path(args.input)

    outdir = Path(args.outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # 检查输入
    # --------------------------------------------------------

    if not input_file.exists():

        raise FileNotFoundError(
            f"Input file not found:\n{input_file}"
        )

    if input_file.stat().st_size == 0:

        raise ValueError(
            f"Input file is empty:\n{input_file}"
        )


    print("=" * 70)
    print("Plot overall OTR")
    print("=" * 70)

    print("Input:")
    print(input_file)

    print()

    print("Output:")
    print(outdir)

    print()


    # --------------------------------------------------------
    # 读取数据
    # --------------------------------------------------------

    df = pd.read_csv(
        input_file
    )


    # --------------------------------------------------------
    # 检查必须列
    # --------------------------------------------------------

    required_cols = {
        "sample_id",
        "total_reads",

        "any_iFP_rate_pct",
        "any_iFP_anchor_rate_pct",

        "any_iFP_spacer_rate_pct",
        "any_iFP_spacer_anchor_rate_pct"
    }

    missing_cols = (
        required_cols
        - set(df.columns)
    )

    if missing_cols:

        raise ValueError(
            "Missing required columns:\n"
            + "\n".join(sorted(missing_cols))
        )


    # --------------------------------------------------------
    # 排序
    # --------------------------------------------------------

    df["sample_id"] = (
        df["sample_id"]
        .astype(str)
    )

    df = (
        df
        .sort_values(
            by="sample_id",
            key=lambda x: x.map(natural_sort_key)
        )
        .reset_index(drop=True)
    )


    # --------------------------------------------------------
    # 输出画图数据
    # --------------------------------------------------------

    plot_data = df[
        [
            "sample_id",
            "total_reads",

            "any_iFP_rate_pct",
            "any_iFP_anchor_rate_pct",

            "any_iFP_spacer_rate_pct",
            "any_iFP_spacer_anchor_rate_pct"
        ]
    ].copy()


    plot_data.columns = [
        "Library",
        "total_reads",

        "FP_pct",
        "FP_X7_pct",

        "iFP_spacer_pct",
        "iFP_spacer_X7_pct"
    ]


    plot_data_file = (
        outdir
        / "overall_OTR_plot_data.csv"
    )

    plot_data.to_csv(
        plot_data_file,
        index=False
    )

    print(
        f"[DONE] Plot data: "
        f"{plot_data_file}"
    )


    # ========================================================
    # Figure 1
    # iFP_ratio
    # ========================================================

    plot_grouped_bar(
        df=df,

        col1="any_iFP_rate_pct",

        col2="any_iFP_anchor_rate_pct",

        label1="FP",

        label2="X7",

        title="iFP_ratio",

        output_prefix=(
            outdir
            / "iFP_ratio"
        ),

        force=args.force
    )


    # ========================================================
    # Figure 2
    # iFP_spacer_ratio
    # ========================================================

    plot_grouped_bar(
        df=df,

        col1="any_iFP_spacer_rate_pct",

        col2="any_iFP_spacer_anchor_rate_pct",

        label1="iFP_spacer",

        label2="X7",

        title="iFP_spacer_ratio",

        output_prefix=(
            outdir
            / "iFP_spacer_ratio"
        ),

        force=args.force
    )


    print()

    print("=" * 70)
    print("All plots done.")
    print("=" * 70)


if __name__ == "__main__":
    main()