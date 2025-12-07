#!/usr/bin/env python3
"""
共识序列与参考序列比对及坐标转换整合脚本   


功能:
- 读取所有样本的共识序列
- 与对应的参考序列进行比对
- 检测差异（突变）
- 将基于旋转参考序列的坐标转换回原始坐标系
- 生成包含原始坐标的差异结果表格

使用方法:
    python consensus_comparison_with_transform.py [--test SAMPLE_ID] [--jobs N] [--output-format {csv|tsv|json}]
"""

import os
import sys
import argparse
import csv
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import logging

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("错误: 需要安装biopython库")
    print("请运行: pip install biopython")
    sys.exit(1)

# 配置路径
CONSENSUS_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/Medaka_consensus_8/r2"
REFERENCE_FILE = "/Users/jck/Desktop/ref/gold/reordered_1824_3221_reference_A2763.fasta"
OUTPUT_DIR = "/Users/jck/Desktop/workflow-qc/raw-sup-accuracy-fastq/git/hbv-nanopore-pipeline/consensus_ref_viewa_9"

# HBV基因组参数
GENOME_LENGTH = 3221
DEFAULT_OFFSET = 1823  # DR1位置对应原始序列1824位置

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(OUTPUT_DIR, 'integrated_analysis.log'))
    ]
)
logger = logging.getLogger(__name__)

def load_reference_sequence(ref_file):
    """加载固定参考序列文件，返回 (序列ID, 序列)"""
    logger.info(f"正在加载固定参考序列文件: {ref_file}")
    
    try:
        for record in SeqIO.parse(ref_file, "fasta"):
            ref_id = record.id
            ref_seq = str(record.seq).upper()
            logger.info(f"成功加载参考序列: {ref_id}, 长度: {len(ref_seq)}")
            return ref_id, ref_seq
        
        logger.error("参考序列文件为空")
        return None, None
    
    except Exception as e:
        logger.error(f"加载参考序列失败: {e}")
        return None, None

def get_sample_list():
    """获取所有样本ID列表"""
    samples = []
    consensus_path = Path(CONSENSUS_DIR)
    
    if not consensus_path.exists():
        logger.error(f"共识序列目录不存在: {CONSENSUS_DIR}")
        return samples
    
    for sample_dir in consensus_path.iterdir():
        if sample_dir.is_dir() and sample_dir.name.isdigit():
            consensus_file = sample_dir / "consensus.fasta"
            if consensus_file.exists():
                samples.append(sample_dir.name)
    
    samples.sort()
    logger.info(f"发现 {len(samples)} 个样本")
    return samples

def load_consensus_sequence(sample_id):
    """加载单个样本的共识序列"""
    consensus_file = Path(CONSENSUS_DIR) / sample_id / "consensus.fasta"
    
    if not consensus_file.exists():
        logger.error(f"样本 {sample_id} 的共识序列文件不存在: {consensus_file}")
        return None, None
    
    try:
        for record in SeqIO.parse(consensus_file, "fasta"):
            # 返回第一个序列记录
            return record.id, str(record.seq).upper()
    except Exception as e:
        logger.error(f"读取样本 {sample_id} 共识序列失败: {e}")
        return None, None

def transform_coordinate(current_pos, offset=DEFAULT_OFFSET, genome_length=GENOME_LENGTH):
    """
    转换坐标：从旋转后坐标转换回原始坐标
    
    Args:
        current_pos: 当前坐标(基于旋转后序列，1-based)
        offset: 偏移量 
        genome_length: 基因组长度
        
    Returns:
        original_pos: 原始坐标(1-based)
    """
    original_pos = ((current_pos - 1) + offset) % genome_length + 1
    return original_pos

def compare_sequences_with_transform(seq1, seq2, offset=DEFAULT_OFFSET):
    """
    比较两个序列并直接转换坐标，返回差异列表
    
    返回格式: [
        {
            'rotated_position': int,   # 旋转后的位置 (1-based)
            'original_position': int,  # 原始位置 (1-based)
            'type': str,              # 'SNP', 'insertion', 'deletion'  
            'reference': str,         # 参考序列碱基
            'consensus': str,         # 共识序列碱基
            'length': int            # 变异长度
        }
    ]
    """
    differences = []
    
    # 确保序列长度一致，如果不一致需要特殊处理
    ref_len = len(seq1)
    con_len = len(seq2)
    max_len = max(ref_len, con_len)
    
    # 填充较短的序列
    if ref_len < max_len:
        seq1 += 'N' * (max_len - ref_len)
    if con_len < max_len:
        seq2 += 'N' * (max_len - con_len)
    
    i = 0
    while i < max_len:
        if seq1[i] != seq2[i]:
            rotated_pos = i + 1
            original_pos = transform_coordinate(rotated_pos, offset)
            
            # 检查是否为简单SNP
            if seq1[i] != '-' and seq2[i] != '-':
                differences.append({
                    'rotated_position': rotated_pos,
                    'original_position': original_pos,
                    'type': 'SNP',
                    'reference': seq1[i],
                    'consensus': seq2[i],
                    'length': 1
                })
            else:
                # 处理插入缺失 (简化版本)
                if seq1[i] == '-':
                    differences.append({
                        'rotated_position': rotated_pos,
                        'original_position': original_pos,
                        'type': 'insertion',
                        'reference': '-',
                        'consensus': seq2[i],
                        'length': 1
                    })
                elif seq2[i] == '-':
                    differences.append({
                        'rotated_position': rotated_pos,
                        'original_position': original_pos,
                        'type': 'deletion',
                        'reference': seq1[i],
                        'consensus': '-',
                        'length': 1
                    })
        i += 1
    
    return differences

def process_single_sample(sample_id, ref_id, reference_seq, output_format='csv'):
    """处理单个样本的比对和坐标转换任务（使用固定参考序列）"""
    logger.info(f"开始处理样本: {sample_id}")
    
    # 读取共识序列
    consensus_header, consensus_seq = load_consensus_sequence(sample_id)
    
    if not consensus_header or not consensus_seq:
        logger.error(f"样本 {sample_id}: 无法读取共识序列")
        return None
    
    # 使用固定参考序列（不需要匹配）
    logger.info(f"样本 {sample_id}: 使用固定参考序列 {ref_id}, 长度 = {len(reference_seq)}, 共识序列长度 = {len(consensus_seq)}")
    
    # 进行序列比对并直接转换坐标
    differences = compare_sequences_with_transform(reference_seq, consensus_seq, DEFAULT_OFFSET)
    
    # 准备结果数据
    result = {
        'sample_id': sample_id,
        'reference_id': ref_id,
        'reference_length': len(reference_seq),
        'consensus_length': len(consensus_seq),
        'offset_used': DEFAULT_OFFSET,
        'total_differences': len(differences),
        'differences': differences,
        'processing_time': datetime.now().isoformat()
    }
    
    # 保存结果文件
    output_file = Path(OUTPUT_DIR) / f"{sample_id}_comparison_transformed.{output_format}"
    save_result(result, output_file, output_format)
    
    logger.info(f"样本 {sample_id}: 检测到 {len(differences)} 个差异（已转换坐标），结果保存至 {output_file}")
    return result

def save_result(result, output_file, output_format='csv'):
    """保存比对和坐标转换结果到文件"""
    try:
        if output_format == 'csv':
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                
                # 写入表头
                writer.writerow(['pos', 'mutat', 'ref', 'alt', 'length'])
                
                # 写入差异详情
                for diff in result['differences']:
                    writer.writerow([
                        diff['original_position'],
                        diff['type'], 
                        diff['reference'],
                        diff['consensus'],
                        diff['length']
                    ])
        
        elif output_format == 'tsv':
            with open(output_file, 'w', encoding='utf-8') as f:
                # 写入头部信息
                f.write(f"# 样本比对结果（含坐标转换）\n")
                f.write(f"样本ID\t{result['sample_id']}\n")
                f.write(f"参考序列ID\t{result['reference_id']}\n")
                f.write(f"参考序列长度\t{result['reference_length']}\n")
                f.write(f"共识序列长度\t{result['consensus_length']}\n")
                f.write(f"坐标偏移量\t{result['offset_used']}\n")
                f.write(f"差异总数\t{result['total_differences']}\n")
                f.write(f"处理时间\t{result['processing_time']}\n")
                f.write(f"\n")
                
                # 写入差异详情
                f.write("原始位置\t旋转后位置\t类型\t参考碱基\t共识碱基\t长度\n")
                for diff in result['differences']:
                    f.write(f"{diff['original_position']}\t{diff['rotated_position']}\t{diff['type']}\t{diff['reference']}\t{diff['consensus']}\t{diff['length']}\n")
        
        elif output_format == 'json':
            with open(output_file.with_suffix('.json'), 'w', encoding='utf-8') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
        
    except Exception as e:
        logger.error(f"保存结果文件失败: {e}")

def generate_summary_report(results, output_format='csv'):
    """生成汇总报告"""
    logger.info("生成汇总报告...")
    
    summary_file = Path(OUTPUT_DIR) / f"integrated_summary.{output_format}"
    
    try:
        if output_format == 'csv':
            with open(summary_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                
                # 写入头部
                writer.writerow(['样本ID', '参考序列ID', '参考长度', '共识长度', '偏移量', '差异总数', '处理时间'])
                
                # 写入结果数据
                for result in results:
                    if result:  # 过滤掉None结果
                        writer.writerow([
                            result['sample_id'],
                            result['reference_id'], 
                            result['reference_length'],
                            result['consensus_length'],
                            result['offset_used'],
                            result['total_differences'],
                            result['processing_time']
                        ])
        
        elif output_format == 'tsv':
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write("样本ID\t参考序列ID\t参考长度\t共识长度\t偏移量\t差异总数\t处理时间\n")
                
                for result in results:
                    if result:
                        f.write(f"{result['sample_id']}\t{result['reference_id']}\t{result['reference_length']}\t{result['consensus_length']}\t{result['offset_used']}\t{result['total_differences']}\t{result['processing_time']}\n")
        
        logger.info(f"汇总报告保存至: {summary_file}")
        
        # 生成详细的变异位置分布报告
        generate_variant_distribution_report(results)
        
    except Exception as e:
        logger.error(f"生成汇总报告失败: {e}")

def generate_variant_distribution_report(results):
    """生成变异位置分布报告"""
    try:
        position_counts = {}
        variant_types = {'SNP': 0, 'insertion': 0, 'deletion': 0}
        
        for result in results:
            if result:
                for diff in result['differences']:
                    pos = diff['original_position']
                    if pos not in position_counts:
                        position_counts[pos] = {'count': 0, 'samples': []}
                    position_counts[pos]['count'] += 1
                    position_counts[pos]['samples'].append(result['sample_id'])
                    variant_types[diff['type']] += 1
        
        # 保存位置分布报告
        dist_file = Path(OUTPUT_DIR) / "variant_position_distribution.txt"
        with open(dist_file, 'w', encoding='utf-8') as f:
            f.write("# 变异位置分布报告\n")
            f.write(f"# 生成时间: {datetime.now().isoformat()}\n\n")
            
            f.write("## 变异类型统计\n")
            for vtype, count in variant_types.items():
                f.write(f"{vtype}: {count}\n")
            f.write(f"总计: {sum(variant_types.values())}\n\n")
            
            f.write("## 高频变异位置 (≥5个样本)\n")
            f.write("位置\t样本数\t样本列表\n")
            
            # 按频率排序
            sorted_positions = sorted(position_counts.items(), 
                                    key=lambda x: x[1]['count'], 
                                    reverse=True)
            
            for pos, data in sorted_positions:
                if data['count'] >= 5:
                    sample_list = ','.join(data['samples'][:10])
                    if len(data['samples']) > 10:
                        sample_list += f"...等{len(data['samples'])}个样本"
                    f.write(f"{pos}\t{data['count']}\t{sample_list}\n")
            
        logger.info(f"变异位置分布报告保存至: {dist_file}")
        
    except Exception as e:
        logger.error(f"生成变异位置分布报告失败: {e}")

def main():
    global DEFAULT_OFFSET  # 声明全局变量（必须在引用前）
    
    parser = argparse.ArgumentParser(description='共识序列与参考序列比对及坐标转换整合脚本')
    parser.add_argument('--test', type=str, help='测试单个样本 (样本ID)')
    parser.add_argument('--jobs', type=int, default=4, help='并行处理任务数 (默认: 4)')
    parser.add_argument('--output-format', choices=['csv', 'tsv', 'json'], default='csv', help='输出文件格式')
    parser.add_argument('--list', action='store_true', help='列出所有可用样本')
    parser.add_argument('--offset', type=int, default=DEFAULT_OFFSET, help=f'坐标转换偏移量 (默认: {DEFAULT_OFFSET})')
    
    args = parser.parse_args()
    
    # 如果指定了自定义偏移量，更新全局变量
    if args.offset != DEFAULT_OFFSET:
        DEFAULT_OFFSET = args.offset
        logger.info(f"使用自定义偏移量: {DEFAULT_OFFSET}")
    
    logger.info("="*60)
    logger.info("共识序列与参考序列比对及坐标转换整合分析 v1.1 (固定参考序列版)")
    logger.info(f"共识序列目录: {CONSENSUS_DIR}")
    logger.info(f"参考序列文件: {REFERENCE_FILE}")
    logger.info(f"坐标偏移量: {DEFAULT_OFFSET}")
    logger.info("="*60)
    
    # 获取样本列表
    samples = get_sample_list()
    
    if args.list:
        print(f"发现 {len(samples)} 个可用样本:")
        for i, sample in enumerate(samples, 1):
            print(f"{i:3d}. {sample}")
        return
    
    if not samples:
        logger.error("未发现任何可用样本")
        return
    
    # 加载固定参考序列
    ref_id, reference_seq = load_reference_sequence(REFERENCE_FILE)
    if not ref_id or not reference_seq:
        logger.error("无法加载参考序列，程序退出")
        return
    
    # 处理样本
    if args.test:
        # 单样本测试
        if args.test not in samples:
            logger.error(f"样本 {args.test} 不存在")
            return
        
        result = process_single_sample(args.test, ref_id, reference_seq, args.output_format)
        if result:
            print(f"\n样本 {args.test} 处理完成:")
            print(f"  - 参考序列: {result['reference_id']}")
            print(f"  - 参考长度: {result['reference_length']}")
            print(f"  - 共识长度: {result['consensus_length']}")
            print(f"  - 坐标偏移量: {result['offset_used']}")
            print(f"  - 差异总数: {result['total_differences']}")
            
            # 显示前几个变异示例
            if result['differences']:
                print(f"  - 变异示例:")
                for i, diff in enumerate(result['differences'][:5]):
                    print(f"    原始位置 {diff['original_position']} (旋转后 {diff['rotated_position']}): "
                          f"{diff['type']} {diff['reference']}→{diff['consensus']}")
                if len(result['differences']) > 5:
                    print(f"    ... 共 {len(result['differences'])} 个变异")
    else:
        # 批量处理
        print(f"\n准备处理 {len(samples)} 个样本 (并行数: {args.jobs})")
        print(f"使用固定参考序列: {ref_id}")
        print(f"使用坐标偏移量: {DEFAULT_OFFSET}")
        
        # 确认继续
        confirm = input("是否继续? [y/N]: ")
        if confirm.lower() != 'y':
            logger.info("用户取消操作")
            return
        
        results = []
        
        # 并行处理
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {
                executor.submit(process_single_sample, sample, ref_id, reference_seq, args.output_format): sample
                for sample in samples
            }
            
            completed = 0
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1
                    logger.info(f"进度: {completed}/{len(samples)} ({100*completed/len(samples):.1f}%)")
                except Exception as e:
                    logger.error(f"处理样本 {sample} 时发生错误: {e}")
        
        # 生成汇总报告
        generate_summary_report(results, args.output_format)
        
        # 统计信息
        successful_results = [r for r in results if r is not None]
        total_differences = sum(r['total_differences'] for r in successful_results)
        
        print(f"\n处理完成!")
        print(f"  - 成功处理: {len(successful_results)}/{len(samples)} 样本")
        print(f"  - 总差异数: {total_differences}")
        print(f"  - 平均差异: {total_differences/len(successful_results):.1f}" if successful_results else "  - 平均差异: 0")
        print(f"  - 结果目录: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
