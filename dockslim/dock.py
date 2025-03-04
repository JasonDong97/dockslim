from pathlib import Path
import subprocess

from rdkit import Chem
import rdkit2pdbqt
from openbabel import pybel, openbabel as ob
import numpy as np
from vina import Vina

import sys
import logging as log


def setup_logging(log_file: str, log_level: int = log.INFO):
    """
    配置全局日志设置

    参数:
        log_file (str): 日志文件路径（默认：当前目录的 app.log）
        log_level (int): 日志级别（默认：INFO）
    """
    # 创建日志目录（如果不存在）
    log_path = Path(log_file).parent
    log_path.mkdir(parents=True, exist_ok=True)

    # 定义日志格式
    formatter = log.Formatter(
        fmt='%(asctime)s - %(levelname)-8s - %(name)-15s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 获取 root logger
    logger = log.getLogger()
    logger.setLevel(log_level)

    # 避免重复添加 handler
    if logger.hasHandlers():
        logger.handlers.clear()

    # 创建文件 handler 并设置格式
    file_handler = log.FileHandler(log_file, encoding='utf-8')
    file_handler.setFormatter(formatter)

    # 创建控制台 handler 并设置格式
    console_handler = log.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)

    # 添加 handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)


def get_hetatm_docking_box(pdb_file):
    """
    计算PDB文件中异质原子(HETATM)的中心和大小，用于分子对接研究。

    参数:
    pdb_file (str): PDB文件路径, 该文件应包含分子结构信息。

    返回:
    tuple: 包含两个字典，第一个字典描述对接框的中心坐标，第二个字典描述对接框的尺寸。
    """
    # 使用RDKit从PDB文件创建分子对象，不进行结构规范化，保留氢原子
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("无法读取PDB文件")

    coords = []
    # 遍历分子中的所有原子，寻找HETATM记录
    for atom in mol.GetAtoms():
        # 更可靠的HETATM筛选方法
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is not None:
            # 方法1：检查是否为异质原子（推荐）
            if pdb_info.GetIsHeteroAtom():
                pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                coords.append([pos.x, pos.y, pos.z])
            # 方法2：直接检查记录类型（兼容旧版本RDKit）
            # if pdb_info.GetRecordType() == "HETATM":
            #     pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            #     coords.append([pos.x, pos.y, pos.z])

    if not coords:
        # 添加调试信息：列出所有残基名称
        residues = set()
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info:
                residues.add(
                    f"{info.GetResidueName()}({info.GetChainId()}:{info.GetResidueNumber()})")
        raise ValueError(
            f"未找到HETATM记录。文件中存在的残基: {', '.join(residues)}\n"
            "请确认：\n"
            "1. PDB文件是否包含配体/金属离子等异质原子\n"
            "2. 若使用共价配体，尝试将残基名称标记为'HETATM'"
        )

    # 将坐标转换为NumPy数组并计算中心和大小
    coords_array = np.round(np.array(coords), 3)
    min_vals = coords_array.min(axis=0)
    max_vals = coords_array.max(axis=0)

    # 计算并四舍五入对接框的中心坐标
    center = np.round((min_vals + max_vals) / 2, 3)
    # 计算并四舍五入对接框的尺寸
    size = np.round(max_vals - min_vals, 3)

    # 返回对接框的中心坐标和尺寸
    return (
        {'center_x': center[0], 'center_y': center[1], 'center_z': center[2]},
        {'size_x': size[0], 'size_y': size[1], 'size_z': size[2]}
    )


def split_protein_ligand(pdb_file):
    # 读取PDB文件
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("无法读取PDB文件")

    # 分割为多个残基
    residues = Chem.SplitMolByPDBResidues(mol)

    # 收集蛋白和配体的残基
    protein_mols = []
    ligand_mols = []

    for res_name, res_mol in residues.items():
        if res_mol.GetNumAtoms() == 0:
            continue

        # 获取第一个原子的信息
        atom = res_mol.GetAtomWithIdx(0)
        info = atom.GetMonomerInfo()

        if info is None:
            protein_mols.append(res_mol)
            continue

        # 检查是否为HETATM且不是水
        if info.GetIsHeteroAtom():
            res_name_pdb = info.GetResidueName().strip()
            if res_name_pdb != 'HOH':
                ligand_mols.append(res_mol)
        else:
            protein_mols.append(res_mol)

    # 合并分子
    def combine_mols(mol_list):
        if not mol_list:
            return None
        combined = mol_list[0]
        for m in mol_list[1:]:
            combined = Chem.CombineMols(combined, m)
        return combined

    protein_combined = combine_mols(protein_mols)
    ligand_combined = combine_mols(ligand_mols)

    workdir = Path(pdb_file).parent
    protein_output = workdir.joinpath("protein.pdb")
    ligand_output = workdir.joinpath("ligand.pdb")

    # 写入文件
    if protein_combined is not None:
        writer = Chem.PDBWriter(protein_output)
        writer.write(protein_combined)
        writer.close()

    if ligand_combined is not None:
        writer = Chem.PDBWriter(ligand_output)
        writer.write(ligand_combined)
        writer.close()

    return protein_output, ligand_output


def run_lepro(complex_pdb: str, out_pdb_file: str):
    '''
    将复合物PDB文件转换为lepro提供的输入文件
    :param complex_pdb: 复合物PDB文件
    :param out_pdb_file: 输出文件夹
    :return: lepro提供的输入文件
    '''
    file_path = Path(complex_pdb).absolute()
    workdir = file_path.parent

    subprocess.run(['lepro_linux_x86', file_path], cwd=workdir)
    workdir.joinpath("dock.in").unlink(missing_ok=True)

    Path(out_pdb_file).parent.mkdir(parents=True, exist_ok=True)
    workdir.joinpath("pro.pdb").rename(out_pdb_file)


def pdb_to_pdbqt(pdb_file: str, out_pdbqt_file: str):
    '''
    将PDB文件转换为PDBQT文件
    :param pdb_file: PDB文件路径
    :param out_pdbqt_file: 输出PDBQT文件路径
    :return:
    '''
    # 从PDB文件中读取拓扑结构
    receptor_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    # 将拓扑结构转换为PDBQT文件
    with open(out_pdbqt_file, 'w') as f:
        # 将拓扑结构转换为PDBQT文件
        lines = rdkit2pdbqt.MolToPDBQTBlock(receptor_mol, False, False, True)
        # 将PDBQT文件写入文件
        f.write(lines)


def pdbqt_to_sdf(pdbqt_file: str = None, output: str = None):
    '''
    将PDBQT文件转换为SDF文件
    :param pdbqt_file: PDBQT文件路径
    :param output: SDF文件路径
    :return: SDF文件路径
    '''
    # 从pdbqt文件中读取结果
    results = [m for m in pybel.readfile(filename=pdbqt_file, format='pdbqt')]
    # 创建一个输出文件
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)
    # 遍历结果，更新pose的数据
    for pose in results:

        pose.data.update({'Pose': pose.data['MODEL']})
        pose.data.update({'Score': pose.data['REMARK'].split()[2]})
        # 删除pose中的MODEL、REMARK、TORSDO
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']
        # 写入文件
        out.write(pose)
    # 关闭文件
    out.close()


def convert(input_file: str, input_format, out_file: str, out_format: str, add_h=True):

    # 由配体的pdb文件,使用 pybel库 转化格式得到mol2,sdf文件，到 './protein_files' 目录下
    mol = [m for m in pybel.readfile(
        filename=input_file, format=input_format)][0]
    if add_h:
        mol.addh()  # 加氢
    mol.write(format=out_format, filename=out_file, overwrite=True)


def vina_dock(complex_file, ligand_file, outdir):

    ligand_file_path = Path(ligand_file)
    workdir = Path(outdir).joinpath(ligand_file_path.stem)
    workdir.mkdir(parents=True, exist_ok=True)

    protein_file = workdir.joinpath("protein.pdb")
    protein_pdbqt = workdir.joinpath("protein.pdbqt")
    ligand_pdbqt = workdir.joinpath("ligand.pdbqt")
    dockout_pdbqt = workdir.joinpath("dockout.pdbqt")
    dockout_sdf = workdir.joinpath("dockout.sdf")

    log.info('Vina dock starts.\n')
    log.info(f"Workdir: {workdir}")
    # 示例用法
    center, size = get_hetatm_docking_box(complex_file)
    log.info(f"box center coordinate: {center}" )
    log.info(f"box size coordinate: {size}")

    run_lepro(complex_file, protein_file)

    pdb_to_pdbqt(protein_file, protein_pdbqt)
    convert(str(ligand_file), "mol2", str(ligand_pdbqt), "pdbqt")

    v = Vina()
    v.set_receptor(str(protein_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=[center['center_x'],
                                center['center_y'],
                                center['center_z']],
                        box_size=[size['size_x'],
                                  size['size_y'],
                                  size['size_z']])
    # Dock the ligand
    v.dock()
    v.write_poses(str(dockout_pdbqt), overwrite=True)

    # 将pdbqt文件转换为sdf文件
    pdbqt_to_sdf(str(dockout_pdbqt), str(dockout_sdf))
    log.info(f'Vina ends.\n')
    log.info(f'Vina dock output is {dockout_sdf}\n')


def batch_vina_dock(complex_file, ligand_dir, outdir):
    for i, ligand_file in enumerate(Path(ligand_dir).iterdir()):
        if ligand_file.suffix != ".mol2":
            continue

        vina_dock(complex_file, ligand_file, outdir)
        log.info(f'{i+1} ligands have been docked.\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Vina docking')
    parser.add_argument("-c", '--complex', type=str,
                        help='complex file', required=True)
    parser.add_argument("-l", '--ligand', type=str,
                        help='ligand file or ligand directory', required=True)
    parser.add_argument("-o", '--outdir', type=str,
                        help='output directory', required=True)
    args = parser.parse_args()

    setup_logging(f"{args.outdir}/dock.log")

    if Path(args.ligand).is_dir():
        batch_vina_dock(args.complex, args.ligand, args.outdir)
    else:
        vina_dock(args.complex, args.ligand, args.outdir)
