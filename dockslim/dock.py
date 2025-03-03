from pathlib import Path
import rdkit2pdbqt
import subprocess

from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import pybel, openbabel as ob
import numpy as np
from vina import Vina


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
def pdb_to_pdbqt_with_pybel(pdb_file: str, pdbqt_file: str):
    '''
    将PDB文件转换为PDBQT文件
    :param pdb_file: PDB文件路径
    :param out_pdbqt_file: 输出PDBQT文件路径
    :return:
    '''
    
    # 确保输入输出路径为字符串
    pdb_file = str(pdb_file)
    pdbqt_file = str(pdbqt_file)

    # 使用低层API确保兼容性
    conv = ob.OBConversion()
    conv.SetInFormat("pdb")
    conv.SetOutFormat("pdbqt")

    # 读取分子
    mol = ob.OBMol()
    if not conv.ReadFile(mol, pdb_file):
        raise IOError(f"无法读取文件: {pdb_file}")

    # 添加氢原子
    mol.AddHydrogens()

    # 计算Gasteiger电荷
    #charge_model = ob.OBChargeModel.FindType("GASTEIGER")
    #charge_model.ComputeCharges(mol)

    # 写入文件
    if not conv.WriteFile(mol, pdbqt_file):
        raise IOError(f"无法写入文件: {pdbqt_file}")

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


def mol2_to_pdbqt(mol2_path: str, pdbqt_path: str) -> None:

    # 读取mol2文件
    mol = Chem.MolFromMol2File(mol2_path, sanitize=True, removeHs=False)
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    # 初始化PDBQT写入器
    writer = PDBQTWriterLegacy()

    # 处理每个MoleculeSetup实例
    for setup in mol_setups:
        # 写入文件
        pdbqt_data = writer.write_string(setup)

        # 如果是元组，合并为字符串（例如用换行符连接）
        if isinstance(pdbqt_data, tuple):
            pdbqt_str = pdbqt_data[0]  # 假设元组元素均为字符串
        else:
            pdbqt_str = pdbqt_data

        # 写入PDBQT文件
        with open(pdbqt_path, 'w') as f:
            f.write(pdbqt_str)


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


if __name__ == '__main__':
    # 示例用法
    center, size = get_hetatm_docking_box("tests/2F0Z1.pdb")
    print("中心坐标:", center)
    print("盒子尺寸:", size)

    protein_file = "tests/protein.pdb"
    protein_pdbqt = "tests/2F0Z1.pdbqt"
    ligand_file = "tests/D001-200081531.mol2"
    ligand_pdbqt = "tests/D001-200081531.pdbqt"
    dockout_pdbqt = "tests/D001-200081531_dockout.pdbqt"

    run_lepro("tests/2F0Z1.pdb", protein_file)

    # protein_file,_ = split_protein_ligand("tests/2F0Z1.pdb")
    
    #pdb_to_pdbqt_with_pybel(protein_file, protein_pdbqt)
    pdb_to_pdbqt(protein_file, protein_pdbqt)
    mol2_to_pdbqt(ligand_file, ligand_pdbqt)

    v = Vina()
    v.set_receptor(protein_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=[center['center_x'],
                                center['center_y'],
                                center['center_z']],
                        box_size=[size['size_x'],
                                  size['size_y'],
                                  size['size_z']])
    # Dock the ligand
    v.dock()
    v.write_poses(dockout_pdbqt, n_poses=10, overwrite=True)

    pdbqt_to_sdf(dockout_pdbqt, "tests/D001-200081531_dockout.sdf")
