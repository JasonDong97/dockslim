from stabilizer.workflow.quantum_chem import mol_Ising_Ham_res

from rdkit import Chem
from pyscf import gto, scf
from pathlib import Path
import argparse
import json


def read_mol(file_path):
    mol = Chem.MolFromMolFile(file_path, removeHs=False)
    return extract_atom_coordinates(mol)


def read_mol2(file_path):
    mol = Chem.MolFromMol2File(file_path)
    # 提取分子中的原子坐标
    return extract_atom_coordinates(mol)


def extract_atom_coordinates(mol):
    atom_coords = [
        (atom.GetSymbol(), tuple(mol.GetConformer().GetAtomPosition(atom.GetIdx())))
        for atom in mol.GetAtoms()
    ]
    return atom_coords


# 输出文件路径
def calc_method(file_path):
    method = Path(file_path).suffix

    # 读取分子坐标
    if method == ".mol":
        return read_mol(file_path)
    elif method == ".mol2":
        return read_mol2(file_path)


def calc(file_path, unit, basis):
    atom_coords = calc_method(file_path)
    return calc_by_atom_coords(atom_coords, unit, basis)


def calc_by_atom_coords(atom_coords, unit="A", basis="sto-3g"):

    mol = gto.M(
        atom=atom_coords,  # 分子构型
        unit=unit,  # 单位
        basis=basis,  # 基组cc-pVDZ
        spin=None,  # 自旋
        # charge = 1 # 电荷
    )

    nao = mol.nao
    elec = mol.nelectron
    atom_str = mol.elements
    element_counts = {element: atom_str.count(element) for element in set(atom_str)}

    chemical_formula = ""
    for element, count in element_counts.items():
        if count > 1:
            chemical_formula += f"{element}{count}"
        else:
            chemical_formula += f"{element}"

    atom_str_pos = chemical_formula

    mf = scf.RHF(mol)
    mf.verbose = 0  # 取消PySCF的输出
    mf.run()

    # 获取 Hartree Fork 能量值
    ene_hf = mf.e_tot

    ene_Ising, elec_state, elec_spin = mol_Ising_Ham_res(
        mf,  # scf.hf 对象
        method="sg3d",  # 三种方法 sec, sg3r, sg3d 默认使用最后一种，效果最好
        get_state=True,  # 如果只需要能量，设置这两个选项为 False .
        get_spin=True,
    )

    tmp_dict = {}
    tmp_dict["Chemical_Formula"] = chemical_formula
    tmp_dict["Atom_Positions"] = atom_str_pos
    tmp_dict["nao"] = int(nao)
    tmp_dict["nelec"] = int(elec)
    tmp_dict["Energy_HF"] = ene_hf
    tmp_dict["Energy_Ising"] = ene_Ising
    tmp_dict["Electron_State"] = elec_state
    tmp_dict["Energy_Electron_Spin"] = elec_spin

    # with open(result_path, "w") as json_file:
    #     json.dump(tmp_dict, json_file)  # 保存结果到 json 文件
    return tmp_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Molecular potential energy calculation."
    )

    parser.add_argument("--file_path", type=str, help="input file path")
    parser.add_argument("--unit", default="A", type=str, help="calculation unit")
    parser.add_argument("--basis", default="sto3g", type=str, help="calculation basis")

    args = parser.parse_args()

    file_path = args.file_path
    unit = args.unit
    basis = args.basis
    result_path = Path(file_path).stem + ".json"

    result = calc(file_path, unit, basis)
    with open(result_path, "w") as json_file:
        json.dump(result, json_file)
