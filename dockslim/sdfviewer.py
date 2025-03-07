from rdkit import Chem
import py3Dmol
from ipywidgets import interact, Dropdown


def view(file_path):

    # 读取所有构象并修复化合价
    supplier = Chem.SDMolSupplier(file_path, sanitize=False)
    conformers = []
    scores = []

    for mol in supplier:
        if mol:
            try:
                mol = Chem.AddHs(mol)
                conformers.append(mol)
                scores.append(mol.GetProp("Score"))
            except Exception as e:
                print(f"跳过无效分子: {str(e)}")
                continue  # 忽略无法处理的分子

    if not conformers:
        print("⚠️ 未找到有效分子！请检查SDF文件格式。")
        return

    # 创建交互式可视化函数
    def show_conformer(pose_index=0):
        try:
            viewer = py3Dmol.view(width=600, height=400)
            mol_block = Chem.MolToMolBlock(conformers[pose_index])
            viewer.addModel(mol_block, "mol")
            viewer.setStyle({"stick": {}, "sphere": {"radius": 0.3}})
            viewer.zoomTo()
            viewer.show()
            print(f"当前构象: Pose {pose_index+1}\n对接分数: {scores[pose_index]}")
        except:
            print("无法显示3D视图")

    print(f"\n ✅ 使用下拉菜单切换不同构象的 {file_path} 3D视图\n")
    # 生成下拉菜单选项
    pose_options = [
        (f"Pose {i+1} (Score: {scores[i]})", i) for i in range(len(conformers))
    ]

    # 显示交互控件
    interact(
        show_conformer,
        pose_index=Dropdown(options=pose_options, value=0, description="选择构象:"),
    )
