from rdkit import Chem
import py3Dmol
import ipywidgets as widgets
from IPython.display import display

def visualize_sdf_conformations(sdf_file):
    # 读取SDF文件
    suppl = Chem.SDMolSupplier(sdf_file)
    mols = [mol for mol in suppl if mol is not None]
    
    if not mols:
        print("未找到有效分子，请检查SDF文件。")
        return
    
    # 判断构象类型（多个分子或单个分子多构象）
    multi_mol = len(mols) > 1
    total_confs = len(mols) if multi_mol else mols[0].GetNumConformers()
    
    if total_confs == 0:
        print("分子不包含任何构象信息。")
        return
    
    # 处理单构象情况
    if total_confs == 1:
        view = py3Dmol.view(width=400, height=300)
        mol = mols[0] if multi_mol else mols[0]
        view.addModel(Chem.MolToMolBlock(mol), 'sdf')
        view.setStyle({'stick': {}})
        view.zoomTo()
        display(view.to_widget())  # 转换为 Widget
        return
    
    # 多构象处理
    view = py3Dmol.view(width=500, height=400)
    
    # 添加所有构象到视图
    for conf_id in range(total_confs):
        if multi_mol:
            mb = Chem.MolToMolBlock(mols[conf_id])
        else:
            mb = Chem.MolToMolBlock(mols[0], confId=conf_id)
        view.addModel(mb, 'sdf')
    
    view.setStyle({'stick': {}})
    view.setModel(0)  # 初始显示第一个构象
    view.zoomTo()
    
    # 创建交互控件
    slider = widgets.IntSlider(
        value=0,
        min=0,
        max=total_confs-1,
        step=1,
        description='选择构象:',
        continuous_update=False
    )
    
    # 构象更新回调函数
    def update_conformer(change):
        view.setModel(change.new)
        view.zoomTo()
        view.render()
    
    slider.observe(update_conformer, names='value')
    
    # 转换为 Widget 并显示
    view_widget = view.to_widget()
    display(widgets.VBox([slider, view_widget]))  # 使用转换后的 Widget