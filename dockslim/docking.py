#!python
import subprocess
import argparse
import json
import base64
import rdkit2pdbqt
import time
import socket
import shutil

from pathlib import Path
from rdkit import Chem
from openbabel import pybel
from vina import Vina
from logs import BaseLogger, MQLogger
from os import environ as env
from rabbitmq import RabbitClient
from pika import BasicProperties

class VinaDock:
    def __init__(self, complex_pdb, ligand, outdir):
        self.complex_pdb = complex_pdb
        self.ligand = ligand
        self.outdir = outdir

        log.info(f"Complex pdb: {complex_pdb}")
        self.run_lepro(complex_pdb)

        log.info(f"Box center coordinate: {self.center}. ")
        log.info(f"Box size coordinate: {self.box_size}. \n")

    def run_lepro(self, complex_pdb: str):
        """
        将复合物PDB文件转换为lepro提供的输入文件

        :param complex_pdb: 复合物PDB文件
        :param out_pdb_file: 输出文件夹
        :return:
            center: 盒子中心坐标
            size: 盒子的大小
        """
        file_path = Path(complex_pdb).absolute()
        workdir = file_path.parent

        subprocess.run(["lepro_linux_x86", file_path], cwd=workdir)
        dock_in = workdir.joinpath("dock.in")
        self.parse_dock_in(dock_in)
        dock_in.unlink(missing_ok=True)

        # Path(out_pdb_file).parent.mkdir(parents=True, exist_ok=True)
        # workdir.joinpath("pro.pdb").rename(out_pdb_file)
        self.protein_file = workdir.joinpath("pro.pdb")

    def parse_dock_in(self, file_path):
        """
        解析 dock.in 文件，提取 Vina 对接盒子参数
        返回格式: (center_dict, size_dict)
        """
        with open(file_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]

        # 定位 Binding pocket 部分
        try:
            bp_index = lines.index("Binding pocket") + 1
        except ValueError:
            raise ValueError("Binding pocket section not found in input file")

        # 提取三轴坐标范围
        axes = []
        for i in range(3):
            if bp_index + i >= len(lines):
                raise ValueError("Incomplete binding pocket coordinates")

            parts = lines[bp_index + i].split()
            if len(parts) != 2:
                raise ValueError(f"Invalid coordinate format at line {bp_index + i}")

            axes.append((float(parts[0]), float(parts[1])))

        # 计算中心坐标和盒子尺寸
        self.center = [(axis[0] + axis[1]) / 2 for axis in axes]
        self.box_size = [axis[1] - axis[0] for axis in axes]
        # self.center = {"center_x": center[0], "center_y": center[1], "center_z": center[2]}
        # self.size = {"size_x": size[0], "size_y": size[1], "size_z": size[2]},

    def pdb_to_pdbqt(self, pdb_file: str, out_pdbqt_file: str):
        """
        将PDB文件转换为PDBQT文件
        :param pdb_file: PDB文件路径
        :param out_pdbqt_file: 输出PDBQT文件路径
        :return:
        """
        # 从PDB文件中读取拓扑结构
        receptor_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        # 将拓扑结构转换为PDBQT文件
        with open(out_pdbqt_file, "w") as f:
            # 将拓扑结构转换为PDBQT文件
            lines = rdkit2pdbqt.MolToPDBQTBlock(receptor_mol, False, False, True)
            # 将PDBQT文件写入文件
            f.write(lines)

    def pdbqt_to_sdf(self, pdbqt_file: str = None, output: str = None):
        """
        将PDBQT文件转换为SDF文件
        :param pdbqt_file: PDBQT文件路径
        :param output: SDF文件路径
        :return: SDF文件路径
        """
        # 从pdbqt文件中读取结果
        results = [m for m in pybel.readfile(filename=pdbqt_file, format="pdbqt")]
        # 创建一个输出文件
        out = pybel.Outputfile(filename=output, format="sdf", overwrite=True)
        # 遍历结果，更新pose的数据
        for pose in results:

            pose.data.update({"Pose": pose.data["MODEL"]})
            pose.data.update({"Score": pose.data["REMARK"].split()[2]})
            # 删除pose中的MODEL、REMARK、TORSDO
            del pose.data["MODEL"], pose.data["REMARK"], pose.data["TORSDO"]
            # 写入文件
            out.write(pose)
        # 关闭文件
        out.close()

    def convert(
        self,
        input_file: str,
        input_format,
        out_file: str,
        out_format: str,
        add_h=True,
    ):

        # 由配体的pdb文件,使用 pybel库 转化格式得到mol2,sdf文件，到 './protein_files' 目录下
        mol = [m for m in pybel.readfile(filename=input_file, format=input_format)][0]
        if add_h:
            mol.addh()  # 加氢
        mol.write(format=out_format, filename=out_file, overwrite=True)

    def docking(self, ligand_file):

        ligand_file_path = Path(ligand_file)
        workdir = Path(self.outdir).joinpath(ligand_file_path.stem)
        workdir.mkdir(parents=True, exist_ok=True)

        # protein_file = workdir.joinpath("protein.pdb")
        protein_pdbqt = workdir.joinpath("protein.pdbqt")
        ligand_pdbqt = workdir.joinpath("ligand.pdbqt")
        dockout_pdbqt = workdir.joinpath("dockout.pdbqt")
        dockout_sdf = workdir.joinpath("dockout.sdf")

        protein_name = Path(self.complex_pdb).name
        ligand_file_name = ligand_file_path.name

        log.info(f"{ligand_file_name} is docking with {protein_name}, Please wait...")
        start_time = time.time()

        # log.info(f"Workdir: {workdir}")

        # 将pdb文件转换为pdbqt文件
        self.pdb_to_pdbqt(self.protein_file, protein_pdbqt)
        # 将mol2文件转换为pdbqt文件
        self.convert(str(ligand_file), "mol2", str(ligand_pdbqt), "pdbqt")

        v = Vina()
        v.set_receptor(str(protein_pdbqt))
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.compute_vina_maps(center=self.center, box_size=self.box_size)
        # Dock the ligand

        v.dock()
        v.write_poses(str(dockout_pdbqt), overwrite=True)

        # 将pdbqt文件转换为sdf文件
        self.pdbqt_to_sdf(str(dockout_pdbqt), str(dockout_sdf))

        log.info(
            f"{ligand_file_name} is docked with {protein_name}, time comsuming: {time.time() - start_time} s. "
        )
        # log.info(f"Vina dock output is {dockout_sdf}")
        log.info(
            {
                "ligand_name": ligand_file_path.stem,
                "dockout_sdf": file_to_base64(dockout_sdf),
            }
        )

    def docking_batch(self):
        ligand_files = [p for p in Path(self.ligand).iterdir() if p.suffix == ".mol2"]
        for i in range(len(ligand_files)):
            ligand_file = ligand_files[i]
            log.info(f"{'-' * 50} ({i+1}/{len(ligand_files)}) {'-' * 50}")
            self.docking(ligand_file)


class DockingListener(RabbitClient):
    def __init__(self, host, port, username, password):
        super().__init__(host, port, username, password)
        self.exchange_name = "pipeline"
        self.queue_name = "pipeline.dock"
        self.routing_key = "pipeline.dock"

    def listen(self):
        self.direct_listener(
            exchange_name=self.exchange_name,
            routing_key=self.routing_key,
            queue_name=self.queue_name,
            on_message_callback=self.on_message,
        )

    def on_message(self, channel, method, props: BasicProperties, body):
        global log
        log.set_reply_to(props.reply_to)
        log.set_correlation_id(props.correlation_id)
        
        try:
            body = json.loads(body)
            log.info(f"node: {socket.gethostname()}, received message: {body}")
            if isinstance(body, dict) and body["start"]:
                self.run_demo()
        except Exception as e:
            log.info(f"node: {socket.gethostname()}, error: {e}")
            
        log.info("done")

    def run_demo(self):
        examples_path = Path(__file__).parent.parent.joinpath("examples")
        complex_pdb = examples_path.joinpath("2F0Z1.pdb")
        ligand = examples_path.joinpath("ligands")
        outdir = examples_path.joinpath("outdir")
        VinaDock(complex_pdb, ligand, outdir).docking_batch()
        shutil.rmtree(outdir)


def cli(complex_pdb, ligand, outdir):
    global log
    log = BaseLogger(f"{outdir}/docking.log")

    vina = VinaDock(complex_pdb, ligand, outdir)
    if Path(ligand).is_dir():
        vina.docking_batch()
    else:
        vina.docking(complex_pdb, ligand, outdir)


def mq_listener():
    dock_listener = DockingListener(
        host=env["RABBITMQ_HOST"],
        port=env["RABBITMQ_PORT"],
        username=env["RABBITMQ_USERNAME"],
        password=env["RABBITMQ_PASSWORD"],
    )

    global log
    log = MQLogger(rabbit_client=dock_listener)

    dock_listener.listen()


def main():
    parser = argparse.ArgumentParser(description="Vina docking")
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Run mode")

    cli_parser = subparsers.add_parser("cli", help="CLI mode")
    cli_parser.add_argument(
        "-c", dest="complex", type=str, required=True, help="complex file"
    )
    cli_parser.add_argument(
        "-l",
        dest="ligand",
        type=str,
        required=True,
        help="ligand file or ligand directory",
    )
    cli_parser.add_argument(
        "-o", dest="outdir", type=str, required=True, help="output directory"
    )

    subparsers.add_parser("mq", help="RabbitMQ mode")
    args = parser.parse_args()

    if args.mode == "mq":
        mq_listener()
    elif args.mode == "cli":
        cli(args.complex, args.ligand, args.outdir)


def file_to_base64(file_path):
    """
    将文件转换为 Base64 编码的字符串。

    :param file_path: 文件路径
    :return: Base64 编码的字符串
    """
    try:
        with open(file_path, "rb") as file:
            file_content = file.read()
            base64_encoded = base64.b64encode(file_content)
            return base64_encoded.decode("utf-8")
    except Exception as e:
        print(f"Error: {e}")
        return None


if __name__ == "__main__":
    main()
