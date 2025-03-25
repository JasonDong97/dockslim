from pathlib import Path

import sys

import rdkit


sys.path.append(Path(__file__).parent.parent.as_posix())

from dockslim.rabbitmq import RabbitClient
from dockslim.logs import MQLogger
from os import environ as env

import json
import ecal
from rdkit import Chem

client = RabbitClient(
    host=env["RABBITMQ_HOST"],
    port=env["RABBITMQ_PORT"],
    username=env["RABBITMQ_USERNAME"],
    password=env["RABBITMQ_PASSWORD"],
)
log = MQLogger(rabbit_client=client)


def mq_listener(ch, method, props, body):
    try:
        msg = json.loads(body)
        if isinstance(msg, dict):
            log.set_reply_to(props.reply_to)
            log.set_correlation_id(props.correlation_id)
            log.info(f"Receive message: {msg}")
            if msg["start"]:
                log.info("Starting molecular potential energy calculation ...")
                if msg["type"] and msg["type"] == 1:
                    run_demo1()
                elif msg["type"] and msg["type"] == 2:
                    run_demo2()

    except Exception as e:
        log.info(f"Error: {e}")
        log.info("done")


def run_demo1():
    cwd = Path(__file__).parent.parent.joinpath("examples/ecal/h2")
    ligand_files = sorted([p for p in Path(cwd).iterdir()])
    results = []
    for i in range(len(ligand_files)):
        ligand_file = ligand_files[i]
        log.info(f"{'-' * 50} ({i+1}/{len(ligand_files)}) {'-' * 50}")
        log.info(
            f"Calculating {ligand_file.name} molecular potential energy, please wait..."
        )

        result = ecal.calc(ligand_file, "A", "sto-3g")
        log.info(f"{ligand_file.name} Result: {result}")
        results.append((ligand_file.name, result))

    log.info({"results": results})
    log.info("done")


def run_demo2():
    file_path = Path("examples/ecal/sdf/ZMR.sdf")
    mols = Chem.SDMolSupplier(file_path, sanitize=False)

    index = 0
    results = []
    for mol in mols:
        index += 1
        pose_name = f"{file_path.stem} Pose {index}"
        log.info(f"{'-' * 50} ({index}/{len(mols)}) {'-' * 50}")
        log.info(f"Calculating {pose_name} molecular potential energy, please wait...")
        atom_coordinate = ecal.extract_atom_coordinates(mol)
        result = ecal.calc_by_atom_coords(atom_coordinate)
        log.info(f"{pose_name} Result: {result}")
        results.append((pose_name, result))

    log.info({"results": results})
    log.info("done")


if __name__ == "__main__":
    client.direct_listener(
        exchange_name="pipeline",
        routing_key="pipeline.ecal",
        queue_name="pipeline.ecal",
        on_message_callback=mq_listener,
    )
