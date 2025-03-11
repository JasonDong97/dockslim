import cmd
import sys
from pathlib import Path

sys.path.append(Path(__file__).parent.parent.as_posix())
from dockslim.rabbitmq import RabbitClient
from dockslim.logs import MQLogger
from os import environ as env

import json
import subprocess, threading


client = RabbitClient(
    host=env["RABBITMQ_HOST"],
    port=env["RABBITMQ_PORT"],
    username=env["RABBITMQ_USERNAME"],
    password=env["RABBITMQ_PASSWORD"],
)

log = MQLogger(rabbit_client=client)


def mq_listener(ch, method, props, body):
    msg = json.loads(body)
    if isinstance(msg, dict):
        log.info(f"Receive message: {msg}")
        if msg["start"]:
            log.info("Starting demo...")
            run_demo()

    log.info("done")


def run_demo():
    cwd = Path(__file__).parent.parent.joinpath("examples/MD/md0")
    p = subprocess.Popen(
        ["gmx", "mdrun", "-v", "-deffnm", "md"],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    threading.Thread(target=log_subprocess, args=(p.stdout,)).start()
    threading.Thread(target=log_subprocess, args=(p.stderr,)).start()
    p.wait()
    


def log_subprocess(stream):
    for line in iter(stream.readline, b""):
        log.info(line.decode("utf-8").strip())


if __name__ == '__main__' :
    client.direct_listener(
        exchange_name="pipeline",
        routing_key="pipeline.md",
        queue_name="pipeline.md",
        on_message_callback=mq_listener,
    )