import sys
from pathlib import Path


posix = Path(__file__).resolve().parent.parent.as_posix()
sys.path.append(posix)
from dockslim.rabbitmq import RabbitClient


def md():
    print("md")


if __name__ == "__main__":
    client = RabbitClient(host="localhost", port=5672, username="era", password="era")
    client.direct_listener(
        exchange_name="pipeline",
        routing_key="pipeline.md",
        queue_name="pipeline.md",
        on_message_callback=md,
    )
