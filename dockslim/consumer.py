# 服务端代码（consumer_direct.py）
import subprocess
import pika
from pydantic import BaseModel
from pathlib import Path


class DockingConsumer:
    host: str
    port: int
    username: str
    password: str

    class Message(BaseModel):
        uuid: str
        start_dock: bool = False

    def __init__(self, host, port, username, password):
        self.host = host
        self.port = port
        cred = pika.PlainCredentials(username, password)
        params = pika.ConnectionParameters(host, port, "/", cred)
        conn = pika.BlockingConnection(params)
        self.channel = conn.channel()
        self.exchange_name = "pipeline"
        self.queue_name = "pipeline.dock"
        self.routing_key = "pipeline.dock"
        self.channel.exchange_declare(
            exchange=self.exchange_name, exchange_type="direct", durable=True
        )
        self.channel.queue_declare(queue=self.queue_name, durable=True)
        self.channel.queue_bind(
            exchange=self.exchange_name,
            queue=self.queue_name,
            routing_key=self.routing_key,
        )

        print(" [*] Connecting to {}:{}".format(self.host, self.port))
        print(" [*] Awaiting RPC requests on {}".format(self.queue_name))

    def start_consuming(self):

        self.channel.basic_qos(prefetch_count=1)
        self.channel.basic_consume(
            queue=self.queue_name,
            on_message_callback=self.on_message_callback,
            auto_ack=False,
        )
        print(" [*] Waiting for messages. To exit press CTRL+C")
        self.channel.start_consuming()

    def on_message_callback(self, ch, method, properties, body):
        try:
            body = body.decode("utf-8")
            print("接收的消息体类型为：", type(body))
            print(f"原始消息体位：{body}")
            message = self.Message.model_validate_json(body, strict=False)
            print(f"接收的消息体为：{message}")

            if message.start_dock:
                demo_shell = Path(__file__).parent.parent.joinpath("examples/run.sh")
                subprocess.run(
                    ["bash", demo_shell.absolute().as_posix()], stdout=subprocess.PIPE
                )
        except:
            print("Received non-UTF-8 message")
        ch.basic_ack(delivery_tag=method.delivery_tag)


if __name__ == "__main__":
    from os import environ as env

    consumer = DockingConsumer(
        host=env["RABBITMQ_HOST"],
        port=env["RABBITMQ_PORT"],
        username=env["RABBITMQ_USERNAME"],
        password=env["RABBITMQ_PASSWORD"],
    ).start_consuming()
