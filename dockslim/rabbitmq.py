# 服务端代码（consumer_direct.py）
import pika
import uuid
import json


class RabbitClient:
    host: str
    port: int
    username: str
    password: str
    responses: dict = {}

    def __init__(self, host, port, username, password):
        """
        初始化RabbitMQ客户端。

        参数:e
        - host: RabbitMQ服务器的主机地址。
        - port: RabbitMQ服务器的端口号。
        - username: 用于连接RabbitMQ服务器的用户名。
        - password: 用于连接RabbitMQ服务器的密码。
        """

        self.host = host
        self.port = port

        print(" [*] RabbitMQ Connecting to {}:{}".format(self.host, self.port))
        self.channel = pika.BlockingConnection(
            pika.ConnectionParameters(
                host=host,
                port=port,
                virtual_host="/",
                credentials=pika.PlainCredentials(username, password),
            )
        ).channel()

    def direct_listener(
        self,
        exchange_name,
        routing_key,
        queue_name,
        on_message_callback,
        auto_ack=True,
        durable=True,
    ):
        self.channel.exchange_declare(exchange_name, durable=True)
        self.channel.queue_declare(queue_name, durable=durable)
        self.channel.queue_bind(queue_name, exchange_name, routing_key)
        self.channel.basic_qos(prefetch_count=1)
        self.channel.basic_consume(queue_name, on_message_callback, auto_ack)

        print(
            f" [*] Waiting for messages. excange name: {exchange_name}, routing key: {routing_key}, queue name: {queue_name}, To exit press CTRL+C"
        )
        self.channel.start_consuming()

    def send(self, exchange_name, routing_key, body, properties=None):
        self.channel.basic_publish(
            exchange_name, routing_key, json.dumps(body), properties
        )

    def send_and_wait_reply(
        self, exchange_name, routing_key, body, on_message_process_callback
    ):
        correlation_id = str(uuid.uuid4())
        callback_queue = self.channel.queue_declare("", auto_delete=True).method.queue
        self.send(
            exchange_name,
            routing_key,
            body,
            properties=pika.BasicProperties(
                reply_to=callback_queue, correlation_id=correlation_id
            ),
        )

        def reply_callback(channel, method, props: pika.BasicProperties, body: bytes):
            message = json.loads(body)

            if props.correlation_id == correlation_id:
                if isinstance(message, str):
                    print(message)

                elif isinstance(message, dict) and on_message_process_callback:
                    on_message_process_callback(message)

                if message == "done":
                    self.channel.stop_consuming(consumer_tag=correlation_id)
                    method_frame = self.channel.queue_delete(callback_queue)
                    print(f"message is done, stop consuming.")
                    print(
                        f"message is done, delete queue, method_frame: {method_frame}"
                    )

        self.channel.basic_consume(
            queue=callback_queue,
            on_message_callback=reply_callback,
            auto_ack=True,
            consumer_tag=correlation_id,
        )
        self.channel.start_consuming()
