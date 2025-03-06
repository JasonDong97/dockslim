# 服务端代码（consumer_direct.py）
import pika
import uuid
import json
import logging as log


class RabbitClient:
    host: str
    port: int
    username: str
    password: str
    responses: dict = {}

    def __init__(self, host, port, username, password):
        """
        初始化RabbitMQ客户端。

        参数:
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

        print(" [*] Waiting for messages. To exit press CTRL+C")
        self.channel.start_consuming()

    def send(self, exchange_name, routing_key, body, properties=None):
        log.info(
            f" [x] ExchangeName: {exchange_name}, RoutingKey: {routing_key}, Sending message: {body}"
        )

        self.channel.basic_publish(
            exchange_name, routing_key, json.dumps(body), properties
        )

    def send_and_receive(self, exchange_name, routing_key, body, time_limit=None):
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

        def receive_callabck(channel, method, props: pika.BasicProperties, body: bytes):
            message = json.loads(body)

            if props.correlation_id == correlation_id :
                if isinstance(message,str):
                    print(message)

                if message == "done":
                    self.channel.stop_consuming(consumer_tag=correlation_id)
                    method_frame = self.channel.queue_delete(callback_queue)
                    print(f"message is done, stop consuming.")
                    print(
                        f"message is done, delete queue, method_frame: {method_frame}"
                    )

        self.channel.basic_consume(
            queue=callback_queue,
            on_message_callback=receive_callabck,
            auto_ack=True,
            consumer_tag=correlation_id,
        )
        self.channel.start_consuming()
        # self.channel.connection.process_data_events(time_limit=time_limit)


if __name__ == "__main__":
    from os import environ as env

    client = RabbitClient(
        host=env["RABBITMQ_HOST"],
        port=env["RABBITMQ_PORT"],
        username=env["RABBITMQ_USERNAME"],
        password=env["RABBITMQ_PASSWORD"],
    )

    client.send_and_receive(
        exchange_name="pipeline",
        routing_key="pipeline.dock",
        body={"start": True},
    )
