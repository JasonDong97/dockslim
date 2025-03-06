# 服务端代码（consumer_direct.py）
import pika


class RabbitClient:
    host: str
    port: int
    username: str
    password: str

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
        print(
            f" [x] ExchangeName: {exchange_name}, RoutingKey: {routing_key}, Sending message: {body}"
        )
        self.channel.basic_publish(
            exchange_name, routing_key, body, properties
        )
