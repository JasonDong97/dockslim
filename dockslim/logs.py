from pathlib import Path
import logging
import sys

import pika

from rabbitmq import RabbitClient


class BaseLogger:
    def __init__(self, log_file: str, log_level: int = logging.INFO):
        self.logger = self.setup_logging(log_file, log_level)

    def setup_logging(self, log_file: str, log_level: int = logging.INFO):
        """
        配置全局日志设置

        参数:
            log_file (str): 日志文件路径（默认：当前目录的 app.log)
            log_level (int): 日志级别(默认:INFO)
        """
        # 创建日志目录（如果不存在）
        log_path = Path(log_file).parent
        log_path.mkdir(parents=True, exist_ok=True)

        # 定义日志格式
        formatter = logging.Formatter(
            fmt="%(asctime)s - %(levelname)-4s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        # 获取 root logger
        logger = logging.getLogger()
        logger.setLevel(log_level)

        # 避免重复添加 handler
        if logger.hasHandlers():
            logger.handlers.clear()

        # 创建文件 handler 并设置格式
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)

        # 创建控制台 handler 并设置格式
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)

        # 添加 handlers
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        return logger

    def info(self, text):
        if isinstance(text, str):
            self.logger.info(text)


class MQLogger(BaseLogger):

    def set_reply_to(self, reply_to):
        self.reply_to = reply_to

    def set_correlation_id(self, correlation_id):
        self.correlation_id = correlation_id

    def __init__(self, rabbit_client: RabbitClient):
        super().__init__("/var/logs/docking.log", logging.INFO)
        self.rabbit_client = rabbit_client

    def info(self, text):
        if isinstance(text, str):
            self.logger.info(text)

        if self.reply_to:
            self.rabbit_client.send(
                exchange_name="",
                routing_key=self.reply_to,
                body=text,
                properties=pika.BasicProperties(correlation_id=self.correlation_id),
            )
