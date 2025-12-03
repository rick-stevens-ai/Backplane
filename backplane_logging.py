#!/usr/bin/env python3
"""
Centralized logging configuration for the Scientific AI Backplane.
Provides color-coded, detailed logging for screen recording demonstrations.
"""

import logging
import sys
from typing import Optional


# ANSI color codes for terminal output
class Colors:
    """ANSI color codes for colored terminal output"""
    RESET = '\033[0m'
    BOLD = '\033[1m'

    # Regular colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    # Bright colors
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class ColoredFormatter(logging.Formatter):
    """Custom formatter with color-coded log levels"""

    # Color mapping for each log level
    LEVEL_COLORS = {
        logging.DEBUG: Colors.BRIGHT_BLACK,
        logging.INFO: Colors.CYAN,
        logging.WARNING: Colors.YELLOW,
        logging.ERROR: Colors.RED,
        logging.CRITICAL: Colors.BRIGHT_RED + Colors.BOLD
    }

    # Component colors
    COMPONENT_COLOR = Colors.MAGENTA
    MESSAGE_COLOR = Colors.WHITE

    def format(self, record):
        """Format log record with colors"""
        # Color the log level
        level_color = self.LEVEL_COLORS.get(record.levelno, Colors.WHITE)
        record.levelname = f"{level_color}{record.levelname:8s}{Colors.RESET}"

        # Color the logger name (component)
        record.name = f"{self.COMPONENT_COLOR}{record.name}{Colors.RESET}"

        # Format the message
        formatted = super().format(record)

        return formatted


def setup_logging(level: int = logging.INFO, log_file: Optional[str] = None) -> None:
    """
    Configure logging for the Backplane system.

    Args:
        level: Logging level (default: INFO)
        log_file: Optional file path to also log to a file
    """
    # Create root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove existing handlers
    root_logger.handlers = []

    # Console handler with colors
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)

    # Format: [LEVEL] [Component] Message
    console_format = ColoredFormatter(
        fmt='[%(levelname)s] [%(name)s] %(message)s',
        datefmt='%H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    root_logger.addHandler(console_handler)

    # Optional file handler (no colors in files)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_format = logging.Formatter(
            fmt='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        root_logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for a specific component.

    Args:
        name: Component name (e.g., 'agent', 'wrapper.cp2k', 'mace')

    Returns:
        Configured logger instance
    """
    return logging.getLogger(name)


# Convenience function for logging section headers
def log_section(logger: logging.Logger, title: str, char: str = '=', width: int = 80):
    """
    Log a prominent section header.

    Args:
        logger: Logger instance
        title: Section title
        char: Character to use for the border
        width: Width of the section header
    """
    border = char * width
    logger.info("")
    logger.info(border)
    logger.info(f"{title:^{width}}")
    logger.info(border)


# Convenience function for logging subsection headers
def log_subsection(logger: logging.Logger, title: str, char: str = '-', width: int = 80):
    """
    Log a subsection header.

    Args:
        logger: Logger instance
        title: Subsection title
        char: Character to use for the border
        width: Width of the subsection header
    """
    border = char * width
    logger.info("")
    logger.info(title)
    logger.info(border)


# Convenience function for logging key-value pairs
def log_dict(logger: logging.Logger, data: dict, title: Optional[str] = None):
    """
    Log a dictionary in a formatted way.

    Args:
        logger: Logger instance
        data: Dictionary to log
        title: Optional title for the data
    """
    if title:
        logger.info(f"{title}:")

    for key, value in data.items():
        logger.info(f"  {key}: {value}")


# Initialize logging with default settings
# This will be called when the module is imported
setup_logging(level=logging.INFO)
