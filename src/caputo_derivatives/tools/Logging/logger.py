import io
import os
from enum import Enum
from typing import Optional, TextIO

class LoggingStatus(Enum):
  fatal   = 0
  error   = 1
  warning = 2
  info    = 3
  print   = 4


class logger:
  def __init__(self, file:Optional[TextIO] = None, status:LoggingStatus = LoggingStatus.error) -> None:
    if isinstance(file, io.IOBase) and (file is not os.devnull):
      self.f = file
    else: # Backup write to devnull if gotten past log level check
      self.f = os.devnull
    self.status = status
  def display(self, *messages:str, end='\n'):
    if self.status.value > LoggingStatus.info.value:
      print(*messages, end=end)
  def info(self, *messages:str, end='\n'):
    if self.status.value > LoggingStatus.warning.value:
      self.f.write(" ".join([f'{m}' for m in messages]))
      self.f.write(end)
  def print(self, *messages:str, end='\n'):
    self.info(*messages, end=end)
    self.display(*messages, end=end)
