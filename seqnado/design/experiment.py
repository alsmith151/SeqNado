from pydantic import BaseModel
from .fastq import FastqSetIP



class IPExperiment(BaseModel):
    ip: FastqSetIP
    control: FastqSetIP | None = None

    @property
    def has_control(self) -> bool:
        return self.control is not None

    @property
    def ip_set_fullname(self) -> str:
        return self.ip.sample_name

    @property
    def control_fullname(self) -> str:
        return self.control.sample_name

    @property
    def ip_performed(self) -> str:
        return self.ip.ip_or_control_name

    @property
    def control_performed(self) -> str:
        return self.control.ip_or_control_name

    @property
    def fastqs_are_paired(self) -> bool:
        ip = self.ip.is_paired
        control = self.control.is_paired if self.control else True
        return ip and control
