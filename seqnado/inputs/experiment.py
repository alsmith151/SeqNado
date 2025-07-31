from pydantic import BaseModel
from .fastq import FastqSetIP

class ExperimentIP(BaseModel):
    ip: FastqSetIP
    control: FastqSetIP | None = None

    @property
    def has_control(self) -> bool:
        """Check if the experiment has a control sample."""
        return self.control is not None

    @property
    def ip_set_fullname(self) -> str:
        """Get the full sample name for the IP sample including the antibody used."""
        return self.ip.full_sample_name

    @property
    def control_fullname(self) -> str:
        """Get the full sample name for the control sample including the control performed."""
        return self.control.full_sample_name if self.control else None

    @property
    def ip_performed(self) -> str:
        """Get the antibody used for the IP sample."""
        return self.ip.antibody

    @property
    def control_performed(self) -> str:
        """Get the antibody used for the control sample."""
        return self.control.antibody if self.control else None

    @property
    def fastqs_are_paired(self) -> bool:
        """Check if both IP and control samples are paired-end."""
        ip = self.ip.is_paired
        control = self.control.is_paired if self.control else True
        return ip and control
