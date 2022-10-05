from typing import *

import dataclasses
import enum


class ValidationError(Exception):
    pass


@dataclasses.dataclass
class FASTASequence:
    sequence: str
    description: Optional[str] = None
    comments: Optional[List[str]] = None


class Reader:
    def __init__(self, fasta_file: TextIO):
        self.fasta_file = fasta_file

    def __iter__(self):
        return self

    def __next__(self) -> FASTASequence:
        comments = []
        description = None

        while (line := self.fasta_file.readline()) != None:
            if line[0] == ";":
                comments.append(line[1:])
            elif line[0] == ">":
                if description != None:
                    raise ValidationError(
                        "Description already defined for this sequence."
                    )
                description = line[1:]
            else:
                return FASTASequence(line, description, comments)

        raise StopIteration
