from typing import Dict


class Configurable:
    def __init__(self, *args):
        self.configs = list(args)

    def __repr__(self):
        return f"Configurable({', '.join(map(str, self.configs))})"


class Config:
    def __init__(self, **kwargs):
        raise NotImplementedError

    @classmethod
    def from_dict(cls, kwargs: Dict):
        updated_kwargs = {}
        for key, value in kwargs.items():
            if key in cls.__dataclass_fields__.keys():
                type_is_union = hasattr(cls.__dataclass_fields__[key].type, "__args__")
                if isinstance(value, dict):
                    if type_is_union:
                        dataclass_type = cls.__dataclass_fields__[key].type.__args__[0]
                    else:
                        dataclass_type = cls.__dataclass_fields__[key].type
                    if issubclass(
                        dataclass_type, Config
                    ):  # NOTE: Using Optional[] can break this
                        value = dataclass_type.from_dict(value)
                updated_kwargs[key] = value
        return cls(**updated_kwargs)
